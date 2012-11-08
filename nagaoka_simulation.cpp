/** @file ising_simulation.cpp
 * @brief Tutto il necessario per calcolare time_series() di grandezze calcolate da partizioni, implementa la classe ising_simulation::
 */
#include <vector>
#include <cmath>
#include "adj_handler.h"
#include "strutture.h"
#include "rand_mersenne.h"
#include "ising_simulation.h"
#include "distance.h"
#include "smart_data_types.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "dopon_problem.h"
using std::cout;
using std::endl;

extern options opts;

const int J = -1;

template <typename data_t> void write_binary_array(const data_t *array, int N, const char *filename, const char * mode = "wb");

using std::vector;


/**
 * @brief Scrive (appende) l'array nel file specificato
 * @param array Array da scrivere
 * @param N Numero di elementi
 * @param filename Nome del file
 */
class nagaoka_simulation : public ising_simulation {
public:
    dopon_problem<config_t> hamiltonian;

    nagaoka_simulation(const adj_struct & NN1) : ising_simulation(NN1), hamiltonian(1, 0.01, 300, NN1) {
        hamiltonian.set_spin_array(config.data());
        hamiltonian.set_t(opts.hopping);
        hamiltonian.set_lambda(300 * opts.hopping);
        hamiltonian.set_J(opts.J);
        hamiltonian.set_L(opts.lato);
        hamiltonian.set_V(0);
        setup_bordering_links();
    }

    void step_wh(int steps=1);

    //measurements
    double average_position();
    vector<double> print_temperature_profile();

    //quantity
    vector<vector<int> > bordering_links;

private:
    void metropolis_wh();
    void metropolis_gradient();
    void creutz_wh();
    void microcanonical_wh();

    void setup_bordering_links();
};

template <typename data_t> void write_binary_array(const data_t *array, int N, const char *filename, const char * mode){
    FILE *out;

    if(filename==0)
        return;

    out = fopen(filename, mode);
    if (out == 0) {
        fprintf(stderr, "Error opening file %s for writing\n", filename);
        exit(1);
    }

    fwrite(array, sizeof (data_t), N, out);
    fclose(out);
}

void nagaoka_simulation::setup_bordering_links(){
    bordering_links.resize(borders.size());

    for (int link = 0; link < NN.n_link; link++) {
        int s1 = NN.positive_links[link].first;
        int s2 = NN.positive_links[link].second;

        /* easy method for squares
         * if(s1 < L || s2 < L)
             bordering_links[0].push_back(link);
           else if( s1 >= L*(L-1) || s2 >= L*(L-1))
             bordering_links[1].push_back(link);
         */
        for (size_t b = 0; b < borders.size(); b++)
            if (std::binary_search(borders[b].begin(), borders[b].end(), s1) ||
                    std::binary_search(borders[b].begin(), borders[b].end(), s2))
                bordering_links[b].push_back(link);
    }
}

void nagaoka_simulation::step_wh(int steps){
    if (update_rule == METROPOLIS)
        for (int i = 0; i < steps; i++)
            if(opts.beta.size()==1)
                metropolis_wh();
            else
                metropolis_gradient();

    if (update_rule == CREUTZ)
        for (int i = 0; i < steps; i++) {
            creutz_wh();

            //heat up the sites making up the borders
            for (int b = 0; b < n_borders_thermalize; b++)
                for (size_t j = 0; j < borders[b].size(); j++)
                        link_energies[borders[b][j]] = -1 / opts.beta[b] * std::log(random.get_double());
        }

    if(update_rule == MICROCANONICAL)
        for (int i = 0; i < steps; i++) {
            microcanonical_wh();

            //heat up the links adjacent to the borders
            for (int b = 0; b < n_borders_thermalize; b++)
                for (size_t j = 0; j < bordering_links[b].size(); j++)
                         link_energies[bordering_links[b][j]] = -1 / opts.beta[b] * std::log(random.get_double());

        }
}

void nagaoka_simulation::metropolis_wh() {
    int s, z, somma_vicini;
    double dH;

    /* Dinamica di Metropolis, 1 passo temporale */
    for (int j = 0; j < NN.N; j++) {
        s = random.get_int() % NN.N;

        z = NN.fetch(s);
        somma_vicini = 0;
        for (int m = 0; m < z; m++)
            somma_vicini += config[NN.vicini[m]];
        dH = J * config[s] * somma_vicini;

        double old_gs_energy = hamiltonian.last_energy;
        config[s] = -config[s];

        dH += hamiltonian.lanczos_lowest_energy() - old_gs_energy;
        //dH += quantum_energy.calculate_lowest_energy() - old_gs_energy;
        if (dH > 0 && random.get_double() > exp(-2 * opts.beta[0] * dH)) {
            //if the energy difference is unwanted, restore the previous state
            //of the spins
            config[s] = -config[s];
            //of the recorded energy
            hamiltonian.last_energy = old_gs_energy;
        }

    }
}

void nagaoka_simulation::metropolis_gradient() {
    int s, z, somma_vicini;
    double dH;

    /* Dinamica di Metropolis, 1 passo temporale */
    for (int j = 0; j < NN.N; j++) {
        s = random.get_int() % NN.N;

        //position along the gradient expressed in the interval [0..1]
        double position = ((s / opts.lato) + 0.0) / (opts.lato-1);
        //local T
        double T = 1/opts.beta[0] + position * (1/opts.beta[1] - 1/opts.beta[0]);
        //local beta
        double local_beta = 1/T;

        z = NN.fetch(s);
        somma_vicini = 0;
        for (int m = 0; m < z; m++)
            somma_vicini += config[NN.vicini[m]];
        dH = J * config[s] * somma_vicini;

        double old_gs_energy = hamiltonian.last_energy;
        config[s] = -config[s];

        dH += hamiltonian.lanczos_lowest_energy() - old_gs_energy;
        if (dH > 0 && random.get_double() > exp(-2 * local_beta * dH)) {
            //if the energy difference is unwanted, restore the previous state
            //of the spins
            config[s] = -config[s];
            //of the recorded energy
            hamiltonian.last_energy = old_gs_energy;
        }
    }
}

void nagaoka_simulation::creutz_wh() {
    int s, z, sum_neigh;
    double dH;
    for (int j = 0; j < NN.N; j++) {
        s = random.get_int() % NN.N;

        z = NN.fetch(s);
        sum_neigh = 0;
        for (int m = 0; m < z; m++)
            sum_neigh += config[NN.vicini[m]];
        dH = J * 2 * config[s] * sum_neigh;

        double old_gs_energy = hamiltonian.last_energy;
        config[s] = -config[s];
        dH += hamiltonian.lanczos_lowest_energy() - old_gs_energy;

        if (dH <= 0 || link_energies[s] >= dH)
            link_energies[s] -= dH;
        else {
            config[s] = -config[s];
            hamiltonian.last_energy = old_gs_energy;
        }

    }
}

void nagaoka_simulation::microcanonical_wh() {
    double dH;
    int somma_vicini;
    int z;

    int s1,s2;
    /* Dinamica Microcanonica, 1 passo temporale */
    for (int j = 0; j < NN.N;) {
        // il generatore di numeri casuali influisce per un 5% sulla performance
        // totale, in realta' è l'accesso disordinato alla memoria che uccide

        uint32_t randnum = random.get_int();
        int link = randnum % NN.n_link;
        int flip1 = (randnum & (1 << 29)) != 0;
        int flip2 = (randnum & (1 << 30)) != 0;
        int flip12 = flip1 && flip2;

        s1 = NN.positive_links[link].first;
        s2 = NN.positive_links[link].second;

        energy_t &linkenergy = link_energies[link];

        if (!(flip1 || flip2) )
            continue;
	j++;

        if (flip1) {
            // z numero dei vicini di s1
            z = NN.fetch(s1);
            somma_vicini = 0;
            for (int m = 0; m < z; m++)
                somma_vicini += config[NN.vicini[m]];
            dH = J * 2 * config[s1] * somma_vicini;
        } else
            dH = 0;

        //energia dai vicini di s2, se cambia
        if (flip2) {
            z = NN.fetch(s2);
            somma_vicini = 0;
            for (int m = 0; m < z; m++)
                somma_vicini += config[NN.vicini[m]];
            dH += J * 2 * config[s2] * somma_vicini;
        }

        //energia dal flip simultaneo di s1 e s2
        if (flip12)
            dH -= J * 4 * config[s1] * config[s2];

        double old_gs_energy = hamiltonian.last_energy;
        if (flip1) config[s1] = -config[s1];
        if (flip2) config[s2] = -config[s2];

        dH += hamiltonian.lanczos_lowest_energy() - old_gs_energy;

        if (dH <= 0 || linkenergy >= dH)
            linkenergy -= dH;
        else {
            if (flip1) config[s1] = -config[s1];
            if (flip2) config[s2] = -config[s2];
            hamiltonian.last_energy = old_gs_energy;
        }
    }
}

double nagaoka_simulation::average_position(){
    int side = opts.lato;
    double position=0.0;
    if(hamiltonian.ground_state.empty())
        return NAN;

    for(int i=0; i < N; i++){
        double rho = hamiltonian.ground_state[i]*hamiltonian.ground_state[i];
        position+= (i/side) * rho;
    }
    return position;
}

vector<double> nagaoka_simulation::print_temperature_profile() {
    vector<double> avg_local_energy(opts.lato);
    vector<int> accepted(opts.lato);

    if(hamiltonian.ground_state.empty())
        return avg_local_energy;

    accepted.assign(opts.lato, 0);
    double epsilon = 0.005;

    if (update_rule == MICROCANONICAL) {
        for (int i = 0; i < NN.n_link; i++) {
            int s1 = NN.positive_links[i].first;
            int s2 = NN.positive_links[i].second;
            int pos1 = (s1 / opts.lato);
            int pos2 = (s2 / opts.lato);
            double rho1 = hamiltonian.ground_state[s1] * hamiltonian.ground_state[s1];
            double rho2 = hamiltonian.ground_state[s2] * hamiltonian.ground_state[s2];

            if (pos1 == pos2 && rho1 < epsilon && rho2 < epsilon) {
                accepted[pos1]++;
                avg_local_energy[pos1] += link_energies[i];
            }

        }
        for (int i = 0; i < opts.lato; i++)
            avg_local_energy[i] /= accepted[i];

    }
    return avg_local_energy;
}

template <typename T> int signof(T number){
    return 2*(number >= 0)-1;
}

template <typename T> void init_AFM(std::vector<T> &config) {
    int side = opts.lato;
    for (int col = 0; col < side; col++)
        for (int row = 0; row < side; row++)
            config[col * side + row ] = 2 * ((row + col) % 2) - 1;
}

template <typename T> void make_hole(std::vector<T> &config) {
    int side = opts.lato;
    for (int col = side / 2 - 2; col < side / 2 + 2; col++)
        for (int row = side / 2 - 2; row < side / 2 + 2; row++)
            config[col * side + row ] = 1;
}

void nagaoka_run(const adj_struct &adj) {
    auto_stats<double> mag("Magnetization");
    auto_stats<double> radius("Radius_bubble");
    auto_stats<double> beta_est("Beta_estimate");
    auto_stats<double> E_hamiltonian("Energy: H");
    auto_stats<double> E_bubble("Energy: B");
    auto_stats<double> E_kin("Energy: K"), E_mag("Energy: M");
    auto_stats<double> E_tot("Energy: H+M");
    auto_stats<double> E_micro("Energy: H+M+K");
    auto_stats<double> position("Position");
    auto_stats<double> V("Voltage");
    double voltage_step = 0.050;
    int reached_border = 0;
    nagaoka_simulation sim(adj);
    const double L = opts.lato;

    init_AFM(sim.config);
    double afm_gs_energy = - sim.energia_magnetica() * opts.J * 0.25;
    sim.link_energies.assign(sim.link_energies.size(), 0.0);

    std::clock_t start = std::clock();
    double time_diff, completed_ratio;
    fprintf(stderr, "\n");

    //init the hole confined to about L/4
    sim.hamiltonian.set_confining(2* L /4);
    sim.step_wh(opts.skip);
    sim.hamiltonian.set_confining(0);

    sim.hamiltonian.set_V(opts.V);
    V = opts.V;
    if (opts.verbose)
            write_binary_array(sim.config_reference(), adj.N, opts.config_out.c_str(), "ab");
    vector<double> avg_local_energy = sim.local_energy();
    std::ofstream out0("output.txt", std::ios::app);
    out0 << std::fixed << std::setprecision(4)
            << "%J,\t\tR,\t\te_H+M+K,\t\te_bub,\t\te_H+M,\t\tM,\t\tbt,\t\tbtEST" << endl;
    for (int i = 1; i < opts.n_seq + 1; i++) {
        //drive the voltage periodically, with T=30
        //const double 2pi = 6.283; sim.quantum_energy.set_V(sin(2pi* i/30));

        //calcolo quantita' da stampare
        position = sim.average_position();

        if(position > 0.85*L || position < 0.15*L){
            //recenter
            sim.hamiltonian.set_confining(2* L /4);
            sim.step_wh(3);
            sim.hamiltonian.set_confining(0);

            if(position > 0.5 * L)
                reached_border += 1;
            else
                reached_border -= 1;
        }
        if(abs(reached_border)>= 2){
            //change the voltage
            if(reached_border > 0){
                printf("\nLowering to %f, by %f\n",V - voltage_step,voltage_step);
                V = V - voltage_step;
            }
            else{
                printf("\nRaising to %f, by %f\n",voltage_step + V,voltage_step);
                V = V + voltage_step;
            }
            voltage_step /= 2;
            sim.hamiltonian.set_V(V);
            reached_border = 0;
        }

        E_kin = sim.energia_cinetica() * opts.J;
        E_mag = - sim.energia_magnetica() * opts.J * 0.25;
        E_hamiltonian = sim.hamiltonian.last_energy * opts.J;
        E_bubble = E_hamiltonian + E_mag - afm_gs_energy;
        E_tot = E_mag + E_hamiltonian;
        E_micro = E_tot + E_kin;

        mag = sim.magnetizzazione();
        radius = sqrt(mag / 3.1415);
        beta_est = sim.link_energies.size() / (E_kin / opts.J);//0.25 * std::log(1. + 4. * sim.link_energies.size() / E_kin);


        //simulation step
        sim.step_wh(opts.sweeps);
        if (opts.verbose) {
            write_binary_array(sim.config_reference(), adj.N, opts.config_out.c_str(), "ab");
            if (opts.simulation_type != METROPOLIS)
                write_binary_array(sim.print_temperature_profile().data(), opts.lato, "avg_energies.bin", "ab");
//            //calcolo energie medie per ogni iterazione
//            vector<double> local_energy = sim.local_energy();
//            for (int k = 0; k < adj.N; k++)
//                avg_local_energy[k] += local_energy[k];
        }
        if (opts.verbose)
            out0 << opts.J << "\t\t"
            << radius << "\t\t"
            << E_micro << "\t\t"
            << E_bubble << "\t\t"
            << E_tot << "\t\t"
            << mag << "\t\t"
            << sim.avg_beta << "\t\t"
            << beta_est << "\t\t"
            << position << "\t\t"
            << V * 1000 << "\t\t"
            << E_hamiltonian  << endl;

        //progress bar
        fprintf(stderr, "\r");
        time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC;
        completed_ratio = (i + opts.skip + 0.0) / (opts.n_seq+opts.skip);
        fprintf(stderr, "%.1f%% done, ETA %.0fs    ",
                completed_ratio * 100, ceil(time_diff * (1 / completed_ratio - 1)));
        fflush(stderr);
    }
    fprintf(stderr, "\n");
    if (opts.simulation_type != METROPOLIS)
            write_binary_array(sim.energy_reference(), sim.energy_size(), "energies_end.bin", "wb");
    write_binary_array(sim.hamiltonian.ground_state.data(),sim.N, "electron_eigenstate.bin","wb");


    //stampa medie
    std::ofstream out1("averages.txt", std::ios::app);
    out1 << std::fixed << std::setprecision(4);
    out1 << opts.J << "\t\t"
            << radius.mean() << "\t\t"
            << E_micro.mean() << "\t\t"
            << E_bubble.mean() << "\t\t"
            << E_tot.mean() << "\t\t"
            << mag.mean() << "\t\t"
            << sim.avg_beta << "\t\t"
            << beta_est.mean() << endl;

    //best results
    std::ofstream out2("best.txt", std::ios::app);
    out2 << std::fixed << std::setprecision(4);
    out2 << opts.J << "\t\t"
            << radius.max() << "\t\t"
            << E_micro.min() << "\t\t"
            << E_bubble.min() << "\t\t"
            << E_tot.min() << "\t\t"
            << mag.max() << "\t\t"
            << sim.avg_beta << "\t\t"
            << beta_est << endl;

    //le variabili auto_stat qui stampano le loro medie, lasciamo una riga di spazio
    printf("\n");
}
