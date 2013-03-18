/** @file nagaoka_simulation.cpp
 * @brief Tutto il necessario per calcolare time_series() di grandezze calcolate da partizioni, implementa la classe ising_simulation::
 */
#include <vector>
#include <cmath>
#include <ctime>
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
using std::vector;
using std::string;

extern options opts;

template <typename T> void init_AFM(std::vector<T> &config);
template <typename T> void make_hole(std::vector<T> &config, int radius = 2, int set_to = 1);

class nagaoka_simulation{
public:
    ///tipo di variabile utilizzata, un semplice byte
    typedef char config_t;
    typedef double energy_t;
    dopon_problem<config_t> hamiltonian;
    
    
public:
    ///update micranonico(con bordi?) o metropolis
    simulation_t update_rule;
    ///link alle informazioni read-only sulla topologia
    const adj_struct & NN;
    ///volume della struttura (in siti)
    int N;
    ///configurazione/stato del sistema
    vector<config_t> config;
    ///energie dei link, per update microcanonico
    vector<energy_t> link_energies;
    ///Nr. di sweep del sistema per intervallo di tempo
    int steps_per_time;
    ///Quanti istanti di tempo saltare all'inizio per termalizzare
    int skip;
    ///Generatore di numeri casuali, uno per simulazione, per poter procedere parallelamente
    RandMT random;
    ///Nr. di bordi da termalizzare (0 per non fare nulla)
    int n_borders_thermalize;

    ///temperatura inversa per il sistema (o i suoi bordi)
    double avg_beta;
    ///J
    double J;
    /// absolute value of J (for the energy)
    double J_abs;
    int energy_size(){
        return link_energies.size();
    }
    ///Calcola l'energia cinetica media per sito
    double energia_cinetica();
    ///Calcola l'energia magnetica media per sito
    double energia_magnetica();
    ///Calcola la magnetizzazione media per sito
    double magnetizzazione();

    nagaoka_simulation(const adj_struct & NN1) : hamiltonian(0.01, NN1), NN(NN1) {
        avg_beta = 0;
        for (size_t i = 0; i < opts.beta.size(); i++)
            avg_beta += 1 / opts.beta[i];
        avg_beta = 1 / (avg_beta / opts.beta.size());
        steps_per_time = opts.sweeps;
        skip = opts.skip;

        N = NN.N;
        update_rule = opts.dynamics;
        J = - opts.J * 0.25; //we use spin +-1 internally, but the physical ones are +- 1/2
        J_abs = std::abs(J);

        if (opts.beta.size() == 1 || NN.borders.size() <= opts.beta.size())
            n_borders_thermalize = 0;
        else
            n_borders_thermalize = opts.beta.size();

        config.resize(NN.N);
        init_AFM(config);

        if (update_rule == CREUTZ)
            link_energies.resize(NN.N);
        else if (update_rule == MICROCANONICAL)
            link_energies.resize(NN.n_link);
        for (size_t i = 0; i < link_energies.size(); i++)
            //distribuzione esponenziale inversa
            link_energies[i] = 4 * std::ceil(-J_abs / 4. / avg_beta * std::log(1 - random.get_double()) - 1);

        hamiltonian.ground_state.assign(NN.N, 1.0);
        for (int i = 0; i < NN.N; i++)
            hamiltonian.ground_state[i] = random.get_double();

        hamiltonian.set_spin_array(config.data());
        hamiltonian.set_J(J_abs); //since the spins are +-1 internally, also here J/4
        hamiltonian.set_L(opts.lato);
        hamiltonian.set_V(0);
    }

    void step_wh(int steps=1);

    //measurements
    double average_position();
    vector<double> print_temperature_profile(double epsilon = 0.005);
    
    bool is_electron_near(int k);

private:
    void metropolis_wh();
    void metropolis_gradient();
    void creutz_wh();
    void microcanonical_wh();
};

template <typename data_t> void write_binary_array(const data_t *array, int N, std::string filename, const char * mode){
    FILE *out;

    out = fopen(filename.c_str(), mode);
    if (out == 0) {
        fprintf(stderr, "Error opening file %s for writing\n", filename.c_str());
        exit(1);
    }

    fwrite(array, sizeof (data_t), N, out);
    fclose(out);
}

double nagaoka_simulation::energia_cinetica() {
    double totale = 0.0;
    if (update_rule == CREUTZ || update_rule == MICROCANONICAL) {
        for (size_t i = 0; i < link_energies.size(); i++)
            totale += static_cast<double>(link_energies[i]);
        return totale;
    }
    else
        return 4.0/(exp(4*avg_beta)-1);
}

double nagaoka_simulation::magnetizzazione(){
    int totale = 0;
    for (int i = 0; i < N; i++)
        totale += config[i];
    double media = totale;
    //media /= N;
    return (std::abs(media));
}

double nagaoka_simulation::energia_magnetica() {
    double dH=0;
    for (int i = 0; i < NN.n_link; i++)
        dH += - J * config[NN.positive_links[i].first] * config[NN.positive_links[i].second];
    
    return dH;
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
                for (size_t j = 0; j < NN.borders[b].size(); j++)
                        link_energies[NN.borders[b][j]] = -J_abs / opts.beta[b] * std::log(random.get_double());
        }

    if(update_rule == MICROCANONICAL)
        for (int i = 0; i < steps; i++) {
            microcanonical_wh();

            //heat up the links adjacent to the borders
            for (int b = 0; b < n_borders_thermalize; b++)
                for (size_t j = 0; j < NN.bordering_links[b].size(); j++)
                         link_energies[NN.bordering_links[b][j]] = -J_abs / opts.beta[b] * std::log(random.get_double());
                                 //4 * std::ceil(-J_abs / 4. / opts.beta[b] * std::log(1 - random.get_double()) - 1);

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
        dH = 2 * J * config[s] * somma_vicini;

        double old_gs_energy = hamiltonian.last_energy;
        config[s] = -config[s];

        if(is_electron_near(s))
                dH += hamiltonian.lanczos_lowest_energy() - old_gs_energy;
        if (dH > 0 && random.get_double() > exp(- opts.beta[0] / J_abs * dH)) {
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
        dH = 2 * J * config[s] * somma_vicini;

        double old_gs_energy = hamiltonian.last_energy;
        config[s] = -config[s];

        if(is_electron_near(s))
                dH += hamiltonian.lanczos_lowest_energy() - old_gs_energy;
        
        if (dH > 0 && random.get_double() > exp(- local_beta / J_abs * dH)) {
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

        if(is_electron_near(s1) || is_electron_near(s2))
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
        double rho = square(hamiltonian.ground_state[i]);
        position+= (i/side) * rho;
    }
    return position;
}

bool nagaoka_simulation::is_electron_near(int k) {
    //return false;
    int z;
    z = NN.fetch(k);
    double total_rho = 0.0;

    total_rho += square(hamiltonian.ground_state[k]);
    for (int m = 0; m < z; m++)
        total_rho += square(hamiltonian.ground_state[NN.vicini[m]]);

    return (total_rho > 1e-5);
}

template <typename T> void shift_lattice(vector<T> &conf, int by, int L){
    static vector<T> temp;
    if(by==0)
        return;
    temp.reserve(conf.size());
    std::copy(conf.begin(),conf.end(),temp.begin());    
    
    //we map i -> i + by * L
    for(int from = 0; from < L; from++){
        //destination column
        int to = (from + by + 2*L ) % L;
        for(int i=0; i<L; i++)
            conf[to*L+i] = temp[from*L+i];
    }
}

vector<double> nagaoka_simulation::print_temperature_profile(double epsilon) {
    vector<double> avg_local_energy(opts.lato);
    vector<int> accepted(opts.lato);

    if(hamiltonian.ground_state.empty())
        return avg_local_energy;

    accepted.assign(opts.lato, 0);

    if (update_rule == MICROCANONICAL) {
        for (int i = 0; i < NN.n_link; i++) {
            int s1 = NN.positive_links[i].first;
            int s2 = NN.positive_links[i].second;
            int pos1 = (s1 / opts.lato);
            int pos2 = (s2 / opts.lato);
            double rho1 = square(hamiltonian.ground_state[s1]);
            double rho2 = square(hamiltonian.ground_state[s2]);

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

template <typename T> void make_hole(std::vector<T> &config, int radius, int set_to) {
    int side = opts.lato;
    for (int col = 0; col < side; col++)
        for (int row = 0; row < side; row++)
            if(square(col-side/2) + square(row-side/2) <= square(radius))
                config[col * side + row ] = set_to;
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
    auto_stats<int> shifts("Polaron shifts"), displacement("Shift per time");
    auto_stats<double> dist_travelled("Instant speed");
    nagaoka_simulation sim(adj);
    const double L = opts.lato;

    init_AFM(sim.config);
    double afm_gs_energy = - sim.energia_magnetica();    
    //sim.link_energies.assign(sim.link_energies.size(), 0.0);

    std::clock_t start = std::clock();
    double time_diff, completed_ratio;
    fprintf(stderr, "\n");

    //we only shift the configuration, so only metropolis is allowed
    if(!sim.update_rule == METROPOLIS){
        fprintf(stderr,"Only supporting Metropolis dynamics at the moment, sorry\n");
        //exit(2);
    }
    
    V = opts.V;
    sim.hamiltonian.set_confining(L/4);
    //sim.step(opts.skip * 100);
    sim.step_wh(opts.skip);
    sim.hamiltonian.set_confining(0);
    sim.hamiltonian.set_V(V);
    shifts = 0;


    std::ofstream out0(string("output" + opts.suffix_out + ".txt").c_str(), std::ios::out);
    out0 << std::fixed << std::setprecision(4)
            << "%J,\t\tR,\t\te_H+M+K,\t\te_bub,\t\te_H+M,\t\tM,\t\tbt,\t\tbtEST" << endl;
    for (int i = 1; i < opts.n_seq + 1; i++) {
        //calcolo quantita' da stampare
        sim.hamiltonian.lanczos_lowest_energy();
        double new_position = sim.average_position();
        dist_travelled = new_position - position;
      
        if ((new_position > (0.75 * L)) || (new_position < (0.25 * L))) {
            displacement = static_cast<int>(round((L/2) - new_position));
            shift_lattice(sim.config, displacement, L);
            shift_lattice(sim.hamiltonian.ground_state, displacement, L);
            sim.hamiltonian.lanczos_lowest_energy();
            new_position = sim.average_position();
        } else 
            displacement = 0;
        
        shifts = shifts - displacement;
        position = new_position;
        E_kin = sim.energia_cinetica();
        E_mag = - sim.energia_magnetica();
        E_hamiltonian = sim.hamiltonian.last_energy;
        E_bubble = E_hamiltonian + E_mag - afm_gs_energy;
        E_tot = E_mag + E_hamiltonian;
        E_micro = E_tot + E_kin;

        mag = sim.magnetizzazione();
        radius = sqrt(mag / 3.1415);
        beta_est = sim.link_energies.size() / (E_kin / opts.J);//0.25 * std::log(1. + 4. * sim.link_energies.size() / (E_kin / opts.J));

        //simulation step
        sim.step_wh(opts.sweeps);
        if (opts.verbose > 1) {
            write_binary_array(sim.config.data(), adj.N, ("states" + opts.suffix_out + ".bin"), "ab");
            if (opts.verbose > 2) {
                write_binary_array(sim.hamiltonian.ground_state.data(), sim.N, "electron_eigenstate" + opts.suffix_out + ".bin", "ab");
                if (opts.dynamics != METROPOLIS){
                    //average temperatures outside the bubble
                    write_binary_array(sim.print_temperature_profile().data(), opts.lato, ("avg_energies" + opts.suffix_out + ".bin"), "ab");
                    //average temperatures everywhere
                    write_binary_array(sim.print_temperature_profile(1).data(), opts.lato, ("raw_avg_energies" + opts.suffix_out + ".bin"), "ab");
                    write_binary_array(sim.link_energies.data(), sim.energy_size(), ("link_energies" + opts.suffix_out + ".bin"), "ab");         
                }
            }
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
            << V << "\t\t"
            << E_hamiltonian  << "\t\t" 
            << shifts << "\t\t"
            << displacement << "\t\t"
            << dist_travelled << endl;

        //progress bar
        fprintf(stderr, "\r");
        time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC;
        completed_ratio = (i + opts.skip + 0.0) / (opts.n_seq+opts.skip);
        fprintf(stderr, "%.1f%% done, ETA %.0fs    ",
                completed_ratio * 100, ceil(time_diff * (1 / completed_ratio - 1)));
        fflush(stderr);
    }
    time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    fprintf(stderr, "\rDone in %.1fs CPU time\n",time_diff);


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
            << beta_est.mean() << "\t\t"
            << position.mean() << "\t\t"
            << V.mean() << "\t\t"
            << shifts  << "\t\t"
            << -displacement.mean() << endl;

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
            << beta_est << "\t\t"
            << position.mean() << "\t\t"
            << V.mean() << "\t\t"
            << shifts  << "\t\t"
            << -displacement.mean() << endl;

    //le variabili auto_stat qui stampano le loro medie, lasciamo una riga di spazio
    printf("\n");
}

int main(int argc, char** argv) {
    set_program_options(opts, argc, argv);

    adj_struct adj;
    switch (opts.topologia) {
        case(RETICOLO_2D):
            adj = adiacenza_square_lattice(opts.lato);
            break;        
        case(TORO_2D):
        default:
            adj = adiacenza_toroidal_lattice(opts.lato);
            break;
    }
    
    nagaoka_run(adj);
}

/*
nagaoka_simulation sim(adj);
    const double L = opts.lato;
    //last snapshot!
    vector<nagaoka_simulation::config_t> config_backup;
    vector<nagaoka_simulation::energy_t> energy_backup;
    vector<double> groundstate_backup;
    double known_position = 0;
    displacement = 0;
    shifts = 0;
    

    init_AFM(sim.config);
    double afm_gs_energy = - sim.energia_magnetica() * opts.J * 0.25;    
    //sim.link_energies.assign(sim.link_energies.size(), 0.0);

    std::clock_t start = std::clock();
    double time_diff, completed_ratio;
    fprintf(stderr, "\n");

    V = opts.V;
    sim.step(opts.skip * 100);
    if(opts.dynamics == METROPOLIS){
        //make_hole(sim.config);
        sim.hamiltonian.set_confining(0.1*L);
    }
    sim.step_wh(opts.skip);
    sim.hamiltonian.set_confining(0);
    sim.step_wh(opts.skip);
    sim.hamiltonian.set_V(V);

    if (opts.dynamics == METROPOLIS) {
        config_backup = sim.config;
        energy_backup = sim.link_energies;
        groundstate_backup = sim.hamiltonian.ground_state;
        known_position = sim.average_position();
    }

    std::ofstream out0(("output" + opts.suffix_out + ".txt").c_str(), std::ios::out);
    out0 << std::fixed << std::setprecision(4)
            << "%J,\t\tR,\t\te_H+M+K,\t\te_bub,\t\te_H+M,\t\tM,\t\tbt,\t\tbtEST" << endl;
    for (int i = 1; i < opts.n_seq + 1; i++) {
        //calcolo quantita' da stampare
        sim.hamiltonian.lanczos_lowest_energy();
        position = sim.average_position();

        if (!config_backup.empty() && (position > 0.8 * L || position < 0.2 * L)) {
            //recenter
            //sim.hamiltonian.set_confining(2* L /4);
            //sim.step_wh(5);
            //sim.hamiltonian.set_confining(0);
            if (position > 0.5 * L)
                shifts = shifts + 1;         
            else
                shifts = shifts - 1;
            displacement = displacement + (position - known_position);
            sim.config = config_backup;
            sim.link_energies = energy_backup;
            sim.hamiltonian.ground_state = groundstate_backup;
            sim.hamiltonian.set_V(V);
        }
        if (opts.dynamics != METROPOLIS && abs(position - 0.5 * L) < 0.3) {
            config_backup = sim.config;
            energy_backup = sim.link_energies;
            groundstate_backup = sim.hamiltonian.ground_state;
            known_position = position;
        }
        */