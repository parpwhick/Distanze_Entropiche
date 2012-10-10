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
    dopon_problem<config_t> quantum_energy;

    nagaoka_simulation(const adj_struct & NN1) : ising_simulation(NN1), quantum_energy(1, 0.01, 300, NN1) {
        quantum_energy.set_spin_array(config.data());
        quantum_energy.set_t(opts.hopping);
        quantum_energy.set_lambda(300 * opts.hopping);
        quantum_energy.set_J(opts.J);
    }

    void metropolis_wh();
    void creutz_wh();
    void step_wh(int steps=1);
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

void nagaoka_simulation::step_wh(int steps){
    if (update_rule == METROPOLIS)
        for (int i = 0; i < steps; i++)
            metropolis_wh();

    if (update_rule == CREUTZ)
        for (int i = 0; i < steps; i++) {
            creutz_wh();
            for (int b = 0; b < n_borders_thermalize; b++)
                thermalize_subset(borders[b],opts.beta[b]);

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

        double old_gs_energy = quantum_energy.last_energy;
        config[s] = -config[s];

        dH += quantum_energy.lanczos_lowest_energy() - old_gs_energy;
        if (dH > 0 && random.get_double() > exp(-2 * avg_beta * dH)) {
            //if the energy difference is unwanted, restore the previous state
            //of the spins
            config[s] = -config[s];
            //of the recorded energy
            quantum_energy.last_energy = old_gs_energy;
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

        double old_gs_energy = quantum_energy.last_energy;
        config[s] = -config[s];
        dH += quantum_energy.lanczos_lowest_energy() - old_gs_energy;

        if (dH <= 0 || link_energies[s] >= dH)
            link_energies[s] -= dH;
        else {
            config[s] = -config[s];
            quantum_energy.last_energy = old_gs_energy;
        }

    }
}
/* microcanonical piece of code

 if (check_quantum_energy) {
            double old_gs_energy = quantum_energy.last_energy;
            if (accept1) config[s1] = -config[s1];
            if (accept2) config[s2] = -config[s2];

            dH += quantum_energy.lanczos_lowest_energy() - old_gs_energy;

            if (dH <= 0 || linkenergy >= dH)
                linkenergy -= dH;
            else {
                if (accept1) config[s1] = -config[s1];
                if (accept2) config[s2] = -config[s2];
                quantum_energy.last_energy = old_gs_energy;
            }
        }
 */

template <typename T> void init_AFM(std::vector<T> &config) {
    int side = opts.lato;
    for (int col = 0; col < side; col++)
        for (int row = 0; row < side; row++)
            config[col * side + row ] = 2 * ((row + col) % 2) - 1;
}

void nagaoka_run(const adj_struct &adj) {
    auto_stats<double> mag("Magnetization");
    auto_stats<double> radius("Radius_bubble");
    auto_stats<double> beta_est("Beta_estimate");
    auto_stats<double> E_hamiltonian("Energy: H");
    auto_stats<double> E_kin("Energy: K"), E_mag("Energy: M");
    auto_stats<double> E_tot("Energy: H+M");

    nagaoka_simulation sim(adj);

    init_AFM(sim.config);
    sim.link_energies.assign(sim.link_energies.size(), 0.0);
    sim.quantum_energy.lanczos_lowest_energy(true);

    std::clock_t start = std::clock();
    double time_diff, completed_ratio;
    fprintf(stderr, "\n");

    init_AFM(sim.config);
    sim.step_wh(opts.skip);
    vector<double> avg_local_energy = sim.local_energy();
    for (int i = 1; i < opts.n_seq + 1; i++) {
        //calcolo quantita' da stampare
        E_kin = sim.energia_cinetica();
        E_mag = - sim.energia_magnetica() * opts.J * 0.25;
        E_hamiltonian = sim.quantum_energy.last_energy * opts.J;
        E_tot = E_mag + E_hamiltonian;
        mag = sim.magnetizzazione();
        radius = sqrt(mag / 3.1415);
        beta_est = 0.25 * std::log(1. + 4. * sim.link_energies.size() / E_kin);

        //simulation step
        sim.step_wh(opts.sweeps);
        //aggiornamento partizione
        write_binary_array(sim.config_reference(), adj.N, "configurations.bin", "ab");

        if (opts.verbose) {
            //calcolo energie medie per ogni iterazione
            vector<double> local_energy = sim.local_energy();
            for (int k = 0; k < adj.N; k++)
                avg_local_energy[k] += local_energy[k];
        }

        //progress bar
        fprintf(stderr, "\r");
        time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC;
        completed_ratio = (i + opts.skip + 0.0) / (opts.n_seq+opts.skip);
        fprintf(stderr, "%.1f%% done, ETA %.0fs    ",
                completed_ratio * 100, ceil(time_diff * (1 / completed_ratio - 1)));
        fflush(stderr);
    }
    fprintf(stderr, "\n");


    //stampa medie
    std::ofstream out1("medie.txt", std::ios::app);
    out1 << std::fixed << std::setprecision(4)
            << "%L,\t\tJ,\t\tbt,\t\tbtEST,\t\te_M,\t\te_H+M,\t\te_H,\t\tM,\t\tR" << endl;
    out1 << opts.lato << "\t\t"
            << opts.J << "\t\t"
            << sim.avg_beta << "\t\t"
            << beta_est.mean() << " \t\t"
            << E_mag.mean() << " \t\t"
            << E_tot.mean() << " \t\t"
            << E_hamiltonian.mean() << " \t\t"
            << mag.mean() << " \t\t"
            << radius.mean() << endl;

    //best results
    std::ofstream out2("best.txt", std::ios::app);
    out2 << std::fixed << std::setprecision(4)
            << "%L,\t\tJ,\t\tbt,\t\tbtEST,\t\te_M,\t\te_H+M,\t\te_H,\t\tM,\t\tR" << endl;
    out2 << opts.lato << "\t\t"
            << opts.J << "\t\t"
            << sim.avg_beta << "\t\t"
            << beta_est.mean() << "\t\t"
            << E_hamiltonian << "\t\t"
            << E_mag << "\t\t"
            << E_tot << "\t\t"
            << mag << "\t\t"
            << radius << endl;

    //le variabili auto_stat qui stampano le loro medie, lasciamo una riga di spazio
    printf("\n");
}
