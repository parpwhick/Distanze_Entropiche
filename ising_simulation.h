/* 
 * File:   ising_simulation.h
 * Author: fake
 *
 * Created on June 27, 2012, 12:45 PM
 */

#ifndef ISING_SIMULATION_H
#define	ISING_SIMULATION_H

#include "rand_mersenne.h"

class ising_simulation {
    simulation_t update_rule;

    adj_struct NN;
    int N;
    int *config;
    int *link_energies;
    int steps_per_time;
    int skip;
    int max_link_energy;
    bool running;
    double beta;
    RandMT random;
    
    void metropolis_subset(int *subset);
public:
    void metropolis_step();
    void microcanonical_step();
    ising_simulation(adj_struct NN1, simulation_t TT, int time_length=1,int initial_time_skip=0);
    ~ising_simulation(){
        if(config) delete config;
        if(link_energies) delete link_energies;        
    }
    void set_beta(double bt) { beta = bt;}
    void set_max_energy (int m){ max_link_energy = m;}
    void measure();
    void step(int steps = 1);
    void test_run(int T);
    void init_config();
    int *copy();
    
    const int *config_reference() { return config;}
};

void time_series(adj_struct adj);

#endif	/* ISING_SIMULATION_H */

