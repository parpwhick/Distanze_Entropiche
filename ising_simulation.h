/* 
 * File:   ising_simulation.h
 * Author: fake
 *
 * Created on June 27, 2012, 12:45 PM
 */

#ifndef ISING_SIMULATION_H
#define	ISING_SIMULATION_H

#include "rand_mersenne.h"
#include <vector>
using std::vector;

class ising_simulation {
    
    typedef char config_t;
    
    simulation_t update_rule;

    const adj_struct & NN;
    int N;
    vector<config_t> config;
    vector<config_t> link_energies;
    int steps_per_time;
    int skip;
    int max_link_energy;
    
    bool running;
    double beta;
    RandMT random;
    
    void metropolis_subset(std::vector<int> subset);
public:
    void metropolis_step();
    void microcanonical_step();
    ising_simulation(const adj_struct & NN1, simulation_t TT, int time_length=1,int initial_time_skip=0);
    void set_beta(double bt) { beta = bt;}
    void set_max_energy (int m){ max_link_energy = m;}
    void measure();
    void step(int steps = 1);
    void test_run(int T);
    void init_config();
    vector<config_t> copy();
    
    vector<int> border1;
    vector<int> border2;
    vector<int> border3;
    //int border_size;
    const config_t *config_reference() { return config.data();}
    double energia_cinetica();
    double energia_magnetica();
    double magnetizzazione();
};

void time_series(const adj_struct & adj);

#endif	/* ISING_SIMULATION_H */

