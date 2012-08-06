/**
 * @file  ising_simulation.h
 * @brief Header per la classe @ref ising_simulation
 */

#ifndef ISING_SIMULATION_H
#define	ISING_SIMULATION_H

#include <vector>
#include <string>
#include <stdexcept>
#include "adj_handler.h"
#include "rand_mersenne.h"
using std::vector;

///Simulazione di sistemi di Ising su topologie arbitrarie
class ising_simulation {
public:
    ///tipo di variabile utilizzata, un semplice byte
    typedef char config_t;
private:
    ///update micranonico(con bordi?) o metropolis
    simulation_t update_rule;
    ///informazioni read-only sulla topologia
    const adj_struct & NN;
    ///volume della struttura (in siti)
    int N;
    ///configurazione/stato del sistema
    vector<config_t> config;
    ///energie dei link, per update microcanonico
    vector<config_t> link_energies;
    ///Nr. di sweep del sistema per intervallo di tempo
    int steps_per_time;
    ///Quanti istanti di tempo saltare all'inizio per termalizzare
    int skip;
    ///Nel caso di una distribuzione random uniforme di energie, l'energia massima
    int max_link_energy;
    ///temperatura inversa per il sistema (o i suoi bordi)
    double beta;
    ///Generatore di numeri casuali, uno per simulazione, per poter procedere parallelamente
    RandMT random;

    ///Esegue un update di tipo metropolis, con temperatura beta, per i siti elencati nel vettore @param subset
    void metropolis_subset(std::vector<int> subset);
public:
    ///Esegue uno sweep con Metropolis
    void metropolis_step();
    ///Uno sweep microcanonico
    void microcanonical_step();
    ///Costruttore per impostare i primi valori
    ising_simulation(const adj_struct & NN1, simulation_t TT, int time_length=1,int initial_time_skip=0);
    ///Imposta beta per il sistema e ricalcola le esponenziali
    void set_beta(double bt) { beta = bt;}
    ///Imposta la massima energia per i link
    void set_max_energy (int m){ max_link_energy = m;}
    ///Stampa alcune misure (non utilizzato)
    void measure();
    ///Esegue il giusto numero di sweep con l'update selezionato
    void step(int steps = 1);
    ///Funzione di test
    void test_run(int T);
    ///Per inizializzare la configurazione in modo particolare
    void init_config();
    ///Restituisce una copia della configurazione
    vector<config_t> copy();
    ///Restituisce un puntatore read-only alla configurazione (che puo' cambiare nel frattempo!)
    const config_t *config_reference() {
        if(!config.empty())
            return config.data();
        else
            throw(std::runtime_error(std::string("Usata configurazione non inizializzata")));
    }
    
    vector<int> border1;
    vector<int> border2;
    vector<int> border3;
    ///Calcola l'energia cinetica media per sito
    double energia_cinetica();
    ///Calcola l'energia magnetica media per sito
    double energia_magnetica();
    ///Calcola la magnetizzazione media per sito
    double magnetizzazione();
};
///Funzione globale che esegue il calcolo delle distanze e altre statistiche tra configurazioni successive
void time_series(const adj_struct & adj);

#endif	/* ISING_SIMULATION_H */

