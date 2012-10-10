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
    typedef double energy_t;
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
    ///Esegue un update di tipo metropolis, con temperatura beta, per i siti elencati nel vettore @param subset
    void metropolis_subset(std::vector<int> subset, double local_beta);
    ///Pone i link attigui ai siti elencati in subset, alla temperatura data
    void thermalize_subset(vector<int> subset, double local_beta);
    ///Nr. di bordi da termalizzare (0 per non fare nulla)
    int n_borders_thermalize;

    ///Esegue uno sweep con Metropolis
    void metropolis_step();
    ///Uno sweep microcanonico
    void microcanonical_step();
    ///Uno sweep di Creutz
    void creutz_step();
    ///Costruttore per impostare i primi valori
    ising_simulation(const adj_struct & NN1);
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
    ///temperatura inversa per il sistema (o i suoi bordi)
    double avg_beta;
    ///Restituisce un puntatore read-only alla configurazione (che puo' cambiare nel frattempo!)
    const config_t *config_reference() {
        if(!config.empty())
            return config.data();
        else
            throw(std::runtime_error(std::string("Usata configurazione non inizializzata")));
    }
    const energy_t *energy_reference() {
        if(!link_energies.empty())
	    return link_energies.data();
        else
            throw(std::runtime_error(std::string("Usate energie non inizializzate")));
    }
    int energy_size(){
        return link_energies.size();
    }
    ///Vettore contenente i vettori che rappresentano i bordi
    vector<vector<int> > borders;
    ///Calcola l'energia cinetica media per sito
    double energia_cinetica();
    ///Calcola l'energia magnetica media per sito
    double energia_magnetica();
    ///Calcola la magnetizzazione media per sito
    double magnetizzazione();
    ///Restituisce un vettore con l'energia media attorno ad ogni sito
    vector<double> local_energy();


};
///Funzione globale che esegue il calcolo delle distanze e altre statistiche tra configurazioni successive
void time_series(const adj_struct & adj);
///Simulazione con variabili ed output specializzati per il problema t-J(z) e vedere la bolla
void nagaoka_run(const adj_struct & adj);

#endif	/* ISING_SIMULATION_H */

