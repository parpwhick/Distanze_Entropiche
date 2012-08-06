/**
 * @file   distance.h
 * @brief Header per la classe @ref distance
 */

#ifndef DISTANCE_H
#define	DISTANCE_H

#include <vector>
#include <algorithm>
#include "partizioni.h"

///Classe per il calcolo di molte distanze ripetutamente - permette di risparmiare l'allocazione di vettori temporanei, permette il multithreading del calcolo, esegue il minimo numero di operazioni necessarie
class distance {
private:
    ///Alloca le variabili interne nel modo opportuno
    void allocate(int n);

    //Variabili preallocate (risparmio tempo e memoria) per il calcolo delle varie distanze
    //partizioni lineari:
    ///Partizione binaria corrispondente al fattore comune
    std::vector<char> common_factor;
    ///Prima partizione binaria ridotta
    std::vector<char> reduced1;
    ///Seconda partizione binaria ridotta
    std::vector<char> reduced2;
    ///Prodotto ridotto
    std::vector<char> product_reduced;
    ///Prodotto non ridotto
    std::vector<char> binary_product;
    //partizioni generiche:
    ///Vettore dei label del prodotto, per il calcolo delle entropie
    std::vector<product_t> product;
    ///Partizione buffer per la riduzione di 1 con 2
    general_partition ridotto1;
    ///Partizione buffer per la riduzione di 2 con 1
    general_partition ridotto2;
    ///Partizione comune nel caso in cui serva
    general_partition partizione_comune;

    ///Dimensione delle partizioni di cui calcolare la distanza
    int N;
    ///Funzione che calcola solo l'effettiva distanza tra p1 e p2, senza sapere cosa sono
    void calc_distance(const general_partition &p1, const general_partition &p2);

public:
    //i risultati dei calcoli:
    ///Risultato della distanza tra due partizioni
    double dist_shan;
    ///Risultato della distanza tra due partizioni dopo la riduzione
    double dist_shan_r;
    ///Risultato della distanza topologica tra due partizioni
    double dist_top;
    ///Risultato della distanza topologica ridotta tra due partizioni
    double dist_top_r;
    ///Risultato della distanza di Hamming tra due vettori
    double dist_ham;

    ///Operatore distanza tra due partizioni lineari
    void dist(const linear_partition& p1, const linear_partition& p2);
    ///Operatore che calcola tutte le distanze, con e senza riduzione, tra due partizioni generali
    void dist(const general_partition& p1, const general_partition& p2);
    ///Sinonimo di dist(p1,p2)
    void inline operator()(const general_partition& p1, const general_partition& p2){
        dist(p1,p2);
    }
    ///Distanza di Hamming tra due vettori
    template <typename T> void hamming_distance(const T* seq1, const T* seq2);

    ///Costruttore per partizioni lunghe N
    distance(int N);
    ///Operatore di copia - crea un oggetto @c distance delle dimensioni giuste, non copia niente.
    distance(const distance &d1);

};

#endif	/* DISTANCE_H */

