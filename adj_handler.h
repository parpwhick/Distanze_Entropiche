/* 
 * File:   adj_handler.h
 * Author: fake
 *
 * Created on May 27, 2012, 8:31 PM
 */

#ifndef ADJ_HANDLER_H
#define	ADJ_HANDLER_H

#include <cstdio>
#include <cstdlib>
#include <memory>

#define LEAST -1

/** \brief Struttura equivalente a una matrice di adiacenza sparsa, ma con accesso rapido.
 *
 * Dato il link k-esimo, che collega i <--k--> j,
 * gli array adi[k] e adj[k] restituiscono i e j.
 * 
 * index[i] restituisce il primo link che coinvolge i, fino a index[i+1] sono
 * tutti i link del sito i-esimo, ordinati.
 * 
 * z=fetch(i) restituisce il nr. di coordinazione di 'i' e carica i link opportuni
 * nell'array vicini[0 .. z-1]
 */
class adj_struct {
public:
    ///Volume dello spazio della configurazione
    int N;
    ///Numero di link
    int n_link;
    ///Massimo numero di coordinazione
    int zmax;
    ///Indici delle righe degli elementi non-nulli di adiacenza
    std::unique_ptr<int[]> adi;
    ///Indici delle colonne degli elementi non-nulli di adiacenza
    std::unique_ptr<int[]> adj;
    ///Indica dove iniziano i link per il sito richiesto
    std::unique_ptr<int[]> index;
    ///Array per il veloce accesso ai vicini, tramite l'indice
    mutable int *vicini;

    /**\brief Trova e carica i vicini per il sito richiesto
     * Restituisce il nr. di coordinazione di 'i' e carica i link opportuni
     * nell'array vicini[0 .. z-1], se ci sono.
     *
     * @param site Sito di cui si chiedono i vicini
     * @return Il numero di coordinazione \e z
     */
    int fetch(int site) const {
        //per evitare di sforare gli array, controlliamo
        if(site>=N){
            vicini = nullptr;
            return 0;
        }
        //il numero di coordinazione e' dato dalla differenza tra indici successivi dell'indice
        int z=index[site+1]-index[site];
        //assegna il giusto indirizzo in memoria
        vicini=&adj[index[site]];
        //precarica il vettore, per un minimo aumento di performance
        __builtin_prefetch(vicini,0,0);
        return(z);
    }

    ///Costruttore per la corretta inizializzazione dei std::unique_ptr, per evitare memory leaks
    adj_struct(int* _adi = nullptr, int* _adj = nullptr, int* _index = nullptr)
            : adi(_adi), adj(_adj), index(_index) {}
};

adj_struct adiacenza_fuzzy_line(int N);
adj_struct adiacenza_simple_line(int N);
adj_struct adiacenza_square_lattice(int lato);
adj_struct adiacenza_from_file(const char *name_vec1, const char *name_vec2, int & N);
adj_struct adiacenza_sierpinski(int GEN, int &total_size);
#endif	/* ADJ_HANDLER_H */

