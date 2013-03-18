/**
 * @file   adj_handler.h
 * @brief Definisce la classe adj_struct:: e funzioni utili per creare tali strutture nei casi noti
 */

#ifndef ADJ_HANDLER_H
#define	ADJ_HANDLER_H

#include <cstdio>
#include <cstdlib>
#include <vector>

using std::pair;
using std::vector;

#define LEAST -1

/** @brief Struttura equivalente a una matrice di adiacenza sparsa, ma con accesso rapido.
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
    ///Numero di link (n. elementi nonnulli della matrice di adiacenza)
    int n_total_links;
    ///Numero link con i>j
    int n_link;
    ///Massimo numero di coordinazione
    int zmax;
    ///Indici delle righe degli elementi non-nulli di adiacenza
    std::vector<int> adi;
    ///Indici delle colonne degli elementi non-nulli di adiacenza
    std::vector<int> adj;
    ///Indica dove iniziano i link per il sito richiesto
    std::vector<int> index;
    ///Coppie di elementi adiacenti, i~j, con i>j
    std::vector<pair<int,int> > positive_links;
    ///Array di bordi, dove i bordi sono array con i siti per ogni bordo
    vector<vector<int> > borders;
    //Quali link sono connessi ad ogni bordo?
    vector<vector<int> > bordering_links;
    ///Array per il veloce accesso ai vicini, tramite l'indice
    mutable const int *vicini;

    /**@brief Trova e carica i vicini per il sito richiesto
     * Restituisce il nr. di coordinazione di 'i' e carica i link opportuni
     * nell'array vicini[0 .. z-1], se ci sono.
     *
     * @param site Sito di cui si chiedono i vicini
     * @return Il numero di coordinazione @a z
     */
    int fetch(int site) const {
        //per evitare di sforare gli array, controlliamo
        if(site>=N){
            vicini = 0;
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

    ///Determina link con i<j
    void normalize();
    ///Determina quali link fanno parte del bordo
    void setup_bordering_links();
    ///Costruttore per allocare direttamente i vettori
    adj_struct(std::vector<int> _adi, std::vector<int> _adj , std::vector<int> _index)
            : adi(_adi), adj(_adj), index(_index) {}
    adj_struct() {}
};

adj_struct adiacenza_fuzzy_line(int N,int fuzzy);
adj_struct adiacenza_simple_line(int N);
adj_struct adiacenza_toroidal_lattice(int lato);
adj_struct adiacenza_open_square_lattice(int lato);
adj_struct adiacenza_square_lattice(int lato);
adj_struct adiacenza_from_file(const char *name_vec1, const char *name_vec2);
adj_struct adiacenza_sierpinski(int GEN);
void adiacenza_to_file(const adj_struct & nn);
#endif	/* ADJ_HANDLER_H */

