/**
 * @file partizioni.h
 * @brief Header che definisce le partizioni generali (@ref general_partition) e il loro iteratore (general_partition::Iterator)
 */

#ifndef PARTIZIONI_H
#define	PARTIZIONI_H

#include <vector>
#include <string>
#include <stdint.h>
#include "adj_handler.h"

#ifndef label_t
#define int32_t label_t
#define uint64_t product_t
#endif

using std::vector;
///Informazioni essenziali su una partizione
class basic_partition {
public:
    ///Numero di atomi
    label_t n;
    ///Volume della partizione (in nr. di siti)
    label_t N;

    ///Entropia di Shannon
    double entropia_shannon;
    ///Entropia topologica
    double entropia_topologica;
};

///Partizione lineare semplice, rappresentabile come un vettore binario
class linear_partition : public basic_partition {
public:
    ///Binary atom beginning mark codes {1,0,0,0,1,1,...}
    int *binary;

    ///Costruttore a partire da una configurazione
    template <typename T> linear_partition(const T*seq, int len);
    ///Crea la partizione a partire da una sequenza data
    template <typename T> void fill(const T*seq, int len);
    ///Stampa il contenuto della partizione e informazioni
    void print();

    ///Costruttore
    linear_partition(int len = 0) {
        n = 0;
        N = len;
        binary = new int[len];
    }

    ~linear_partition() {
        if (N)
            delete[] binary;
        n = N = 0;
    }
};

/** @brief Informazioni minime per individuare un atomo in maniera efficiente.
 *
 * Tramite queste informazioni e il vettore prev_site[] di ogni partizione, è possibile
 * percorrere in modo ottimale ogni atomo (in una direzione sola).
 */
class atom {
public:
    ///Dimensione dell'atomo
    label_t size;
    ///Primo sito di appartenenza
    label_t start;
    ///Ultimo sito. Il vettore @c prev_site collega @e end a @e start attraverso tutto l'atomo.
    label_t end;

    bool operator==(const atom & other) const {
        return (size == other.size && start == other.start && end == other.end);
    }
};

/**@brief Partizione di tipo completamente generale.
 *
 * La classe contiene i metodi necessari partizionare un qualunque vettore (fornita l'adiacenza),
 * ottenere le informazioni sulla partizione,
 * eseguire tutte le operazioni necessarie tra due partizioni - nel qual caso la partizione è il risultato del
 * prodotto/intersezione/riduzione.
 */
class general_partition : public basic_partition {
    friend class distance;
private:
    ///Alloca i vettori prima di poter procedere a qualunque calcolo (chiamato automaticamente).
    void allocate(label_t len);
    ///Calcola l'entropia nel caso di una partizione nota solo tramite i suoi @ref labels. Non usata attualemente.
    void entropy_calculation();
    ///Dopo la percolazione, corregge i @ref labels e riempie le strutture necessarie.
    void relabel();
    ///Genera next_site a partire da prev_site
    void generate_forward_linking();


    ///Etichetta (ovvero l'indice dell'atomo partendo da zero), corrispondente al sito i-esimo: labels[i]
    vector<label_t> labels;
    ///Sito precedente (o uguale) facente parte dello stesso atomo del sito i-esimo: prev_site[i]
    vector<label_t> prev_site;
    ///Sito precedente (o uguale) facente parte dello stesso atomo del sito i-esimo: prev_site[i]
    vector<label_t> next_site;
    ///Vettore di atomi, numerati in [0,n). Contiene le informazioni essenziali per manipolare tutti gli atomi della partizione. L'allocazione è ottimizzata, il vettore è lungo @e n.
    vector<atom> atomi;

    //L'utilizzo del copy operator è vietato.
    //general_partition(const general_partition &);
    //L'operatore di assegnazione è vietato.
    //general_partition & operator=(const general_partition &);
public:

    ///Costruttore della classe. Imposta il numero di atomi a 0, alloca la dimensione @e len se necessario.
    general_partition(int len = 0);

    ///Genera la partizione a partire una configurazione e una struttura di adiacenza
    template <typename data_t> void from_configuration(const data_t *configuration, const adj_struct & adj);
    ///A partire da due partizioni p1 e p2, genera nell'oggetto la partizione ridotta
    void reduce(const general_partition &p1, const general_partition &p2);
    ///Genera la partizione intersezione di p1 con p2
    void linear_intersection(const general_partition &p1, const general_partition &p2);
    ///Genera la partizione prodotto di p1 con p2
    void product(const general_partition & p1, const general_partition & p2);

    void print();
    ///Stampa matrice di adiacenza intra-cluster
    void print_cluster_adjacency();
    ///Restituisce un riferimento read-only per accedere alle etichette
    const label_t * show_labels(){ return labels.data();}

    ///A partire da un atomo generico, restituisce l'atomo che ha il primo sito in comune
    const atom& find_atom(const atom &atomo1) const {
        return atomi[labels[atomo1.start]];
    }
    ///Restituisce l'atomo a cui il sito appartiene
    const atom& find_atom(const label_t posizione) const {
        return atomi[labels[posizione]];
    }

    /**@brief L'iteratore serve per scorrere rapidamente e senza errori tutti i siti di un atomo
     *
     * L'iteratore scorre tutti i siti di un atomo nella direzione indicata dal vettore che gli viene passato
     * durante la costruzione. Nel caso di default, e' usato @a prev_site per cui l'iteratore scorre all'indietro.
     * E' possibile usare l'iteratore per scorrere in avanti, usando un @a next_site.
     *
     * Per avere un ciclo for su tutto l'atomo a cui @c sito appartiene, si puo' usare:
     * @code
     * for(Iter_t ii = partizione.begin(partizione.labels[sito]); ii != partizione.end(); ii++)
     *          dosomething(*ii);
     * @endcode
     */
    /// Eredita della classe std::iterator per avere tutti i typedef necessari automaticamente
    class Iterator : public std::iterator<std::forward_iterator_tag, label_t>{

    private:
        ///Sito restituito dall'iteratore
        label_t _site;
        ///Riferimento interno a prev_site, per poter scorrere i siti
        const label_t *_next;

    public:
        //costruttore, chiamabile solo dalla classe partizione
        Iterator(int dove, const label_t *vicini) : _site(dove), _next(vicini) {
        };

        ///Per leggere il sito a cui punta l'iteratore
        int operator*() {
            return (_site);
        }

        ///Confronto tra due iteratori
        bool operator==(Iterator due) {
            return (_site == due._site);
        }

        ///Diseguaglianza tra due iteratori
        bool operator!=(Iterator due) {
            return (_site != due._site);
        }

        ///Incremento intelligente: se esiste un altro sito appartenente all'atomo, allora viene caricato. Altrimenti l'iteratore restituisce @ref end.
        int operator++() {
            if (_site == _next[_site])
                _site = -1;
            else
                _site = _next[_site];
            return (_site);
        }

        int operator++(int) {
            return (operator++());
        }

    };
    ///In questo caso gli iteratori inversi sono uguali agli iteratori
    typedef Iterator reverse_iterator;

    ///Restituisce l'iteratore per tutto l'atomo a cui il sito indicato appartiene
    reverse_iterator begin(const int where) const {
        return reverse_iterator(atomi[where].end, & prev_site[0]);
    }
    ///Restituisce l'iteratore per l'atomo richiesto
    reverse_iterator begin(const atom &where) const {
        return reverse_iterator(where.end, & prev_site[0]);
    }
    ///Iteratore banale per riconoscere la fine di un atomo
    reverse_iterator end() const {
        return reverse_iterator(-1, 0);
    }
};

typedef general_partition::reverse_iterator Iter_t;



#endif	/* PARTIZIONI_H */

