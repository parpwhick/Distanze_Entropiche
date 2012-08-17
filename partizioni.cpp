/**
 * @file partizioni.cpp
 * @brief Definisce tutti i metodi della classe @ref general_partition e l'essenziale funzione ordered_vector_entropy()
 */
#include <cstdlib>
#include <cstdio>

#include <cmath>
#include <vector>
#include <algorithm>

#include "strutture.h"
#include "adj_handler.h"
#include "partizioni.h"

using std::vector;
using std::pair;
using std::make_pair;
using std::sort;

extern options opts;
extern double *mylog;

/************************************************
 ************************************************
 *            PARTIZIONI COMPLESSE              *
 ************************************************
 ************************************************
 */

/** allocate: funzione helper, assicura di avere memoria allocata per tutto quello
 *  che serve. Da chiamare sempre prima di usare la memoria. Se è gia allocata, non fa
 *  nulla a parte un controllo di sicurezza!
 */
void general_partition::allocate(label_t len) {
    //errore se:
    //len == 0 : cerco di allocare una partizione lunga 0
    //labels è pieno, ma N != len -> probabilmente si sta usando partizioni di lunghezze diverse
    if (len==0 || (!labels.empty() && N != len)) {
        fprintf(stderr, "Allocating for a different length, or length 0\n");
        exit(1);
    }
    N = len;

    //preallocazione dei vettori, se ci si riesce
    try {
        labels.resize(N);
        prev_site.resize(N);
    } catch (std::bad_alloc &e) {
        fprintf(stderr, "Error allocating new partition: %s\n", e.what());
        exit(1);
    }
}

general_partition::general_partition(int len) {
    n = 0;
    N = len;
    entropia_topologica = 0;
    entropia_shannon = 0;

    if (N)
        allocate(N);
}

/**@brief Calcola l'entropia di un range ordinato di dimensione N
 * Gli atomi sono individuati da un'etichetta diversa della precedente, il metodo e' lo stesso usato in
 * @ref general_partition::product() e @ref entropy_binary_partition()
 *
 * @param temp[] Vettore ordinato con le etichette di cui si vuole conoscere l'entropia
 * @param N lunghezza del vettore
 * @return La coppia <(double) H, (int) n_atomi> di entropia e n */
template <typename T> entropy_pair ordered_vector_entropy(const T *temp, int N){
    int label_count = 0;
    double H = 0;
    int mu;
    int begin;

    //the first position always starts an atom
    begin = 0;
    label_count = 1;
    T old_val = temp[0];

    for (label_t i = 1; i < N; i++) {
        //whenever we find a new atom
        if (temp[i] != old_val) {
            //a new atom starts

            //the closed (old)atom's length is calculated
            mu = i - begin;
            //the new one is ready to go
            label_count++;
            begin = i;
            //cache the new label to check
            old_val = temp[i];
            //we add the entropy, with the trick mu>0 and when mu=1 the log is 0
            if (mu > 1)
                H += (double) mu * mylog[mu];
        }
    }
    //the last one, so it's not left hanging
    mu = N - begin;
    H += mu * mylog[mu];

    //normalize the result
    H = -H / N + mylog[N];

    return std::make_pair(H,label_count);
}
template entropy_pair ordered_vector_entropy(const label_t *temp, int N);
template entropy_pair ordered_vector_entropy(const product_t *temp, int N);
template entropy_pair ordered_vector_entropy(const char *temp, int N);
template entropy_pair ordered_vector_entropy(const std::pair<label_t,label_t> *temp, int N);
template entropy_pair ordered_vector_entropy(const uint32_t *temp, int N);

/**
 * Crea un array temporaneo in cui copia i labels per usare il sort, che distruggerebbe l'ordine iniziale, modificando la partizione.
 */
void general_partition::entropy_calculation() {
    //create a temp vector, a copy of "labels"
    //during the sort, temp is destructively changed, can't use the labels!
    vector<label_t> temp(labels);
    sort(temp.begin(),temp.end());

    entropy_pair entropie = ordered_vector_entropy(temp.data(),N);
    
    n = entropie.second;
    entropia_topologica = mylog[n];
    entropia_shannon = entropie.first;
}

/** Implementazione di similitudine a meno di 'epsilon', tramite differenza
 * simmetrica di due insiemi (ottenuti tramite iteratori) ordinati.
 * Ritorna un risultato non appena le differenze sono sufficienti.
 *
 * @return 'false' non appena le differenze sono maggiori di epsilon, 'true' altrimenti.
 */
bool is_similar(Iter_t from1, Iter_t from2, Iter_t to, int epsilon) {
    int differenze = 0;

    while (differenze < epsilon) {
        /*if (differenze > tol) {
            // Se il numero delle differenze è elevato, escludo l'uguaglianza
            return (false);
        } else */
        if ((from1 == to) && (from2 == to)) {
            // Se ho raggiunto la fine di entrambi gli atomi senza troppe differenze
            // allora accetto l'uguaglianza
            return (true);
        } else if ((from1 == to) != (from2 == to)) {
            // Nel caso di atomi con disparita' numerica, continuo a scorrere il piu
            // lungo, per vedere se accumulo abbastanza differenze
            differenze++;
            from1++;
            from2++;
        } else if (*from1>*from2) {
            // Se ho elementi diversi nei due atomi, scorro per riallinearli
            differenze++;
            from1++;
        } else if (*from2>*from1) {
            differenze++;
            from2++;
        } else {
            // I due atomi hanno lo stesso elemento - scorro in avanti entrambi
            from1++;
            from2++;
        }
    }
    return(false);
}

/** Differenza simmetrica tra due insiemi modificata per ritornare non appena
 * si trova una differenza -- ottimizzazione per il caso epsilon == 0 .
 *
 * @return 'true' quando gli elementi nei due range [from,last) sono uguali.
 * @param from1 Iteratore del primo insieme
 * @param from2 Iteratore del secondo insieme
 * @param last1 Iteratore 'end' del primo insieme
 * @param last2 Iteratore 'end' del secondo insieme
 */
bool is_equal(Iter_t from1, Iter_t last1, Iter_t from2, Iter_t last2) {
    while ((from1 != last1) && (from2 != last2)) {
        if (*from1 != *from2)
            return false;
        ++from1;
        ++from2;
    }
    return (from1 == last1) && (from2 == last2);
}

/**
 * Genera una partizione ridotta incompleta, generata da p1, scartando gli elementi
 * uguali (a meno di epsilon) ad atomi in p2.\n
 * La partizione cosi generata è completa, contiene tutte le informazioni necessarie, sovrabbondante
 * rispetto a quanto richiesto dal solo calcolo delle distanze.
 *
 * Epsilon è preso dalla variabile globale opts.epsilon
 * @param p1 Partizione da ridurre
 * @param p2 Partizione con cui confrontare
 * @return Partizione ridotta nell'oggetto chiamante
 */
void general_partition::reduce(const general_partition &p1, const general_partition &p2) {
    //inizializzazioni
    const int & epsilon = opts.epsilon;
    if (p1.n == 0 || p2.n == 0) {
        fprintf(stderr, "REDUCE: Trying to use empty partition\n");
        exit(1);
    }
    allocate(p1.N);
    atomi.resize(N);

    labels.assign(N, -1);
    n=0;

    int common_size = N;
    //per ogni atomo
    for (label_t which = 0; which < p1.n; which++) {
        // Considero l'atomo n-esimo della prima partizione
        label_t size = p1.atomi[which].size;

        // Creazione degli iteratori, per ottenere tutti i siti in atomo1
        Iter_t ii1 = p1.begin(which);
        Iter_t end = p1.end();

        /** Nel caso epsilon > 0, scorro ogni atomo della partizione 2 che ha qualche sito coincidente
         *  con l'atomo della partizione 1. Infatti atomo1 puo' essere eliminato da un qualunque
         *  atomo della partizione2, bisogna controllare in maniera esaustiva usando is_similar()*/
        if (epsilon) {
            //salto se l'atomo è sufficientemente piccolo
            bool almost_equal = 2 * size < epsilon;
            Iter_t ii = ii1;
            /* Ottimizzazione: controlliamo la epsilon-uguaglianza solo una volta
             *  per ogni atomo della partizione 2, quindi teniamo un indice degli
             *  atomi incontrati, e facciamo un controllo solo si trova un atomo nuovo.
             */
            label_t atomo2_old = -1;
            //for sul iteratore dell'atomo 1
            for (; ii != end && !almost_equal; ii++) {
                label_t atomo2 = p2.labels[*ii];
                if (atomo2 == atomo2_old)
                    //controlliamo sito successivo
                    continue;
                if (std::abs(p2.atomi[atomo2].size - size) > epsilon) {
                    //se la differenza delle dimensioni è maggiore di epsilon,
                    //sicuramente non possono essere uguali!
                    almost_equal = false;
                    break;
                }
                if (is_similar(ii1, p2.begin(atomo2), end, epsilon)) {
                    //test vero e proprio di similarita' a meno di epsilon, i precedenti
                    //test erano ottimizzazioni
                    almost_equal = true;
                    break;
                }
                atomo2_old = atomo2;
            }
            if (almost_equal)
                //atomo1 è saltato, perche è stato trovato un corrispondente
                //nella partizione 2
                continue;
        } else {
            label_t atomo2 = p2.labels[*ii1];
            /// Nel caso epsilon==0, l'uguaglianza secca è implementata velocemente in is_equal()
            if(is_equal(ii1,end,p2.begin(atomo2),end))
                continue;
        }
        // common_size rappresenta la dimensione degli atomi scartati, diminuisce ogni volta che teniamo un atomo.
        common_size -= size;
        // Per ogni atomo tenuto si tiene conto del suo contributo entropico.
        entropia_shannon += size * mylog[size];

        /** Se l'intersezione è banale, intersechiamo il corrispondente del fattore dicotomico con i precedenti.
         *  Questa intersezione in realta' non è necessaria - si scrive direttamente il risultato.
         *  L'assegnazione di @ref labels della partizione ridotta è rapida: corrispondono
         *  all'indice delle partizioni dicotomiche tenute, @ref n. Il caso dell'atomo di fondo è trattato a parte.\n
         *  Similmente @ref prev_site all'interno degli atomi è uguale, ricopiamo quindi sito per sito.
         */
        for (; ii1 != end; ii1++){
            labels[*ii1] = n;
            prev_site[*ii1] = p1.prev_site[*ii1];
        }
        ///Gli atomi sono "copiati" dalla partizione1 alla partizione ridotta, per questo hanno gli stessi parametri.
        atomi[n]=p1.atomi[which];
        n++;
    }

    /** In fine si tiene conto dell'atomo sconnesso di background, formato dai pezzi comuni (o simili)
     * delle partizioni.
     * La dimensione è semplice da calcolare e nel caso di molti scarti da' un contributo entropico
     * significativo. Le informazioni su questo atomo non sono presenti nella partizione iniziale, ma vanno calcolate appositamente.
     * Ad esempio si cerca sito per sito gli appartenenti all'atomo di fondo e si assegna opportunamente @c prev_site di volta in volta.
     */

    //se la parte comune ha dim > 0...
    if (common_size) {
        atomi[n].size = common_size;
        label_t dove;
        for(dove=0; dove<N; dove++)
            if(labels[dove] == -1)
                break;

        //mettere a posto il primo sito
        labels[dove] = n;
        atomi[n].start = dove;
        label_t prev = dove;
        prev_site[dove] = dove;

        //cercare il resto dei siti, e collegarli ai precenti
        for(dove++; dove < N; dove++){
            if(labels[dove] == -1){
                labels[dove] = n;
                prev_site[dove] = prev;
                prev = dove;
            }
        }
        //chiudiamo con l'ultimo trovato e segnamo di aver trovato un atomo in piu'.
        atomi[n].end = prev;
        n++;
    }
    // entropia della parte comune e normalizzazione
    entropia_shannon += common_size * mylog[common_size];
    entropia_shannon = -entropia_shannon / N + mylog[N];
    entropia_topologica = mylog[n];
}


/** findroot ricorsivamente trova il primo sito del cluster, detto 'radice' o 'root', applicando
 * path compression (in questa versione nonricorsiva usando path-halving), per ottimizzare
 * future chiamate alla stessa funzione.
 */
template <typename data_t> inline int findroot(int i, data_t *ptr) {
    int r, s;
    r = s = i;
    while (ptr[r] >= 0) {
        ptr[s] = ptr[r];
        s = r;
        r = ptr[r];
    }
    return r;
}

/** Calcolo della partizione a partire da un vettore di stato configuration[], di adiacenza adj --
 *  necessaria per interpretare il vettore di stato.
 *
 * Il calcolo è fatto tramite percolazione sulla struttura data da adj, mettendo link
 *  solo laddove i siti vicini hanno gli stessi valori. Il numero dei vicini non è fissato.
 * Alla fine si chiama 'relabel' per riempire tutte le informazioni sulla partizione oltre
 *  al label di ciascun atomo.
 *
 * @param configuration Vettore di configurazione/stato
 * @param adj Struttura di adiacenza, contenente le informazioni sulla topologia
 * @return Partizione completa ottenuta dalla configurazione fornita
 */
template <typename data_t>
void general_partition::from_configuration(const data_t *configuration, const adj_struct & adj) {
    label_t s1, s2;
    label_t r1, r2;
    int z;

    allocate(adj.N);
    /// Ogni sito in partenza è un atomo a sè stante, da unire poi agli altri, con label -1.
    labels.assign(N, -1);

    // percolazione e primi labels
    for (s1 = 0; s1 < N; s1++) {
        ///Per ogni sito si trova il cluster di appartenenza
        r1 = findroot(s1, &labels[0]);

        z = adj.fetch(s1);
        for (int j = 0; j < z; j++) {
            ///Si guarda ogni vicino
            s2 = adj.vicini[j];

            ///e se fanno parte dello stesso cluster
            if (s2 >= s1 || s2 < 0 || configuration[s1] != configuration[s2])
                continue;

            r2 = findroot(s2, &labels[0]);
            ///si unisce l'atomo del vicino al proprio, aggiungendo il "ramo" alla propria radice.
            if (r1 != r2) {
                if (labels[r1] >= labels[r2]) {
                    labels[r2] += labels[r1];
                    labels[r1] = r2;
                    r1 = r2;
                } else {
                    labels[r1] += labels[r2];
                    labels[r2] = r1;
                }
            }
            //next neightbor
        }
        //next site
    }

    this->relabel();

    if (opts.graphics && (opts.topologia == RETICOLO_2D)) {
        static int imagecount = 0;
        char filename[255];
        imagecount++;
        sprintf(filename, "reticolo%03d.ppm", imagecount);
        ppmout2(configuration, &labels[0], opts.lato, filename);
    }
}
template void general_partition::from_configuration(const int *configuration, const adj_struct & adj);
template void general_partition::from_configuration(const char *configuration, const adj_struct & adj);

/** A partire dai label della percolazione, costruisce l'elenco degli atomi, il vettore prev_site;
 *  La partizione finale non usa labels che indicano il cluster con il primo sito di appartenza, ma sono rinumerati,
 *  usando come etichette l'indice progressivo dall'atomo.
 */
void general_partition::relabel() {
    entropia_shannon = 0;
    n = 0;

    // 0-conto nr. atomi diversi per il solo scopo di allocare ottimalmente
    // 1-assicurazione che prev_site punta al cluster di appartenenza
    for(label_t i = 0; i < N; i++){
        n += labels[i] < 0;
        //a questo punto prev_site contiene il primo sito appartenente al cluster
        prev_site[i] = findroot(i, &labels[0]);
    }
    atomi.resize(n);
    n = 0;

    /**
      - creazione array atomi
      - inizializzazione ogni elemento
      - creazione indice (label atomo) <--> root
      - calcolo entropia a partire dai size nei root
     */
#define ATOMO atomi[n]
    for (label_t i = 0; i < N; i++) {
        if (labels[i] < 0) {
            entropia_shannon += -labels[i] * mylog[-labels[i]];
            ATOMO.size = -labels[i];
            ATOMO.end = i;
            ATOMO.start = i;
            labels[i] = n;
            n++;
        }
    }
    entropia_topologica = mylog[n];
    entropia_shannon = -entropia_shannon / N + mylog[N];

    /**
      - relabeling secondo l'indice dell'atomo, non del sito di appartenenza
      - creazione del collegamento prev_site(che usa il precedente, non il primo com'era prima)
         e aggiornamento atom.end, che indica l'ultimo sito trovato apparentenente all'atomo
     */
    for (label_t i = 0; i < N; i++) {
        int atom_label = labels[prev_site[i]];
        labels[i] = atom_label;
        prev_site[i] = std::min(atomi[atom_label].end, i);
        atomi[atom_label].end = i;
        label_t & atom_start = atomi[atom_label].start;
        atom_start = std::min(atom_start, i);
    }
}

/** Funzione come from_configuration (percolazione) specializzata per il caso di una struttura 2-dimensionale.
 * La partizione intersezione è definita come la partizione piu' fine, minore sia di p1 che di p2.
 * Queste caratteristiche si ottengono creando una partizione che ha come vicini l'unione dei vicini
 * di p1 e di p2. \n
 * Il caso ottimizzato è bidimensionale: basta prendere l'unione di p1.prev_site[] e p2.prev_site[].
 *  
 * @return La partizione risultante è l'intersezione (op. simmetrica) di p1 con p2.
 */
void general_partition::common_subfactor(const general_partition &p1, const general_partition &p2) {
    label_t s1, s2;
    allocate(p1.N);

    //set all labels to -1
    labels.assign(N, -1);

    // percolazione e primi labels
    for (s1 = 0; s1 < N; s1++) {
        int r2;
        int r1 = findroot(s1, &labels[0]);

        //usa come primo vicino il sito dalla partizione 1
        s2 = p1.prev_site[s1];
        if (s1 != s2 && s2 >= 0) {
            r2 = findroot(s2, &labels[0]);
            if (r1 != r2) {
                if (labels[r1] >= labels[r2]) {
                    labels[r2] += labels[r1];
                    labels[r1] = r2;
                    r1 = r2;
                } else {
                    labels[r1] += labels[r2];
                    labels[r2] = r1;
                }
            }
        }

        //usa come secondo vicino il sito dalla partizione 2
        s2 = p2.prev_site[s1];
        if (s1 != s2 && s2 >= 0) {
            r2 = findroot(s2, &labels[0]);
            if (r1 != r2) {
                if (labels[r1] >= labels[r2]) {
                    labels[r2] += labels[r1];
                    labels[r1] = r2;
                    r1 = r2;
                } else {
                    labels[r1] += labels[r2];
                    labels[r2] = r1;
                }
            }
        }
    }

    this->relabel();

    // Grafico il reticolo risultante, se richiesto
    if (opts.graphics && (opts.topologia == RETICOLO_2D)) {
        static int imagecount = 0;
        char filename[255];
        imagecount++;
        sprintf(filename, "comune%03d.ppm", imagecount);
        ppmout2(&p1.labels[0], &p2.labels[0], opts.lato, filename);
        imagecount++;
        sprintf(filename, "comune%03d.ppm", imagecount);
        ppmout(&labels[0], opts.lato, filename);
    }
}

/** Costruisce la matrice di adiacenza sparsa intra-cluster (utile per il laplaciano della partizione).
 * La matrice di adiacenza è cosi costruita:\n
 * @verbatim
  A(i,j) =  / 1 se i è nello stesso cluster di j
            \ 0 altrimenti
   @endverbatim
 * Il risultato è una matrice poco sparsa, ma puo' dare info spettrali sui clusters.
 * 
 * @param () Nessun parametro di input
 * @return Stampa la matrice corrispondente alla partizione.
 */
void general_partition::print_cluster_adjacency() {
    //label_t è definito in strutture.h come int32
    std::vector<label_t> riga;
    std::vector<label_t> colonna;
    int totale = 0;
    FILE *vec1 = fopen("vector1.bin", "wb");
    FILE *vec2 = fopen("vector2.bin", "wb");
    //per ogni sito
    for (label_t which = 0; which < N; which++) {
        /** Per ogni sito appartenente al reticolo, recupero le informazioni
           sul cluster di appartenenza (atomo in questa nomenclatura).
           Cio' è necessario per ottenere tutti i siti (ordinati) del cluster
           cercato con efficienza massima. */
        const atom &atomo = atomi[labels[which]];
        int quanti = atomo.size;
        riga.resize(quanti);
        colonna.resize(quanti);

        int sito = 0;
        //scorro tutti i siti appartenenti allo stesso atomo,
        //con l'iteratore ii. Il sito corrispondente è  *ii
        for (Iter_t ii = this->begin(atomo); ii != this->end(); ii++) {
            /* gli elementi nonnulli della matrice di adiacenza A(i,j)
               sono in (which, *ii), salvo i valori delle righe e delle colonne
               corrispondenti in due vettori */

            //salto elemento diagonale
            if (which == *ii)
                continue;
            riga[sito] = which + 1;
            colonna[sito] = *ii + 1;
            sito += 1;
        }
        //stampa i vettori cosi costruiti!!!!!!!!!!
        fwrite(&riga[0], sizeof (label_t), sito, vec1);
        fwrite(&colonna[0], sizeof (label_t), sito, vec2);
        totale += sito;

        /** Creazione vettori di adiacenza per matrice sparse in stile Matlab
         La funzione per caricare i dati cosi creati è:
        @code
        function adiacenza=load_sierpinski()
            indici_riga=fread(fopen('vector1.bin','r'),inf,'int32');
            indici_colonna=fread(fopen('vector2.bin','r'),inf,'int32');
            N=max(max(indici_riga),max(indici_colonna));
            adiacenza=sparse(indici_riga,indici_colonna,1,N,N);
        end
        @endcode
         */
    }
    fclose(vec1);
    fclose(vec2);
    printf("Elementi nonnulli della matrice di adiacenza: %d\n", totale);
}

/** Calcolo della partizione prodotto (simmetrico).
 * Il risultato è un oggetto partizione completo, con tutte le proprieta' definite.
 * Per il calcolo di una distanza è eccessivo - meglio usare la classe distance,
 *  che implementa lo stesso codice troncato e senza temporanei inutili.
 *
 * @param p1 Primo fattore
 * @param p2 Secondo fattore
 * @return L'oggetto partizione prodotto
 */
void general_partition::product(const general_partition & p1, const general_partition & p2) {
    label_t mu;
    label_t begin;
    allocate(p1.N);
    atomi.resize(N);

    /** I label del prodotto sono rappresentati dalle coppie (pair<label_t, label_t>)
     * di label dei fattori. Il vettore @c label_index contiene la coppia <label_prodotto, indice>.\n
     * Per riconoscere gli atomi, ordino il vettore label_index:
     *  il risultato sarà un vettore ordinato, per cui riconosciamo gli atomi come elementi
     *  contigui con lo stesso label. Il label è temporaneo, gli atomi avranno indice
     *  progressivo.
     * Esempio:
     * (0,1) - 1
     * (0,2) - 2
     * (1,0) - 3
     * (1,0) - 3
     * (1,2) - 4
     */

    /** Il vettore @e label_index contiene i seguenti membri:
     *  - label_index[i].first - indice (temp) del prodotto.
     *  - label_index[i].second - indica la posizione del sito nella partizione,
     *     è utilizzato per indicizzare i cambiamenti e gli accessi ai vettori labels[], ecc.
     */
    vector<pair<uint32_t, label_t> > label_index(N);
    n = 0;
    entropia_shannon = 0;

    label_t count = 0;
    label_t old_count = 0;
    int fiddle = 0;
    /** Gli elementi sono ordinati atomo-per-atomo. Per ogni sito dell'atomo1, si controllano
     * tutti i label di appartenenza per la partizione2 - i diversi atomi incontrati
     * formano un diverso atomo della partizione prodotto!
     * Si utilizza sempre il sort, ma stavolta per ogni atomo della partizione 1,
     * riducendo trasticamente i tempi, poiché sono state struttate le informazioni sulla
     * disposizione dei vari elementi.
     */
    for (label_t atom_index = 0; atom_index < p1.n; atom_index++) {
        /* fiddle factor serve per distinguere atomi con il primo indice diverso:
         * infatti salviamo solo le etichette dell'atomo2, quindi se queste sono identiche,
         * ma differiscono per l'atomo1, vanno distinte!
         * Si vuole evitare il seguente caso:
         * (0,1)
         * (1,1)
         * che sarebbero riconosciuti come uguali. Per cui si setta il primo bit del label,
         * alternativamente 0 e 1, per distinguere i vari atomi.
         */
        fiddle = (fiddle) ? 0 : (1<<(sizeof(label_t)*8-1));
        for (Iter_t ii = p1.begin(atom_index); ii != p1.end(); ii++)
            label_index[count++] = make_pair(fiddle | p2.labels[*ii], *ii);

        sort(&label_index[old_count], &label_index[count]);
        old_count = count;
    }
    /* Dopo il sort, riconosciamo gli atomi. Il primo sito sicuramente inizia un nuovo
     * atomo, per i successivi riconosciamo un atomo non appena il label cambia
     */
    begin = 0;
    uint32_t old_val = label_index[0].first;
    
    //we skip over the first value, so initialize here
    labels[label_index[0].second] = n;
    prev_site[label_index[0].second] = label_index[0].second;

    for (label_t i = 1; i < N; i++) {
        label_t & LAST_POS = label_index[i - 1].second;
        label_t & THIS_POS = label_index[i].second;
        //whenever we find a new atom
        if (label_index[i].first != old_val) {
            //a new atom starts

            //the closed (old)atom's length is calculated
            mu = i - begin;
            atomi[n].start = label_index[begin].second;
            atomi[n].size = mu;
            atomi[n].end = LAST_POS;
            //the new one is ready to go
            //prev_site del sito corrente è il sito stesso
            prev_site[THIS_POS] = THIS_POS;
            n++;
            begin = i;
            //cache the new label to check
            old_val = label_index[i].first;

            //we add the entropy, with the trick mu>0 and when mu=1 the log is 0
            if (mu > 1)
                entropia_shannon += (double) mu * mylog[mu];
        } else {
            //prev site del sito corrente è il precedente trovato nello stesso atomo
            prev_site[THIS_POS] = LAST_POS;
            //il sito attuale è sicuramente l'ultima posizione trovata per l'atomo corrente
            atomi[n].end = THIS_POS;
        }
        //relabel anyway con l'indice di atomo corrente
        labels[THIS_POS] = n;
    }
    //the last one, so it's not left hanging
    mu = N - begin;
    atomi[n].start = label_index[begin].second;
    atomi[n].size = mu;
    atomi[n].end = label_index[N - 1].second;
    entropia_shannon += mu * mylog[mu];
    //il numero di atomi trovati è l'indice corrente +1
    n++;

    //normalize the entropy
    entropia_shannon = -entropia_shannon / N + mylog[N];
    entropia_topologica = mylog[n];
}

/**
 * @brief Genera un vettore @c next_site a partire da un vettore @c prev_site
 * @param prev_site Vettore rappresentativo e membro della classe @ref general_partition
 * @return Vettore @c next_site, che fornisce il sito successivo all'indice del vettore
 */
void general_partition::generate_forward_linking(){    
    next_site.resize(prev_site.size());

    for(label_t i=0; i < (label_t)prev_site.size(); i++){
        next_site[i]=i;
        if(prev_site[i] != i)
            next_site[prev_site[i]]=i;
    }
}

bool general_partition::operator<=(const general_partition & p2) const {
    for (label_t atom_label2 = 0; atom_label2 < p2.n; atom_label2++) {
       Iter_t ii = p2.begin(atom_label2);
       const label_t atom_label1 = labels[*ii];
       for (; ii != end(); ii++)
           if(atom_label1 != labels[*ii])
               return false;
    }
    return true;
}

bool general_partition::operator==(const general_partition & p2) const {
    //metodo2 - confronto tra strutture di tipo "atom"
    for (label_t atom_label1 = 0; atom_label1 < n; atom_label1++) {
        const atom & atomo1 = atomi[atom_label1];
        const atom & atomo2 = p2.find_atom(atomo1);
        if(! (atomo1 == atomo2))
            return false;
    }
    return true;
}

bool general_partition::consistency_check()  {
    bool failed = false;

    //counting the real number of atoms by the labels
    vector<label_t> temp(labels);
    sort(temp.begin(),temp.end());
    entropy_pair entropie = ordered_vector_entropy(temp.data(),N);
    int n_calculated = entropie.second;
    double H = entropie.first;
    if(n!=n_calculated){
        printf("Number of atoms is wrong, is %d, should be %d\n",n,n_calculated);
        failed = true;
    }
    if(abs(H-entropia_shannon)>1e-9){
        printf("Entropy is wrong, is %.5f, should be %.5f\n",H,entropia_shannon);
        failed=true;
    }

    //atoms & labels
    for (label_t atom_label1 = 0; atom_label1 < n; atom_label1++) {
        for (Iter_t ii = reverse_iterator(atomi[atom_label1].end, & prev_site[0]); ii != end(); ii++) {
            if(*ii >= N){
                printf("Iterator is wrong, returned next site %d\n",*ii);
                print_array(&prev_site[0], N, "prev_s");
                failed = true;
                break;
            }
            if (labels[*ii] != atom_label1) {
                printf("Self consistency, backward iterator check 1 failed: %d instead of %d, at %d\n", labels[*ii], atom_label1, *ii);
                printf("Atom %d, from %d to %d\n",atom_label1,atomi[atom_label1].start,atomi[atom_label1].end);
                failed = true;
                break;
            }
        }
    }
    //checking that by the iterators, all the partition can be covered
    int count = 0;
    for (label_t atom_label1 = 0; atom_label1 < n; atom_label1++)
        for (Iter_t ii = reverse_iterator(atomi[atom_label1].end, & prev_site[0]); ii != end(); ii++)
            count++;
    if (count != N) {
        printf("Self consistency, backward iterator check 2 failed: count is %d (/%d)\n", count, N);
        for (label_t atom_label1 = 0; atom_label1 < n; atom_label1++)
        printf("Atom %d, from %d to %d\n",atom_label1,atomi[atom_label1].start,atomi[atom_label1].end);
        failed = true;
    }

    generate_forward_linking();

    //checking that the end and start of every atom are well set
    for (label_t atom_label1 = 0; atom_label1 < n; atom_label1++) {
        int inizio = atomi[atom_label1].start;
        if(inizio !=  prev_site[inizio]){
            printf("[%2d] should begin at %d, begins at %d\n",atom_label1,inizio,prev_site[inizio]);
            failed=true;
        }
        int finish = atomi[atom_label1].end;
        if(finish !=  next_site[finish]){
            printf("[%2d] should end at %d, ends at %d\n",atom_label1,finish,next_site[finish]);
            failed=true;
        }
    }
    //operators <= and == check
    bool uguale =  *this == *this;
    bool minore_uguale = *this <= *this;
    if(!uguale || !minore_uguale){
        printf("Failed equality test: == %d, <= %d\n",uguale,minore_uguale);
        failed = true;
    }

    //print summary of the partition
    if (failed) {
        label_t quanto = std::min(N, (label_t) 50);
        printf("n: %d\n", n);
        print_array(&labels[0], quanto, "labels");
        print_array(&prev_site[0], quanto, "prev_s");
        print_array(&next_site[0], quanto, "next_s");
    }
    return failed;
}
