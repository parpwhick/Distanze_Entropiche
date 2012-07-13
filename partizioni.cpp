#include <cstdlib>
#include <cstdio>
#include <cassert>

#include <cmath>
#include <vector>
#include <algorithm>

#include "strutture.h"
#include "adj_handler.h"

using std::vector;
using std::pair;
using std::make_pair;
using std::sort;

extern options opts;
extern double *mylog;
int *colore;

static rand55 gen;

template <typename data_t> void print_array(const data_t *array, int len, const char *nome) {
    printf("%s [%3d", nome, array[0]);
    for (int i = 1; i < len; i++)
        printf(",%3d", array[i]);
    printf("]\n");
}
template void print_array(const label_t *grid, int sz, const char *filename);

template <typename T>
void print_square_lattice(const T* valori, int lato) {
    for (int i = 0; i < lato; i++) {
        for (int j = 0; j < lato; j++) {
            printf("%2d ", valori[j * lato + i]);
        }
        printf("\n");
    }
    printf("\n");
}


#define COL_MAX 1000

template <typename T>
void ppmout(const T *grid1, int sz, const char *filename) {

    if (colore == 0) {
        colore = new int[COL_MAX];
        for (int i = 0; i < COL_MAX; i++) {
            colore[i] = gen.rand_long();
        }
        colore[0] = 0x0F5A3A1F; // blue almost black
    }

    int MULT;
    if (sz > 500)
        MULT = 1;
    else
        MULT = 500 / sz + 1;


    FILE *fout = fopen(filename, "w");
    fprintf(fout, "P6\n %d %d\n 255\n", MULT*sz, MULT * sz);


    for (int rg = 0; rg < sz; rg++) {
        for (int i = 0; i < MULT; i++)
            for (int cl = 0; cl < sz; cl++) {
                int sito = grid1[rg * sz + cl] % COL_MAX;
                int color = colore[sito];
                for (int j = 0; j < MULT; j++)
                    fwrite(&color, 3, 1, fout);
            }
    }
}

template <typename T, typename U>
void ppmout2(const T *grid1, const U* grid2, int sz, const char *filename) {
    //Se l'array dei colori non e' inizializzato - riempiamolo!
    if (colore == 0) {
        colore = new int[COL_MAX];
        for (int i = 0; i < COL_MAX; i++) {
            colore[i] = gen.rand_long();
        }
        colore[0] = 0x0F5A3A1F; // blue almost black
    }

    int MULT;
    if (sz > 500)
        MULT = 3;
    else
        MULT = 500 / sz + 1;


    int black = 0x20202020;
    int intermezzo = 100; //pixels tra i pannelli

    FILE *fout = fopen(filename, "w");
    fprintf(fout, "P6\n %d %d\n 255\n", 2 * MULT * sz + intermezzo, MULT * sz);

    for (int rg = 0; rg < sz; rg++) {
        for (int i = 0; i < MULT; i++) {
            //primo reticolo
            for (int cl = 0; cl < sz; cl++) {
                int sito = grid1[rg * sz + cl] % COL_MAX;
                int color = colore[sito];
                for (int j = 0; j < MULT; j++)
                    fwrite(&color, 3, 1, fout);
            }
            //10 pixels neri in mezzo
            for (int j = 0; j < intermezzo; j++)
                fwrite(&black, 3, 1, fout);

            //secondo reticolo
            for (int cl = 0; cl < sz; cl++) {
                int sito = grid2[rg * sz + cl] % COL_MAX;
                int color = colore[sito];
                for (int j = 0; j < MULT; j++)
                    fwrite(&color, 3, 1, fout);
            }
        }
    }
}

template void ppmout(const int32_t *grid, int sz, const char *filename);
template void ppmout(const uint64_t *grid, int sz, const char *filename);
template void ppmout2(const int *grid1, const int* grid2, int, const char *);
template void ppmout2(const uint64_t *grid1, const unsigned long* grid2, int, const char *);

/************************************************
 ************************************************
 *             PARTIZIONI SEMPLICI              *
 ************************************************
 ************************************************
 */


void linear_partition::print() {
    if (opts.verbose > 3) {
        // {1,0,0,1,0,1,...}
        print_array(binary, N, "Binary");
    }
    fprintf(stdout, "Partitions[n]: %d, Shannon %f, Topological %f\n", n, entropia_shannon,
            entropia_topologica);
}

template <typename T>
linear_partition::linear_partition(const T* seq, int len) {
    this->fill(seq, len);
}

template <typename T>
void linear_partition::fill(const T* seq, int len) {
    //Total length of the partition is equal to the sequence
    N = len;
    binary = new int[len];

    //first one always start an atom
    binary[0] = 1;

    //checking for a different symbol from the one before
    //when that happens => new atom!
    for (int i = 1; i < len; i++)
        binary[i] = seq[i] != seq[i - 1];

    //number of atoms found
    n = 0;
    for (int i = 0; i < len; i++)
        n += binary[i];

    entropia_topologica = mylog[n];
    entropia_shannon = entropy_binary_partition(binary, N);
}
template void linear_partition::fill(const char *, int);

/************************************************
 ************************************************
 *            PARTIZIONI COMPLESSE              *
 ************************************************
 ************************************************
 */

/* allocate: funzione helper, assicura di avere memoria allocata per tutto quello
 *  che serve. Da chiamare sempre prima di usare la memoria. Se e' gia allocata, non fa
 *  nulla a parte un controllo di sicurezza!
 */
void general_partition::allocate(label_t len) {
    //errore se:
    //len == 0 : cerco di allocare una partizione lunga 0
    //labels e' pieno, ma N != len -> probabilmente si sta usando partizioni di lunghezze diverse
    if (len==0 || (!labels.empty() && N != len)) {
        fprintf(stderr, "Allocating for a different length, or length 0\n");
        exit(1);
    }
    N = len;

    //preallocazione dei vettori, se ci si riesce
    try {
        labels.reserve(N);
        prev_site.reserve(N);
    } catch (std::bad_alloc &e) {
        fprintf(stderr, "Error allocating new partition: %s\n", e.what());
        exit(1);
    }
}

general_partition::general_partition(int len) {
    n = 0;
    allocated_n = 0;
    N = len;
    entropia_topologica = 0;
    entropia_shannon = 0;

    if (N)
        allocate(N);
}

/* A partire da una partizione con label assegnati, calcola l'entropia
 *  con il metodo del sort. Calcola entrambi i tipi di entropia.
 * Funzione al momento non utilizzata (rimane per completezza)
 */
void general_partition::entropy_calculation() {
    int label_count = 0;
    double H = 0;
    int mu;
    int begin;

    //create a temp vector, a copy of "labels"
    //during the sort, temp is destructively changed, can't use the labels!
    vector<label_t> temp(labels);
    sort(temp.begin(),temp.end());

    //the first position always starts an atom
    begin = 0;
    label_count = 1;
    int old_val = temp[0];

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

    entropia_topologica = mylog[label_count];
    n = label_count;
    entropia_shannon = H;
}

/* Implementazione di similitudine a meno di 'epsilon', tramite differenza
 * simmetrica di due insiemi (ottenuti tramite iteratori) ordinati.
 * L'implementazione e' custom per ritornare un risultato non appena le differenze
 * sono sufficienti a dare una risposta immediata
 */
bool is_atom_similar(Iter_t from1, Iter_t from2, Iter_t to, int epsilon) {
    int differenze = 0;

    while (differenze < epsilon) {
        /*if (differenze > tol) {
            // Se il numero delle differenze e' elevato, escludo l'uguaglianza
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

/* Versione di differenza simmetrica tra due insiemi per ritornare non appena
 * si trova una differenza -- ottimizzazione per il caso epsilon == 0 .
 */
bool is_atom_equal(Iter_t from1, Iter_t from2, Iter_t to) {
    while (true) {
        if ((from1 == to) != (from2 == to)) {
            return (false);
        } else if ((from1 == to) && (from2 == to)) {
            // Se ho raggiunto la fine di entrambi gli atomi senza troppe differenze
            // allora accetto l'uguaglianza
            return (true);
        } else if (*from1 != *from2) {
            return(false);
        } else {
            // I due atomi hanno lo stesso elemento - scorro in avanti entrambi
            from1++;
            from2++;
        }
    }
}

/* general_partition::reduce (part1, part2)
 * Genera una partizione ridotta incompleta, generata da part1, scartando gli elementi
 *  uguali (a meno di epsilon) nella partizione 2.
 * La partizione non ha l'elenco degli atomi, ne' il vettore prev_site, ma basta per
 *  calcolare il prodotto (una funzione solo dei labels e dell'entropia).
 */
void general_partition::reduce(const general_partition &p1, const general_partition &p2) {
    //inizializzazioni
    const int & epsilon = opts.epsilon;
    int fattori_indipendenti = 0;
    if (p1.n == 0 || p2.n == 0) {
        fprintf(stderr, "REDUCE: Trying to use empty partition\n");
        exit(1);
    }
    allocate(p1.N);

    labels.assign(N, 1);

    int common_size = N;
    //per ogni atomo
    for (label_t which = 0; which < p1.n; which++) {
        // Considero l'atomo n-esimo della prima partizione
        label_t size = p1.atomi[which].size;

        // Se gli atomi sono uguali, l'intersezione delle partizioni dicotomiche
        // non e' banale => salto

        // Creazione degli iteratori, per ottenere tutti i siti in atomo1
        Iter_t ii1 = p1.begin(which);
        Iter_t end = p1.end();

        // uguaglianza "fuzzy" tra atomi, a meno di 'epsilon' siti
        // similmente, se sono 'uguali', salto l'atomo nella partizione risultante

        //nel caso di epsilon > 0, scorrere tutto un atomo per controllare l'altro.
        if (epsilon) {
            //salto se l'atomo e' sufficientemente piccolo
            bool almost_equal = 2 * size < epsilon;
            Iter_t ii = ii1;
            /* Scorro ogni atomo della partizione 2 che ha qualche sito coincidente
             *  con l'atomo della partizione 1.
             * Ottimizzazione: controlliamo la epsilon-uguaglianza solo una volta
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
                    //se la differenza delle dimensioni e' maggiore di epsilon,
                    //sicuramente non possono essere uguali!
                    almost_equal = false;
                    break;
                }
                if (is_atom_similar(ii1, p2.begin(atomo2), end, epsilon)) {
                    almost_equal = true;
                    break;
                }
                atomo2_old = atomo2;
            }
            if (almost_equal)
                //atomo1 e' saltato, perche e' stato trovato un corrispondente
                //nella partizione 2
                continue;
        } else {
            label_t atomo2 = p2.labels[*ii1];
            if (is_atom_equal(ii1, p2.begin(atomo2), end))
                continue;
        }
        // altrimenti interseca il fattore dicotomico con i precedenti
        fattori_indipendenti++;
        common_size -= size;
        entropia_shannon += size * mylog[size];

        // faccio il prodotto rapido - moltiplico le etichette di atomo1 per il numero
        for (; ii1 != end; ii1++)
            labels[*ii1] *= fattori_indipendenti + 1;
    }
    // si tiene conto dell'atomo sconnesso di background, formato dai pezzi comuni
    entropia_shannon += common_size * mylog[common_size];
    entropia_shannon = -entropia_shannon / N + mylog[N];
    entropia_topologica = mylog[fattori_indipendenti + 1];
}


/* findroot ricorsivamente trova il primo sito del cluster, la radice appunto.
 * La formulazione ricorsiva implica path-compression.  Non e' possibile fare
 *  una versione iterativa mantenendo la compression: il compilatore non puo'
 *  fare tail optimization in questo caso, findroot e' responsabile per il 30%
 *  del tempo speso dall'applicazione :(
 */
template <typename data_t> inline int findroot(int i, data_t *ptr) {
    if (ptr[i] < 0) return i;
    return ptr[i] = findroot(ptr[i], ptr);
}

template <typename data_t>
void general_partition::from_configuration(const data_t *configuration, const adj_struct & adj, int N1) {
    label_t s1, s2;
    label_t r1, r2;
    int z;

    if(N1==0) N1=adj.N;
    allocate(N1);

    //set all labels to -1
    labels.assign(N, -1);

    // percolazione e primi labels
    for (s1 = 0; s1 < N; s1++) {
        r1 = findroot(s1, &labels[0]);

        z = adj.fetch(s1);
        for (int j = 0; j < z; j++) {
            //select next neighbor
            s2 = adj.vicini[j];

            //check them for being in the same cluster, skip when they're not
            if (s2 >= s1 || s2 < 0 || configuration[s1] != configuration[s2])
                continue;

            r2 = findroot(s2, &labels[0]);
            // attribution to proper tree root
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
template void general_partition::from_configuration(const int *configuration, const adj_struct & adj, int N1);
template void general_partition::from_configuration(const char *configuration, const adj_struct & adj, int N1);

/* relabel: funzione essenziale per il corretto partizionamento
 * Costruisce l'elenco degli atomi, il vettore prev_site;
 * A partire dai labels, che indicano le posizioni, i siti sono rinumerati,
 *  usando come etichette l'indice dall'atomo, contando dal basso.
 */
void general_partition::relabel() {
    entropia_shannon = 0;
    n = 0;

    // 0-conto nr. atomi diversi per il solo scopo di allocare ottimalmente
    for(label_t i = 0; i < N; i++){
        n += labels[i] < 0;
        //a questo punto prev_site contiene il primo sito appartenente al cluster
        prev_site[i] = findroot(i, &labels[0]);
    }
    atomi.reserve(n);
    n = 0;

    // 1-creazione array atomi
    // 2-inizializzazione ogni elemento
    // 3-creazione indice (label atomo) <--> root
    // 4-calcolo entropia a partire dai size nei root
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

    // 1-relabeling secondo l'indice dell'atomo, non del sito di appartenenza
    // 2-creazione del collegamento prev_site(che usa il precedente, non il primo com'era prima)
    //   e aggiornamento atom.end, che indica l'ultimo sito trovato apparentenente all'atomo
    for (label_t i = 0; i < N; i++) {
        int atom_label = labels[prev_site[i]];
        labels[i] = atom_label;
        prev_site[i] = std::min(atomi[atom_label].end, i);
        atomi[atom_label].end = i;
    }
}

/* linear_intersection: funzione come from_configuration (percolazione) specializzata
 *  per il caso di una struttura 2dimensionale, con argomenti due partizioni (invece di
 *  una configurazione e una struttura adiacenza!).
 */
void general_partition::linear_intersection(const general_partition &p1, const general_partition &p2) {
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
        if (s1 == s2 || s2 < 0)
            continue;
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

        //usa come secondo vicino il sito dalla partizione 2
        s2 = p2.prev_site[s1];
        if (s1 == s2 || s2 < 0)
            continue;
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

/* Funzione utile per costruire il laplaciano dei clusters.
 * La matrice di adiacenza e' cosi costruita:
 * A(i,j) = { 1 se i e' nello stesso cluster di j, 0 altrimenti}
 * Il risultato e' una matrice poco sparsa, ma puo' dare info spettrali sui clusters.
 */
void general_partition::print_cluster_adjacency() {
    //label_t e' definito in strutture.h come int32
    std::vector<label_t> riga;
    std::vector<label_t> colonna;
    int totale = 0;
    FILE *vec1 = fopen("vector1.bin", "wb");
    FILE *vec2 = fopen("vector2.bin", "wb");
    //per ogni sito
    for (label_t which = 0; which < N; which++) {
        /* Per ogni sito appartenente al reticolo, recupero le informazioni
           sul cluster di appartenenza (atomo in questa nomenclatura).
           Cio' e' necessario per ottenere tutti i siti (ordinati) del cluster
           cercato con efficienza massima. */
        const atom &atomo = atomi[labels[which]];
        int quanti = atomo.size;
        riga.reserve(quanti);
        colonna.reserve(quanti);

        int sito = 0;
        //scorro tutti i siti appartenenti allo stesso atomo,
        //con l'iteratore ii. Il sito corrispondente e'  *ii
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

        /* Creazione vettori di adiacenza per matrice sparse in stile Matlab */
        /* La funzione per caricare i dati cosi creati e':
        -------------------------------
        function adiacenza=load_sierpinski()
            indici_riga=fread(fopen('vector1.bin','r'),inf,'int32');
            indici_colonna=fread(fopen('vector2.bin','r'),inf,'int32');
            N=max(max(indici_riga),max(indici_colonna));
            adiacenza=sparse(indici_riga,indici_colonna,1,N,N);
        end
         ******************************/
    }
    fclose(vec1);
    fclose(vec2);
    printf("Elementi nonnulli della matrice di adiacenza: %d\n", totale);
}

#define ATOM atomi[label_count-1]
#define LAST_POS product[i-1].second
#define THIS_POS product[i].second
/* Calcolo della partizione prodotto a partire da due altre.
 * Il risultato e' un oggetto partizione completo, con tutte le proprieta' definite.
 * Per il calcolo di una distanza e' eccessivo - meglio usare la classe distance,
 *  che implementa lo stesso codice troncato e senza temporanei inutili.
 */
void general_partition::product(const general_partition & p1, const general_partition & p2) {
    label_t label_count = 0;
    double H = 0;
    label_t mu;
    label_t begin;
    allocate(p1.N);

    /* I label del prodotto sono rappresentati dalle coppie (pair<label_t, label_t>)
     *  di label dei fattori. Il vettore product contiene la coppia <label_prodotto, indice>.
     * Per riconoscere gli atomi, ordino il vettore product rispetto al label del prodotto:
     * -il risultato sara un vettore ordinato, per cui riconosciamo gli atomi come elementi
     *  contigui con lo stesso label. Il label e' temporaneo, gli atomi avranno indice
     *  progressivo.
     *
     * product[i].first - indice (temp) del prodotto.
     * product[i].second, da' l'indice corrispondente del sito nella partizione,
     *  e' utilizzato per indicizzare i cambiamenti e gli accessi ai vettori labels[], ecc.
     */
    vector<pair<pair<label_t, label_t>, label_t> > product(N);

    for (label_t i = 0; i < N; i++)
        product[i] = make_pair(make_pair(p1.labels[i], p2.labels[i]), i);

    sort(product.begin(), product.end());

    //the first position always starts an atom
    begin = 0;
    label_count = 1;
    pair<label_t, label_t> old_val = product[0].first;

    //we skip over the first value, so initialize here
    labels[product[0].second] = label_count;
    prev_site[product[0].second] = product[0].second;

    for (label_t i = 1; i < N; i++) {
        //whenever we find a new atom
        if (product[i].first != old_val) {
            //a new atom starts

            //the closed (old)atom's length is calculated
            mu = i - begin;
            ATOM.start = product[begin].second;
            ATOM.size = mu;
            ATOM.end = LAST_POS;
            //the new one is ready to go
            prev_site[THIS_POS] = THIS_POS;
            label_count++;
            begin = i;
            //cache the new label to check
            old_val = product[i].first;

            //we add the entropy, with the trick mu>0 and when mu=1 the log is 0
            if (mu > 1)
                H += (double) mu * mylog[mu];
        } else
            //prev site del sito corrente e' il precedente trovato nello stesso atomo
            prev_site[THIS_POS] = LAST_POS;

        labels[THIS_POS] = label_count;
    }
    //the last one, so it's not left hanging
    mu = N - begin;
    ATOM.start = product[begin].second;
    ATOM.size = mu;
    ATOM.end = product[N - 1].second;
    H += mu * mylog[mu];

    //normalize the result
    H = -H / N + mylog[N];

    entropia_topologica = mylog[label_count];
    n = label_count;
}
