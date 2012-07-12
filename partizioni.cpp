#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

#include "strutture.h"
#include "adj_handler.h"

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


void general_partition::allocate(label_t len) {
    assert(len != 0);
    
    if (!labels.empty() && N != len) {
        fprintf(stderr, "Allocating again for a different length without freeing first\n");
        exit(1);
    }
    N = len;
    
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

void general_partition::sort_entropy() {
    int label_count = 0;
    double H = 0;
    int mu;
    int begin;

    std::vector<label_t> temp(labels);
    std::sort(temp.begin(),temp.end());
    
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

bool symmetric_difference(Iter_t from1, Iter_t from2, Iter_t to, int tol = 10) {
    int differenze = 0;

    while (true) {
        // Se il numero delle differenze e' elevato, escludo l'uguaglianza
        if (differenze > tol) {
            return (false);
            // Se ho raggiunto la fine di entrambi gli atomi senza troppe differenze
            // allora accetto l'uguaglianza
        } else if ((from1 == to) && (from2 == to)) {
            return (true);
            // Nel caso di atomi con disparita' numerica, continuo a scorrere il piu
            // lungo, per vedere se accumulo abbastanza differenze
        } else if ((from1 == to) != (from2 == to)) {
            differenze++;
            from1++;
            from2++;
            // Se ho elementi diversi nei due atomi, scorro per riallinearli
        } else if (*from1>*from2) {
            differenze++;
            from1++;
        } else if (*from2>*from1) {
            differenze++;
            from2++;
            // I due atomi hanno lo stesso elemento - scorro in avanti entrambi
        } else {
            from1++;
            from2++;
        }
    }
}

void general_partition::reduce(const general_partition &p1, const general_partition &p2) {
    //inizializzazioni
    int fattori_indipendenti = 0;
    allocate(p1.N);

    for (label_t i = 0; i < N; i++)
        labels[i] = 1;

    int common_size = N;
    //per ogni atomo
    for (label_t which = 0; which < p1.n; which++) {
        // Considero l'atomo n-esimo del primo
        // Trovo l'atomo che corrisponde nella seconda partizione
        //  attraverso il primo sito in comune
        const atom &atomo1 = p1.atomi[which];
        const atom &atomo2 = p2.find_atom(atomo1);
        int size = atomo1.size;

        // Se gli atomi sono uguali, l'intersezione delle partizioni dicotomiche
        // non e' banale => salto

        // Creazione degli iteratori, per ottenere tutti i siti in atomo1
        Iter_t ii = p1.begin(atomo1);
        Iter_t end = p1.end();

        // uguaglianza "fuzzy" tra atomi, a meno di 'tol' siti
        // similmente, se sono 'uguali', salto l'atomo nella partizione risultante
        
        //nel caso di epsilon > 0, scorrere tutto un atomo per controllare l'altro... peccato!
        if (symmetric_difference(ii, p2.begin(atomo2), end, 0))
            continue;

        // altrimenti interseca il fattore dicotomico con i precedenti
        fattori_indipendenti++;
        common_size -= size;
        entropia_shannon += size * mylog[size];

        // faccio il prodotto rapido - moltiplico i siti di atomo1 per numero
        for (; ii != end; ii++)
            labels[*ii] *= fattori_indipendenti + 1;
    }
    // si tiene conto dell'atomo sconnesso di background, formato dai pezzi comuni
    entropia_shannon += common_size * mylog[common_size];
    entropia_shannon = -entropia_shannon / N + mylog[N];
    entropia_topologica = mylog[fattori_indipendenti + 1];
    //printf("Stima: %g, sicura: %g\n",entropia,entropia_shannon);
}

template <typename data_t> int findroot(int i, data_t *ptr) {
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

void general_partition::relabel() {
    static std::vector<label_t> new_label;
    new_label.reserve(N);
    entropia_shannon = 0;
    n = 0;

    // 0-conto nr. atomi diversi per il solo scopo di allocare ottimalmente
    for(label_t i = 0; i < N; i++)
        n += labels[i] < 0;
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
            new_label[i] = n;
            prev_site[i] = i;
            n++;
        } else {
            prev_site[i] = findroot(i, &labels[0]);
            new_label[i] = prev_site[i];
        }
    }

    entropia_topologica = mylog[n];
    entropia_shannon = -entropia_shannon / N + mylog[N];
    // 1-relabeling secondo l'indice dell'atomo, non del sito di appartenenza
    // 2-hashing
    // 3-creazione del collegamento prev_site e atom.end
    for (label_t i = 0; i < N; i++) {
        int atom_pos = new_label[prev_site[i]];
        labels[i] = atom_pos;
        prev_site[i] = std::min(atomi[atom_pos].end, i);
        atomi[atom_pos].end = i;
    }
}

void general_partition::linear_intersection(const general_partition &p1, const general_partition &p2) {
    label_t s1, s2;
    allocate(p1.N);

    //set all labels to -1
    labels.assign(N, -1);

    // percolazione e primi labels
    for (s1 = 0; s1 < N; s1++) {
        int r2;
        int r1 = findroot(s1, &labels[0]);

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