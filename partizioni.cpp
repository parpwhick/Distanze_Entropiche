#include <cstdlib>
#include <cstdio>
#include <string>
#include <cassert>

#include "strutture.h"


extern options opts;
extern double *mylog;
extern int *colore;


void print_array(int *array, int len, const char *nome) {
    printf("%s [%3d", nome, array[0]);
    for (int i = 1; i < len; i++)
        printf(",%3d", array[i]);
    printf("]\n");
}

template <typename T>
void print_square_lattice(const T* valori, int lato){
    for(int i=0; i < lato; i++){
        for(int j=0; j < lato; j++){
            printf("%2d ",valori[j*lato+i]);
        }
        printf("\n");
    }
    printf("\n");
}


#define COL_MAX 1000
template <typename T>
void ppmout(const T *grid1, int sz, const char *filename) {

    if(colore==0){
        colore = new int[COL_MAX];
        for (int i = 0; i < COL_MAX; i++) {
            colore[i] = xrand();
        }
        colore[0] = 0x0F5A3A1F; // blue almost black
    }
    
    int MULT;
    if(sz > 500)
        MULT=1;
    else
        MULT=500/sz + 1;
    
    
    FILE *fout = fopen(filename,"w");
    fprintf(fout, "P6\n %d %d\n 255\n", MULT*sz, MULT*sz);


    for (int rg = 0; rg < sz; rg++) {
        for (int i = 0; i < MULT; i++)
            for (int cl = 0; cl < sz; cl++) {
                int sito=grid1[rg * sz + cl] % COL_MAX;
                int color = colore[sito];
                for (int j = 0; j < MULT; j++)
                    fwrite(&color, 3, 1, fout);
            }
    }
}
template <typename T, typename U>
void ppmout2(const T *grid1, const U* grid2, int sz, const char *filename) {
    //Se l'array dei colori non e' inizializzato - riempiamolo!
    if(colore==0){
        colore = new int[COL_MAX];
        for (int i = 0; i < COL_MAX; i++) {
            colore[i] = xrand();
        }
        colore[0] = 0x0F5A3A1F; // blue almost black
    }
    
    int MULT;
    if(sz > 500)
        MULT=3;
    else
        MULT=500/sz + 1;
    
    
    int black=0x20202020;
    int intermezzo=100; //pixels tra i pannelli
    
    FILE *fout = fopen(filename,"w");
    fprintf(fout, "P6\n %d %d\n 255\n", 2*MULT*sz+intermezzo, MULT*sz);

    for (int rg = 0; rg < sz; rg++) {
        for (int i = 0; i < MULT; i++){
            //primo reticolo
            for (int cl = 0; cl < sz; cl++) {
                int sito=grid1[rg * sz + cl] % COL_MAX;
                int color = colore[sito];
                for (int j = 0; j < MULT; j++)
                    fwrite(&color, 3, 1, fout);
            }
            //10 pixels neri in mezzo
            for (int j = 0; j < intermezzo; j++)
                    fwrite(&black, 3, 1, fout);
            
            //secondo reticolo
            for (int cl = 0; cl < sz; cl++) {
                int sito=grid2[rg * sz + cl] % COL_MAX;
                int color = colore[sito];
                for (int j = 0; j < MULT; j++)
                    fwrite(&color, 3, 1, fout);
            }
        }
    }
}

template void ppmout(const int32_t *grid, int sz, const char *filename);
template void ppmout(const uint64_t *grid, int sz, const char *filename);
template void ppmout2(const int *grid1, const int* grid2, int , const char *) ;
template void ppmout2(const uint64_t *grid1, const unsigned long* grid2, int , const char *) ;


/************************************************
 ************************************************
 *             PARTIZIONI SEMPLICI              *
 ************************************************
 ************************************************
 */


void linear_partition::print() {
   if (opts.verbose > 3) {
        // {1,0,0,1,0,1,...}
       print_array(binary,N,"Binary");
    }
    fprintf(stdout, "Partitions[n]: %d, Shannon %f, Topological %f\n", n,entropia_shannon,
            entropia_topologica);
}

template <typename T>
linear_partition::linear_partition(const T* seq, int len) {
    this->fill(seq,len);
}

template <typename T>
void linear_partition::fill(const T* seq, int len) {
    int i, j;
    
    //Total length of the partition is equal to the sequence
    N = len;
        
    binary = new int[len];
    //atom_positions = new int[len];

    //first one always start an atom
    binary[0] = 1;

    //checking for a different symbol from the one before
    //when that happens => new atom!
    for (i = 1; i < len; i++)
        binary[i] = seq[i] != seq[i - 1];

    //now going back and writing the atom positions
    j=0;
    for (i = 0; i < len; i++)
        if (binary[i]) {
            //atom_positions[j++] = i + 1;
            j++;
        }
    //number of atoms found - index++ of last the atom in the array (zero addressing)
    n = j;
    entropia_topologica=mylog[n];
    entropia_shannon=entropy_binary_partition(binary,N);
}
template void linear_partition::fill(const char *, int);


/************************************************
 ************************************************
 *            PARTIZIONI COMPLESSE              *
 ************************************************
 ************************************************
 */


void general_partition::allocate(int len){
    assert(len !=0);
    N=len;
    
    if (!labels)
        labels = new int[N];
    if (!prev_site)
        prev_site = new int[N];
    if (!atomi)
        atomi = new atom[N];
}

general_partition::general_partition(int len) {
    n = 0;
    N = len;
    entropia_topologica = 0;
    entropia_shannon = 0;

    if (N)
        allocate(N);
    else {
        labels = 0;
        prev_site = 0;
        atomi = 0;
        NNB = 0;
    }    
}

general_partition::~general_partition(){
    if (N) {
        delete []labels;
        delete []prev_site;
        delete []atomi;
    }    
    if(NNB)
        delete []NNB;
}


template <typename T>
general_partition::general_partition(const T* seq, int len) {
    labels=0;
    N=0;
    prev_site=0;
    this->from_linear_sequence(seq,len);
}

void general_partition::trivial(int len){
    N=len;
    int *vuoto=new int[N];
    for(int i=0; i< N; i++)
        vuoto[i]=0;
    
    from_linear_sequence(vuoto,N);
    delete []vuoto;
    
}

template <typename T>
void general_partition::from_linear_sequence(const T* seq, int len) {
    int i, j;
    int last_good;
    
    allocate(len);    
    dim=1;
    
    //presetting labels for this partition to 0
    for (i=0;i<N;i++)
        prev_site[i]=-1;
     
    for (i=0;i<N;i++){
        //if the site was already "colored" check the next one
        if (prev_site[i] != -1)
            continue;     
        last_good=i;
        prev_site[i]=i;
        
        for(j=i+1; j<last_good+opts.fuzzy+2 && j<N; j++)
            //if it belongs to the same cluster...
            if(seq[j]==seq[i]){                
                prev_site[j]=last_good;
                last_good=j;
            }
    }
    NNB=new int*[dim];
    NNB[0]=prev_site;
    
    from_nnb(NNB,dim); 
    
    if (opts.graphics && (opts.topologia & LINEARE)){
        static int imagenr=0; 
        char filename[255];
        imagenr++;
        sprintf(filename, "sequenza%03d.ppm", imagenr);
        ppmout(labels, N, filename);
    }
}
template void general_partition::from_linear_sequence(const char *,int);
template void general_partition::from_linear_sequence(const int *,int);


#define nnu (i - (i % lato)+ ((i+lato-1)%lato))
#define nnd ((i/lato)*lato + ((i+lato+1)%lato))
#define nnl (i+N-lato)%N
#define nnr (i+N+lato)%N

template <typename T>
void general_partition::from_square_lattice(const T* valori, int L,int) {    
    lato=L;
    N=L*L;
    dim=2;
    allocate(N);

    NNB=new int*[dim];
    for(int j=0;j<dim;j++)
        NNB[j]=new int[N];

    // Ho due tipi di vicini, quello Up e quello Left
    for (int i = 0; i < N; i++) {
        NNB[0][i] = (valori[i] == valori[nnu]) ? nnu : i;
        NNB[1][i] = (valori[i] == valori[nnl]) ? nnl : i;
    }  
    //Calcola la partizione, a partire dalla mappa dei vicini
    from_nnb(NNB, dim);

    if (opts.graphics) {
        static int imagenr = 0;
        char filename[255];
        imagenr++;
        sprintf(filename, "reticolo%03d.ppm", imagenr);
        ppmout2(valori, labels, L, filename);
    }

}
template void general_partition::from_square_lattice(const int*, int, int);
template void general_partition::from_square_lattice(const char*, int, int);



inline int compare (const void * a, const void * b){
  return ( *(int*)a - *(int*)b );
}

void general_partition::sort_entropy(){
    int i;
    int label_count = 0;
    double H = 0;
    int mu;
    int begin;
    
    int *temp=new int32_t[N];
    allocate(N);
    
    for(int i=0; i<N;i++)
        temp[i]=labels[i];
    
    qsort(temp,N,sizeof(temp[0]),compare);
    
    //the first position always starts an atom
    begin = 0;
    label_count=1;
    int old_val=temp[0];
        
    for (i = 1; i < N; i++) {      
        //whenever we find a new atom
        if (temp[i]!=old_val) {
            //a new atom starts
            
            //the closed (old)atom's length is calculated
            mu = i - begin;           
            //the new one is ready to go
            label_count++;
            begin = i;
            //cache the new label to check
            old_val=temp[i];
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
    
    entropia_topologica=mylog[label_count];
    n=label_count;
    entropia_shannon=H;
    
    delete []temp;
}

bool symmetric_difference(Iter_t from1, Iter_t from2, Iter_t to, int tol=10){               
    int differenze=0;
    
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
        }else if ((from1 == to) != (from2 == to)) {
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

void general_partition::reduce(const general_partition &p1, const general_partition &p2){
    //inizializzazioni
    int fattori_indipendenti=0;
    lato=p1.lato;
    allocate(p1.N);
  
    for(int i=0;i<N;i++)
        labels[i]=1;
    
    //per ogni atomo
    for (int which = 0; which < p1.n; which++) {
        // Considero l'atomo n-esimo del primo
        // Trovo l'atomo che corrisponde nella seconda partizione
        //  attraverso il primo sito in comune
        const atom &atomo1 = p1.atomi[which];
        const atom &atomo2 = p2.find_atom(atomo1);          
        
        // Se gli atomi sono uguali, l'intersezione delle partizioni dicotomiche
        // non e' banale => salto
        
        // Creazione degli iteratori, per ottenere tutti i siti in atomo1
        Iter_t ii = p1.begin(atomo1);
        Iter_t end = p1.end();       
        
        // uguaglianza "fuzzy" tra atomi, a meno di 'tol' siti
        // similmente, se sono 'uguali', salto l'atomo nella partizione risultante
        if(symmetric_difference(ii,p2.begin(atomo2),end,0)     )
            continue;
        
        // altrimenti interseca il fattore dicotomico con i precedenti
        fattori_indipendenti++;               
        
        // faccio il prodotto rapido - moltiplico i siti di atomo1 per numero
        for (; ii != end; ii++) 
            labels[*ii] *= fattori_indipendenti + 1;           
    }   
    this->sort_entropy();
}



int findroot(int i,int *ptr)
{
  if (ptr[i]<0) return i;
  return ptr[i] = findroot(ptr[i],ptr);
}

class neigh_factory{
public:
    int n;    
    int buffer[1000];
    int *vicinato[2];
    int *configuration;
    int **adiacenza;
    int siti;
    int z;

    void f1(int site){
        n=2;
        buffer[0]=vicinato[0][site];
        buffer[1]=vicinato[1][site];
    }
    
    void f2(int site){
        n=0;
        for (int i=0; i<z; i++){
            int s1=adiacenza[site][i];
            if(s1 < 0 || s1==site || configuration[site]!=configuration[s1])
                continue;
            buffer[n++]=s1;
        }            
    }
    void (neigh_factory::*fetch)(int);
    
    int operator[](int i){
        return buffer[i];
    }
    
    void operator()(int site){
        (this->*fetch)(site);
    }
    
    void init(const general_partition &p1, const general_partition &p2) {
        n = 2;
        vicinato[0] = p1.prev_site;
        vicinato[1] = p2.prev_site;
        fetch=&neigh_factory::f1;
    }
    
    void init(int *valori_siti, int **adj, int N, int nmax) {
        configuration=valori_siti;
        adiacenza=adj;
        siti=N;
        z=nmax;
        fetch=&neigh_factory::f2;
    }
};


void general_partition::linear_intersection(const general_partition &p1, const general_partition &p2){
    lato=p1.lato;
    allocate(p1.N);
    
    neigh_factory get_neigh;
    get_neigh.init(p1,p2);
    // Calcolo a partire dall'insieme dei vicini    
    from_nnb(get_neigh); 

    // Grafico il reticolo risultante, se richiesto
    if (opts.graphics && (opts.topologia & RETICOLO)) {
        static int imagecount=0;     
        char filename[255];
        imagecount++;
        sprintf(filename, "comune%03d.ppm", imagecount);
        ppmout2(p1.labels, p2.labels, lato, filename);
        imagecount++;
        sprintf(filename, "comune%03d.ppm", imagecount);
        ppmout(labels, lato, filename);
    }
}


void general_partition::from_nnb(neigh_factory get_neigh){
    int s1, s2;
    int i;
    
    for (i = 0; i < N; i++) {
        labels[i] = -1;
    }  
 
    // percolazione e primi labels
    for (s1 = 0; s1 < N; s1++) {
        int r2;
        int r1=findroot(s1,labels);

        get_neigh(s1);
        for (int j = 0; j < get_neigh.n; j++) {            
            s2 = get_neigh[j];
            r2 = findroot(s2, labels);
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
    int *new_label=new int[N];
    entropia_shannon=0;
    n=0;
    
    // 1-creazione array atomi
    // 2-inizializzazione ogni elemento
    // 3-creazione indice (label atomo) <--> root
    // 4-calcolo entropia a partire dai size nei root
    #define ATOMO atomi[n]
    for(i=0;i<N;i++){
        if(labels[i]<0){
            entropia_shannon+= -labels[i]*mylog[-labels[i]];
            ATOMO.size=-labels[i];
            ATOMO.end=i;
            ATOMO.start=i;
            new_label[i]=n;
            prev_site[i]=i;
            n++;
        } else{
            prev_site[i]=findroot(i,labels);
            new_label[i]=prev_site[i];
        }
    }
 
    entropia_topologica=mylog[n];
    entropia_shannon= -entropia_shannon/N+mylog[N];    
    // 1-relabeling secondo l'indice dell'atomo, non del sito di appartenenza
    // 2-hashing
    // 3-creazione del collegamento prev_site e atom.end
    for(i=0;i<N;i++){
        int atom_pos=new_label[prev_site[i]];
        labels[i]=atom_pos;
        prev_site[i]=std::min(atomi[atom_pos].end,i);
        atomi[atom_pos].end=i;      
    }
    
    delete []new_label;        
}


void general_partition::from_nnb(int **neighbors, int dimensione){
    int s1, s2;
    int i;
    dim=dimensione;
    
    for (i = 0; i < N; i++) {
        labels[i] = -1;
    }  
 
    // percolazione e primi labels
    for (s1 = 0; s1 < N; s1++) {
        int r2;
        int r1=findroot(s1,labels);

        for (int j = 0; j < dim; j++) {            
            s2 = neighbors[j][s1];
            if (s1 == s2 || s2 < 0)
                continue;
            r2 = findroot(s2, labels);
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
    int *new_label=new int[N];
    entropia_shannon=0;
    n=0;
    
    // 1-creazione array atomi
    // 2-inizializzazione ogni elemento
    // 3-creazione indice (label atomo) <--> root
    // 4-calcolo entropia a partire dai size nei root
    #define ATOMO atomi[n]
    for(i=0;i<N;i++){
        if(labels[i]<0){
            entropia_shannon+= -labels[i]*mylog[-labels[i]];
            ATOMO.size=-labels[i];
            ATOMO.end=i;
            ATOMO.start=i;
            new_label[i]=n;
            prev_site[i]=i;
            n++;
        } else{
            prev_site[i]=findroot(i,labels);
            new_label[i]=prev_site[i];
        }
    }
 
    entropia_topologica=mylog[n];
    entropia_shannon= -entropia_shannon/N+mylog[N];    
    // 1-relabeling secondo l'indice dell'atomo, non del sito di appartenenza
    // 2-hashing
    // 3-creazione del collegamento prev_site e atom.end
    for(i=0;i<N;i++){
        int atom_pos=new_label[prev_site[i]];
        labels[i]=atom_pos;
        prev_site[i]=std::min(atomi[atom_pos].end,i);
        atomi[atom_pos].end=i;      
    }
    
    delete []new_label;        
}
