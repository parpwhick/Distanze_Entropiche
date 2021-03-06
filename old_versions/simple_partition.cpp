#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <string>
#include "strutture.h"

extern options opts;
extern double *mylog;


void print_array(int *array, int len, const char *nome);


#define	FORCE_INLINE __attribute__((always_inline))


void DJBHASH_STEP(u_int32_t &hash, u_int32_t value) ;


template <class T> class my_allocator;

// specialize for void:
template <> class my_allocator<void> {
public:
    typedef void*       pointer;
    typedef const void* const_pointer;
    // reference to void members are impossible.
    typedef void value_type;
    template <class U> struct rebind { typedef my_allocator<U>    other; };
};

template <typename T> class my_allocator : public std::allocator<T> {
public:
    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;
    typedef T*        pointer;
    typedef const T*  const_pointer;
    typedef T&        reference;
    typedef const T&  const_reference;
    typedef T         value_type;
    static size_type total_memory;

    template <class U> 
    struct rebind { 
        typedef my_allocator<U> other; 
    };

    my_allocator() throw() 
    {
        
    }

    my_allocator(const my_allocator& to_copy) throw() 
    { 
    }

    template <class U> 
    my_allocator(const my_allocator<U>& to_copy) throw()
    {
    }

    ~my_allocator() throw()
    {
    }

    pointer address(reference x) const
    {
        return std::allocator<T>::address(x);
    }

    const_pointer address(const_reference x) const
    {
        return std::allocator<T>::address(x);
    }

    pointer allocate(size_type s1, typename std::allocator<void>::const_pointer hint = 0)
    {
        size_t block_size = s1 * sizeof (T);
        total_memory += block_size;
        std::cout << typeid(T).name() << " allocated, MiB: " <<  (block_size / (1024.0 * 1024.0)) << ", total: " << (total_memory>>20) << "\n";
        
        return std::allocator<T>::allocate(s1, hint);
    }

    void deallocate(pointer p, size_type n)
    {
        size_t block_size = n * sizeof (T);
        total_memory -= block_size;
        std::cout << typeid(T).name() << " deallocated, MiB: " <<   (block_size/ (1024.0 * 1024.0)) << ", total: " << (total_memory>>20) << "\n";
        std::allocator<T>::deallocate(p, n);
    }

    size_type max_size() const throw()
    {
        return std::allocator<T>::max_size();
    }

    void construct(pointer p, const T& val)
    {
        std::allocator<T>::construct(p, val);
    }

    void destroy(pointer p)
    {
        std::allocator<T>::destroy (p);
    }
};
template<typename T> std::size_t my_allocator<T>::total_memory = 0;

class comparatore_coppie{
public:
    comparatore_coppie(const general_partition &p1,const general_partition &p2) :
    l1(p1.labels.data()), l2(p2.labels.data()) {}

    bool inline operator<(const label_t &a, const label_t &b) const  {
           return (l1[a] < l1[b]) ||
                  (!(l1[a] < l1[b]) && l2[a] < l2[b]);
       }
private:
    label_t *l1;
    label_t *l2;
};

template <typename data_t> inline int findroot_recursive(int i, data_t *ptr) { //143 secondi
    if (ptr[i] < 0) return i;
    return ptr[i] = findroot_recursive(ptr[i], ptr);
}

template <typename data_t> inline int findroot_stack(int i, data_t *ptr) { // 142 secondi
    int n = 0;
    data_t stack[50];
    while (ptr[i] >= 0) {
        stack[n++] = i;
        i = ptr[i];
    }
    for (int k = 0; k < n; k++)
        ptr[stack[k]] = i;
    return i;
}


#define HASHVAR (atomi[label_count-1].hash)
template <typename T>
void linear_partition::fill(const T* seq, int len,int fuzzy) {
    int i, j;
    //this run's  atoms count
    int label_count;
    //position of the last site belonging to the atom we're trying to add to
    int last_good;
    //Total length of the partition is equal to the sequence
    N = len;
    entropia_shannon=0;
    entropia_topologica=0;
    
    if(fuzzy==-1)
        fuzzy=opts.fuzzy;
     
    if(!labels)
        allocate(N);
    
    //presetting labels for this partition to 0
    for (i=0;i<N;i++)
        labels[i]=0;
     
    
    //beginning count of identified atoms
    label_count=0;
    //we check every site for belonging to a partition set
    for (i=0;i<N;i++){
        //if the site was already "colored" check the next one
        if (labels[i])
            continue;
        //we give a new label to the site
        label_count++;
        //atom_positions[label_count-1]=i;
        labels[i]=label_count;
        
        //we know that it belongs (first, trivial)
        last_good=i;
        
        DJBHASH_STEP(HASHVAR,i);
        
        int this_atom_size=1;
        //il primo a sinistra non ha vicino
        nnL[i]=i;
        //now we add all the other possible sites to it, starting from i+1!
        //staying below N
        //and allowing a spacing of at most opts.fuzzy spaces between them
        for(j=i+1; j<last_good+fuzzy+2 && j<N; j++){
            
            //if it belongs to the same cluster...
            if(seq[j]==seq[i]){
                //..we label it accordingly
                labels[j]=label_count;
                //the last site's neighbor is the current one
                nnL[j]=last_good;
                next_site[last_good]=j;
                DJBHASH_STEP(HASHVAR,j);
                //and marking it as "good" (since we're on an ordered line)
                last_good=j;
                this_atom_size++;
                
            }
        }
        //l'ultimo a destra non ha vicino
        next_site[last_good]=last_good;
        
        //tutti i dati necessari per non cercare gli atomi
        atomi[label_count-1].size=this_atom_size;
        atomi[label_count-1].start=i;
        atomi[label_count-1].end=last_good;
        
        entropia_shannon+=this_atom_size * mylog[this_atom_size];
        
        
    }
    //this many different atoms were found
    entropia_shannon=-entropia_shannon/N+mylog[N];
    n=label_count;
    entropia_topologica=mylog[n];

}
//template linear_partition::linear_partition(const int *, int, int);
template void linear_partition::fill(const int *, int, int);
template void linear_partition::fill(const char *, int, int);



int compare_labelled_array (const void * a, const void * b) {
  return ( ((labelled_array*)a)->value - ((labelled_array*)b)->value );
}

void sort_entropy(linear_partition &p, labelled_array *product) {
    int i;
    int label_count = 0;
    double H = 0;
    int mu;
    int begin;
    
    qsort(product,p.N,sizeof(labelled_array),compare_labelled_array);
    
    //the first position always starts an atom
    begin = 0;
    label_count=1;
    int old_val=product[0].value;
        
    for (i = 1; i < p.N; i++) {      
        //whenever we find a new atom
        if (product[i].value!=old_val) {
            //a new atom starts
            
            //the closed (old)atom's length is calculated
            mu = i - begin;
           
            //the new one is ready to go
            label_count++;
            begin = i;
            //cache the new label to check
            old_val=product[i].value;

            //we add the entropy, with the trick mu>0 and when mu=1 the log is 0
            if (mu > 1)
                H += (double) mu * mylog[mu];   
        }
    }
    //the last one, so it's not left hanging
    mu = p.N - begin;
    H += mu * mylog[mu];
    
    //normalize the result
    H = -H / p.N + mylog[p.N];
    
    p.entropia_topologica=mylog[label_count];
    p.n=label_count;
    p.entropia_shannon=H;
}



//bool quickest_nontrivial_intersection(const linear_partition &ridurre, const linear_partition &common,const atom &atomo, const atom &com_atom){
//
//    int start=atomo.start;
//    int end=atomo.end;
//    int size=atomo.size;
//    bool ret=true;
//    
//    bool canary= com_atom.hash == atomo.hash;
//    
//    ret = ret && com_atom.start == start;
//    ret = ret && com_atom.end==end;
//    ret = ret && com_atom.size==size;
//   
//    
//    if(ret != canary){
//    printf("atomo_ridurre: %lx, atomo_common: %lx\n",atomo.hash,com_atom.hash);
//    printf("atom label: %d, common label:%d\n",ridurre.labels[start],common.nnR[start]);
//    printf("atom start: %d, common start: %d\n", start, common.nnL[start]);
//    printf("atom end: %d, common end: %d\n", end, com_atom.end);
//    printf("atom size: %d, common size: %d\n", size, -common.labels[common.nnL[start]]);
//    printf("return:%d, canary:%d\n",ret,canary);
//    printf("\n");  
//    }
//    return(ret);
//}


//quickest





void linear_partition::print() {
    if (opts.verbose > 1) {
        print_array(labels,N,"Lbl");

        print_array(nnL,N,"Nnb");
        
    }
    //adding the entropy information
    fprintf(stdout, "Partitions[f]: %d, Shannon %f, Topological %f\n", n, entropia_shannon,
            entropia_topologica);
    
    
}


void distance::linear_product_pmatrix(const linear_partition &p1, const linear_partition &p2) {
    int i;
    int pmax;
    int label_count=0;
    double H = 0;
    

    //limitiamo la matrice al numero minimo di righe*colonne
    pmax=(p1.n+1) * (p2.n+1);
        
    for (i=0;i<pmax;i++)
        pmatrix[i]=0;
    
    for (i = 0; i < N; i++)
        pmatrix[(p1.labels[i])*p2.n + p2.labels[i]]++;
    
    for(i=0;i<pmax;i++)
        if(pmatrix[i]){
            label_count++;
            H+=pmatrix[i]*mylog[pmatrix[i]];
        }
            
    //normalize the result
    H = -H / N + mylog[N];

    double h1 = p1.entropia_shannon,
            h2 = p2.entropia_shannon,
            t1 = p1.entropia_topologica,
            t2 = p2.entropia_topologica;
    this->dist_fuzzy = 2 * H - h1 - h2;
    this->dist_fuzzy_t = 2 * mylog[label_count] - t1 - t2;
}


//void distance::linear_product_map(const linear_partition &p1, const linear_partition &p2) {
//    int i;
//    int label_count=0;
//    double H = 0;
//
//    mappa.clear();
//        
//    for (i = 0; i < N; i++)
//        mappa[(p1.labels[i] << 16) | (p2.labels[i])]++;
//
//    for (std::unordered_map<int, int>::iterator ii = mappa.begin(); ii != mappa.end(); ++ii) {
//        label_count++;
//        H += ii->second * mylog[ii->second];
//    }
//
//    //normalize the result
//    H = -H / N + mylog[N];
//
//    double h1 = p1.entropia_shannon,
//            h2 = p2.entropia_shannon,
//            t1 = p1.entropia_topologica,
//            t2 = p2.entropia_topologica;
//    this->dist_fuzzy = 2 * H - h1 - h2;
//    this->dist_fuzzy_t = 2 * mylog[label_count] - t1 - t2;
//}



template <typename T>
void general_partition::from_linear_sequence(const T* seq, int len) {
    label_t i, j;
    label_t last_good;
    
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
    NNB=new label_t*[dim];
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
