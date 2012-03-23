/* 
 * File:   strutture.h
 * Author: fake
 *
 * Created on January 13, 2012, 10:51 PM
 */

#include "rand55.h"


#define FORCE_INLINE __attribute__((always_inline))

enum algorithm {
    NORMAL,
    SORTED,
    PMATRIX,
    MAP,
    AUTO
};

enum DISTANZE {
    SHAN = 1,
    SHAN_TOP = 1 << 1,
    RID = 1 << 2,
    RID_TOP = 1 << 3,
    GENERAL = 1 << 4,
    GENERAL_TOP = 1 << 5,
    GENERAL_RID = 1 << 6,
    GENERAL_RID_TOP = 1 << 7,
    HAMM = 1 << 20
};

enum source{
    SEQUENCE=1,
    LATTICE=1<<1,
    RANDOM=1<<2,
    FROM_FILE=1<<3
};


typedef struct {
    int seq_len;
    int n_seq;
    int lato;
    
    int from;
    char filename[255];
    bool write;
    bool distance;

    int fuzzy;
    int seed;
    int verbose;
    int n_symbols;
    bool translate;
    bool graphics;
    
    int threads;
    algorithm alg;
    
    int da_calcolare;

}options;


typedef struct {
    int value;
    int pos;
} labelled_array;

class basic_partition{
public:
    //number of atoms found
    int n;
    //total length of the partition
    int N;
    
    double entropia_shannon;
    double entropia_topologica;
    
    
};

class linear_partition: public basic_partition {

public:
    //array of atom starting positions {1,4,5,...}
    //int *atom_positions;
    //binary atom beginning codes {1,0,0,0,1,1,...}
    int *binary;
        
    template <typename T>  linear_partition(const T*seq, int len);
    template <typename T>  void fill(const T*seq, int len);
    void print();
    linear_partition(int len=0){
        n=0;
        N=len;
        binary=new int[len];
        //atom_positions=new int[len];
    }
    ~linear_partition(){
        if(N){
        //delete[] atom_positions;
        delete[] binary;
        }
        n=N=0;
    }
};

class atom {
public:
    int size;    
    int start;
    int end;    
    u_int32_t hash;
    
    atom(){
        size=0;        
        end=0;
        start=0;
        hash=5381;
    }

    bool operator==(const atom &due) const {
        return (hash==due.hash);
    }
    
    bool operator!=(const atom &due) const {
        return (hash!=due.hash);
    }
};

class general_partition: public basic_partition{
public:
    //labels identify generic atoms across the partition
    int *labels;
    //nearest neighbors
    int *prev_site;
    
    
    atom * atomi;
    int lato;
    
    int dim;
    int **NNB;
    
    void allocate(int len);
    template <typename T> general_partition(const T* seq, int len);
    template <typename T> void fill(const T* seq, int len);
    template <typename T> void from_linear_sequence(const T* seq, int len);
    template <typename T> void from_square_lattice(const T* val, int lato,int dim);
    void trivial(int len);
    void from_nnb(int **NNB, int dim);
    
    void from_atom(int *label, const int which, const int set_to);
    general_partition(int len=0);
    void sort_entropy();
    void reduce(const general_partition &ridurre, const general_partition &common);
    void linear_intersection(const general_partition &p1, const general_partition &p2);
    void print();    
    atom& find_atom(const atom &atomo1) const{ return atomi[labels[atomo1.start]];}
    
    
    ~general_partition();
    
    class Iterator{
    private:
        int _site;
        const int *_next;
        
    public:
        Iterator(int dove, const int *vicini): _site(dove), _next(vicini) {};
        
        int operator*(){
            return(_site);
        }
        
        bool operator==(Iterator due){
            return (_site == due._site);
        }
        
        bool operator!=(Iterator due){
            return (_site != due._site);
        }
        
        int operator++(){
            if(_site==_next[_site])
                _site=-1;
            else
                _site=_next[_site];
            return(_site);
        }
                
        int operator++(int) { return(operator++());}        
        
    };
    
    Iterator begin(const int where) const{
        return Iterator(atomi[where].end,prev_site);
    }
    
    Iterator begin(const atom &where) const{
        return Iterator(where.end,prev_site);
    }
    
    Iterator end() const{
        return Iterator(-1,prev_site);
    }
};

typedef general_partition::Iterator Iter_t;


class distance{
private:
    void allocate(int n);

    int *common_factor;
    int *reduced1;
    int *reduced2;
    u_int64_t *product;
    int *product_reduced;
    int *labels;
    short *pmatrix;
    general_partition ridotto1;
    general_partition ridotto2;
    general_partition partizione_comune;
    //std::unordered_map<int, int> mappa;
    
    int N;
   
    
public:
    double dist_s;
    double dist_s_r;
    double dist_top;
    double dist_top_r;
    double dist_ham;
    double dist_fuzzy;
    double dist_fuzzy_t;
    double dist_fuzzy_r;
    double dist_fuzzy_r_t;
    
    
        
    void binary_partition(const linear_partition &first, const linear_partition &second,int ridotta=1);
    void linear_product_sorted(const general_partition &p1, const general_partition &p2);
    void fill(const linear_partition& e1, const linear_partition& e2) ;
    void fill(const general_partition& e1, const general_partition& e2) ;
    template <typename T> void hamming_distance(const T* seq1,const T* seq2);
    
    ~distance();
    distance(int n=1000);
    distance(const distance &d1);
    

};

void set_program_options(options &opt, int argc, char**argv) ;
template <typename T> double entropy_binary_partition(const T *p, int n);
void fill_entries_randomly(const options opts, std::string *entries);
void fill_seq_from_file(options &opts, std::string* entries);
void generate_next_sequence(std::string &entry);
void load_lattices_from_file(options &opts, int** num_entries);
template <typename T> void ppmout(const T *grid, int sz, const char *filename);
template <typename T, typename U> void ppmout2(const T *grid1, const U* grid2, int sz, const char *filename);

void calcola_matrice_distanze(linear_partition *X, general_partition *Z, std::string *char_entries);
void print_array(int *array, int len, const char *nome);