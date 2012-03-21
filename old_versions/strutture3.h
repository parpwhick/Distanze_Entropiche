/* 
 * File:   strutture.h
 * Author: fake
 *
 * Created on January 13, 2012, 10:51 PM
 */

using namespace std;

typedef struct {
    int seq_len;
    int n_seq;
    bool random_sequences;
    char filename[255];
    bool write;
    bool distance;

    int fuzzy;
    int seed;
    int verbose;
    int n_symbols;
    bool translate;

}options;

struct dist_tuple {
    double dist_s;
    double dist_s_r;
    double dist_top;
    double dist_top_r;
    double dist_ham;
    double dist_fuzzy;
    double dist_fuzzy_t;
};

class partition {

public:
    //array of atom starting positions {1,4,5,...}
    int *atom_positions;
    //binary atom beginning codes {1,0,0,0,1,1,...}
    int *binary;
    //number of atoms found
    int n;
    //total length of the partition
    int N;
    
    double entropia_shannon();
    double entropia_topologica();
    void print();
    
    template <typename T>  partition(const T*seq, int len);
    partition(){
        n=0;
        N=0;
        binary=0;
        atom_positions=0;
    }
    ~partition(){
        delete[] atom_positions;
        delete[] binary;
        n=N=0;
    }
};

class general_partition: public partition{
private :
    partition::binary;
public :
    //labels identify generic atoms across the partition
    int *labels;
    //not to count them every time, we do it at initialization time
    int *atom_sizes;
    
    template <typename T> general_partition(const T* seq, int len);
    double entropia_shannon();
    void print();
    
    ~general_partition(){
        delete []labels;
        delete []atom_sizes;
    }
};



class entry {
public:
    //identifier of the sequence
    char* name;
    //string of symbols making the sequence
    char* seq;
    //length of the sequence
    int n;
    //sequence index, relevant only for printing
    int index;
    //partitioned sequence
    partition *Z;
    general_partition *X;

    void make_partition();
    
    void print();
    ~entry(){
        delete[] name;
        delete[] seq;
        delete Z;
        delete X;
        n=index=0;
    }
};


void set_program_options(options &opt, int argc, char**argv) ;
void fill_entries_randomly(const options opt, entry *entries);
void fill_entries_from_file(options &opts, entry *entries);
void create_temp_arrays(options opts);
dist_tuple distance(entry &e1, entry &e2) ;
dist_tuple general_distance(general_partition *p1, general_partition *p2);