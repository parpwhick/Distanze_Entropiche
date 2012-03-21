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

    int seed;
    int verbose;
    int n_symbols;
    bool translate;

}options;


struct dist_tuple {
    double dist;
    double dist_r;
    double dist_top;
    double dist_top_r;
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
    partition() {
        atom_positions = 0;
        binary = 0;
        n = 0;
        N = 0;
    }
};

typedef char* stringa;
class entry {
public:
    //identifier of the sequence
    stringa name;
    //string of symbols making the sequence
    stringa seq;
    //length of the sequence
    int n;
    //sequence index, relevant only for printing
    int index;
    //partitioned sequence
    partition Z;

    void make_partition();
    void print();
};




void set_program_options(options &opt, int argc, char**argv) ;
void fill_entries_randomly(const options opt, entry *entries);
void fill_entries_from_file(options &opt, entry *entries);
