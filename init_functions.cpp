#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "strutture.h"
#include "rand_mersenne.h"
using namespace std;

options opts;


void error(const char* message) {
    fputs(message, stderr);
    exit(1);
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}


void print_help() {
    const char *message =
            "Usage: distanze [-option1 [arg]] [-option2 [arg]]...\n"
            "\n"
            "General option list for this program:\n"
            "  -random           Turns on random generation [on]\n"
            "  -file FILENAME    Read configurations from FILENAME [off]\n"
            "  -simulation       Configurations from temporal evolution\n"
            "  -num N            Limits the number of configurations to N [2550]\n"
            "  -nodistance       Doesn't calculate any distance matrix\n"
            "  -symbols N        Generate random strings with N symbols [2]\n"
            "  -nowrite          Don't write the distance matrices [off]\n"
            "  -threads N        Use N threads to calculate distance matrix [1]\n"
            //          "  -translate        Translate protein sequence into reduced 10 letter code [off]\n"
            "  -v [-v -v]        Turns on increasingly verbose messages to stderr [off]\n"
            "  -help             Shows this message\n"
            "  -demo             Only shows examples of partition algebra [off]\n"
            ;
    const char *message2 =
            "  -graphics         Make .ppm graphics when using square lattice topology\n"
            "  \n"
            "Reduction options:\n"
            "  -epsilon N        Use epsilon-reduction with unnormalized size N [0]\n"
            "  -common           Reduce using common partition\n"
            "  -direct           Reduce two partitions against each other [default]\n"
            "\n"
            "Options the choice of a topology:\n"
            "  -sequence N       Linear open sequence topology with 2 nearest neighbours long N\n"
            "  -fuzzy N          Linear open sequence with N nearest neighbours\n"
            "  -torus L          The configuration is from a periodic square of side L\n"
            "  -square L         The configuration is from an open square of side L [default, 25]\n"
            "  -sierpinski G     Use the Sierpinski gasket of generation G\n"
            "  -adj file1 file2  The files encode a generic adiacency matrix by its\n"
            "                    nonzero elements, with rows and cols in the files, eg \n"
            "                    -adj rows.bin cols.bin\n"
            "\n"
            "Simulation options:\n"
            "  -microcanonical   Evolution according to microcanonical law [default]\n"
            "  -metropolis       Evolution according to Metropolis rule\n"
            "  -creutz           Evolution according to Creutz rule\n"
            "  -beta B[,B2,..]   Floating point beta parameter for the (many) borders [0.45]\n"
            "  -sweeps N         Number of full sweeps in a time unit [1]\n"
            "  -skip N           Skip N time units to thermalize the system [50000]\n"
            "  -config_out FILE  Write generated configurations to FILE [configurations.bin]"
            "\n"
            "Quantum Hamiltonian:\n"
            "  -electrons        Calculate at each iteration the electron energy\n"
            "  -t X              Set hopping parameter to X [1.0]\n"
            "  -J X              Set J to X [0.01]\n"
            "  -V X              Set V to X [0]\n"
            "\n"
            ;
    fprintf(stderr, "%s", message);
    if (opts.partition_type == GENERAL_PARTITION)
        fprintf(stderr, "%s", message2);
}

void set_program_options(options &opts, int argc, char**argv) {
    opts.seq_len = 625;
    opts.lato = 25;
    opts.n_seq = 2550;
    opts.epsilon = 0;
    opts.n_symbols = 2;
    if (opts.partition_type == GENERAL_PARTITION)
        opts.topologia = RETICOLO_2D;
    else
        opts.topologia = LINEARE;
    opts.letto_da = RANDOM;
    opts.dynamics = MICROCANONICAL;
    opts.skip = 50000;
    opts.sweeps = 1;
    opts.graphics = false;
    opts.verbose = 0;
    opts.write = true;
    opts.distance = true;
    opts.fuzzy = 0;
    opts.riduzione = DIRETTA;
    opts.threads = 2;
    opts.demo = false;
    opts.electrons = false;
    opts.hopping = 1.0;
    opts.J = 0.01;
    opts.V = 0;
    opts.suffix_out = "";
    opts.da_calcolare = 0
            | SHAN | TOP
            | RID | RID_TOP
            // | HAMM
            ;

    int killswitch = 0;
    string command_line;
    for(int i=0; i<argc; i++){
        command_line += argv[i];
        command_line += " ";
    }
    opts.command_line=command_line;
    string input;
    if (argc > 1) {
        int read_argvs = 1;
        do {
            input = argv[read_argvs++];
            if (input == "-random") {
                fprintf(stderr, "Specifying random sequence generation\n");
                opts.letto_da = RANDOM;
            }else if (input == "-simulation") {
                fprintf(stderr, "Evolving configurations to analyze\n");
                opts.letto_da = SIMULATION;
            } else if (input == "-file") {
                if (argc - read_argvs < 1)
                    error("Missing filename to read!\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                strncpy(opts.state_filename, argv[read_argvs++], 255);
                fprintf(stderr, "Reading from filename: %s\n", opts.state_filename);
                opts.letto_da = FROM_FILE;
            } else if (input == "-adj") {
                if (argc - read_argvs < 2)
                    error("Need to specify two vector files\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                strncpy(opts.adj_vec_1, argv[read_argvs++], 255);
                strncpy(opts.adj_vec_2, argv[read_argvs++], 255);
                opts.topologia = FROM_FILE;
            } else if (input == "-sequence") {
                if (argc - read_argvs < 1)
                    error("Need to specify sequence length\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.seq_len = atoi(argv[read_argvs++]);
                opts.topologia = LINEARE;
                fprintf(stderr, "Analyzing 1d sequences long %d\n", opts.seq_len);
            } else if (input == "-epsilon") {
                if (argc - read_argvs < 1)
                    error("Need to specify epsilon\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.epsilon = atoi(argv[read_argvs++]);
                fprintf(stderr, "Reduction uses epsilon algorithm, with size %d\n", opts.epsilon);
            } else if (input == "-direct") {
                fprintf(stderr, "Using reduction against each other\n");
                opts.riduzione = DIRETTA;
            } else if (input == "-common") {
                fprintf(stderr, "Using common-factor reduction\n");
                opts.riduzione = COMUNE;
            } else if (input == "-square") {
                if (argc - read_argvs < 1)
                    error("Need to specify lattice side length\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.topologia = RETICOLO_2D;
                opts.lato = atoi(argv[read_argvs++]);
                fprintf(stderr, "Using square lattice topology, with side %d\n", opts.lato);
            } else if (input == "-torus") {
                if (argc - read_argvs < 1)
                    error("Need to specify lattice side length\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.topologia = TORO_2D;
                opts.lato = atoi(argv[read_argvs++]);
                fprintf(stderr, "Using toroidal lattice topology, with side %d\n", opts.lato);
            } else if (input == "-sierpinski") {
                if (argc - read_argvs < 1)
                    error("Need to specify generation\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.topologia = SIERPINSKI;
                opts.sierpinski_gen = atoi(argv[read_argvs++]);
                fprintf(stderr, "Using Sierpinski generation %d\n",opts.sierpinski_gen);
            }
            else if (input == "-num") {
                if (argc - read_argvs < 1)
                    error("Missing max number of configurations\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.n_seq = atoi(argv[read_argvs++]);
                fprintf(stderr, "Number of configurations limited to: %d\n", opts.n_seq);
            } else if (input == "-fuzzy") {
                if (argc - read_argvs < 1)
                    error("Missing number of partitioning fuzziness\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.epsilon = atoi(argv[read_argvs++]);
                if (opts.epsilon > 0)
                    opts.topologia = FUZZY;
                else
                    opts.topologia = LINEARE;
                fprintf(stderr, "Fuzziness degree set to: %d\n", opts.epsilon);
            }else if (input == "-beta") {
                if (argc - read_argvs < 1)
                    error("Missing beta parameter\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                std::vector<string> input_beta = split(argv[read_argvs++],',');
                for(size_t i=0; i<input_beta.size(); i++)
                    opts.beta.push_back(atof(input_beta[i].c_str()));
                fprintf(stderr, "Beta set to: ");
                for(size_t i=0; i<opts.beta.size(); i++)
                    fprintf(stderr,"%.2f ",opts.beta[i]);
                fprintf(stderr,"\n");
            } else if (input == "-sweeps") {
                if (argc - read_argvs < 1)
                    error("Missing number of sweeps\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.sweeps = atoi(argv[read_argvs++]);
                fprintf(stderr, "New configuration every %d sweeps\n", opts.sweeps);
            } else if (input == "-skip") {
                if (argc - read_argvs < 1)
                    error("Missing number of skipped configurations\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.skip = atoi(argv[read_argvs++]);
                fprintf(stderr, "Skipping %d configurations at beginning\n", opts.skip);
            } else if (input == "-suffix") {
                if (argc - read_argvs < 1)
                    error("Missing file name\n");

                opts.suffix_out = argv[read_argvs++];
                fprintf(stderr, "Writing run data with '%s' appended\n", opts.suffix_out.c_str());
            } else if (input == "-microcanonical") {
                fprintf(stderr, "Temporal series with microcanonical rule\n");
                opts.dynamics = MICROCANONICAL;
                opts.letto_da = SIMULATION;
            } else if (input == "-metropolis") {
                fprintf(stderr, "Temporal with Metropolis rule\n");
                opts.dynamics = METROPOLIS;
                opts.letto_da = SIMULATION;
            } else if (input == "-creutz") {
                fprintf(stderr, "Temporal series with Creutz rule\n");
                opts.dynamics = CREUTZ;
                opts.letto_da = SIMULATION;
            } else if (input == "-electrons") {
                fprintf(stderr, "Calculating Hamiltonian at each iteration\n");
                opts.electrons = true;
            } else if (input == "-t") {
                if (argc - read_argvs < 1)
                    error("Missing t\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.hopping = atof(argv[read_argvs++]);
                fprintf(stderr, "Hoping parameter set to %f\n", opts.hopping);
            } else if (input == "-J") {
                if (argc - read_argvs < 1)
                    error("Missing J\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.J = atof(argv[read_argvs++]);
                fprintf(stderr, "J interaction set to %f\n", opts.J);
            } else if (input == "-V") {
                if (argc - read_argvs < 1)
                    error("Missing J\n");

                opts.V = atof(argv[read_argvs++]);
                fprintf(stderr, "Voltage set to %f\n", opts.V);
            } else if (input == "-v") {
                opts.verbose++;
                fprintf(stderr, "Verbosity at %d\n", opts.verbose);
            } else if (input == "-graphics") {
                opts.graphics = true;
                fprintf(stderr, "Making pretty lattice graphs\n");
            } else if (input == "-hamming") {
                opts.da_calcolare |= HAMM;
                fprintf(stderr, "Hamming distance\n");
            } else if (input == "-nowrite") {
                opts.write = false;
                fprintf(stderr, "Not writing out the distance matrices\n");
            } else if (input == "-nodistance") {
                opts.distance = false;
                fprintf(stderr, "Not calculating the distance matrix\n");
            } else if (input == "-symbols") {
                if (argc - read_argvs < 1)
                    error("Missing max number of random letters to use\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.n_symbols = atoi(argv[read_argvs++]);
                fprintf(stderr, "Number of letters limited to: %d\n", opts.n_symbols);
            } else if (input == "-threads") {
                if (argc - read_argvs < 1)
                    error("Expecting thread number\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.threads = atoi(argv[read_argvs++]);
                fprintf(stderr, "Threads limited to: %d\n", opts.threads);
            } else if (input == "-h" || input == "-help" || input == "--help") {
                print_help();
                killswitch = 1;
            } else if (input == "-demo") {
                fprintf(stderr, "Demo mode on!\n");
                opts.demo = true;
                opts.n_seq = 6;
            } else {
                fprintf(stderr, "Unknown option: %s\n", input.c_str());
                print_help();
                killswitch = 1;
            }
        } while (argc - read_argvs > 0);
    } else {
        print_help();
        killswitch = 1;
    }
    if (killswitch)
        exit(0);
    if(opts.beta.empty())
        opts.beta.push_back(0.44);

    fprintf(stderr, "\n");

}

static RandMT gen;

void fill_entries_randomly(std::string *entries) {
    for (int i = 0; i < opts.n_seq; i++)
        for (int j = 0; j < opts.seq_len; j++)
            entries[i] += (gen.get_int() % opts.n_symbols); //+ 'A');

}

void generate_next_sequence(int *entry) {
    for (int j = 0; j < opts.seq_len; j++)
        entry[j] = (gen.get_int() % opts.n_symbols);
}

int load_config(options &opts, int *num_entry) {
    static FILE* in = 0;
    static int full_configs = 0;
    static int read_count = 0;
    static int finished = 0;

    //if they've already been all loaded, do nothing
    if (finished)
        return 0;

    //on the first call, open the file, calculate how many configurations to read
    if (in == 0) {
        in = fopen(opts.state_filename, "r");
        if (in == 0) {
            fprintf(stderr, "LATTICE LOAD: Error opening file %s\n", opts.state_filename);
            exit(2);
        }
        long M;
        (void) fseek(in, 0L, SEEK_END);
        M = ftell(in);
        rewind(in);

        M /= sizeof (int32_t);
        full_configs = M / opts.seq_len;
        if (M % opts.seq_len)
            fprintf(stderr, "LATTICE LOAD: Warning, found only %d full configurations, %d integers left\n", full_configs, (int) M % opts.seq_len);

        if (full_configs <= opts.n_seq)
            opts.n_seq = full_configs;
        else {
            fprintf(stderr, "LATTICE LOAD: Full %d configurations, %d requested\n", full_configs, opts.n_seq);
            full_configs = opts.n_seq;
        }
    }

    //try reading, and if all goes well, return
    int read = fread(num_entry, sizeof (int32_t), opts.seq_len, in);
    if (feof(in)) { //read < opts.seq_len){
        fprintf(stderr, "LATTICE LOAD: Warning, unexpected end of file\n");
        opts.n_seq = read_count;
        finished = 1;
        return 0;
    } else if (read < opts.seq_len) {
        fprintf(stderr, "LATTICE LOAD: Read too little on configuration %d\n", read_count + 1);
        opts.n_seq = read_count;
        finished = 1;
        return 0;
    } else
        read_count++;

    if (read_count == full_configs) {
        fclose(in);
        finished = 1;
    }

    return (1);
}

void fill_seq_from_file(options &opts, std::string *sequenze) {
    int max_length = 0;
    string buffer;
    buffer.reserve(150);

    ifstream in(opts.state_filename);
    if (!in) {
        fprintf(stderr, "Can't open file %s\n", opts.state_filename);
        exit(2);
    }

    //vector<string> V((istream_iterator<string>(cin)), istream_iterator<string>());
    int cur_entry = -1;
    while (1) {
        //get every line and check if read properly
        getline(in, buffer);
        if (in.eof()) {
            break;
        }
        int len = buffer.length();
        //skip comments in the sequence file, i.e. lines marked with #
        if (buffer[0] == '#' || len == 0) {
            cout << "Skipping comment or whiteline " << endl;
            continue;
        }
        if (buffer[len - 1] == '\r')
            buffer.resize(len - 1);
        //a line beginning with ">" is the name of a new sequence
        if (buffer[0] == '>') {
            if (cur_entry == opts.n_seq - 1)
                break;
            cur_entry += 1;
            sequenze[cur_entry].reserve(opts.seq_len);
        }
        else
            sequenze[cur_entry].append(buffer);
    }

    //adding 1 to set the correct number of entries (0..cur_entry -> makes for cur_entry+1 entries)
    cur_entry += 1;
    opts.n_seq = cur_entry;

    for (int i = 0; i < opts.n_seq; i++)
        max_length = std::max(max_length, (int) sequenze[i].size());

    if (max_length > opts.seq_len) {
        max_length = opts.seq_len;
        for (int i = 0; i < opts.n_seq; i++)
            sequenze[i].resize(opts.seq_len);
    }

    fprintf(stderr, "Read %d sequences, the longest was %d bytes\n\n", cur_entry, max_length);
    opts.seq_len = max_length;

}

int *colore;
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
            colore[i] = gen.get_int();
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
            colore[i] = gen.get_int();
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

template void ppmout(const label_t *grid, int sz, const char *filename);
template void ppmout(const product_t *grid, int sz, const char *filename);
template void ppmout2(const label_t *grid1, const label_t* grid2, int, const char *);
template void ppmout2(const product_t *grid1, const product_t* grid2, int, const char *);
template void ppmout2(const char *grid1, const label_t* grid2, int, const char *);
