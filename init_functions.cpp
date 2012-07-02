#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>

#include "strutture.h"
#include "rand_mersenne.h"
using namespace std;

options opts;


void translate(string type, char *testo, int len);

void error(const char* message) {
    fputs(message, stderr);
    exit(1);
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
            "  -length N         Limits configuration length to N [600]\n"
            "  -nodistance       Doesn't calculate any distance matrix\n"
            "  -symbols N        Generate random strings with N symbols [2]\n"
            "  -nowrite          Don't write the distance matrices [off]\n"
            "  -threads N        Use N threads to calculate distance matrix [1]\n"
            //          "  -translate        Translate protein sequence into reduced 10 letter code [off]\n"
            "  -v [-v -v]        Turns on increasingly verbose messages to stderr [off]\n"
            "  -help             Shows this message\n"
            ;
    const char *message2 =
            "  -graphics         Make .ppm graphics when using square lattice topology\n"
            "  \n"
            "Options the choice of a topology:\n"
            "  -sequence         Linear open sequence topology with 2 nearest neighbours\n"
            "  -fuzzy N          Linear open sequence with N nearest neighbours\n"
            "  -square L         The configuration is from a square of side L [default, 25]\n"
            "  -adj file1 file2  The files encode a generic adiacency matrix by its\n"
            "                    nonzero elements, with rows and cols in the files, eg \n"
            "                    -adj rows.bin cols.bin\n"
            "\n"
            "Simulation options:\n"
            "  -microcanonical   Evolution according to microcanonical law [default]\n"
            "  -link_energy N    Max of the microcanonical kinetic energy [10]\n"
            "  -metropolis       Evolution according to Metropolis rule\n"
            "  -beta B           Floating point beta parameter [0.45]\n"
            ;
    fprintf(stderr, "%s", message);
    if (opts.partition_type == GENERAL_PARTITION)
        fprintf(stderr, "%s", message2);
}

void set_program_options(options &opts, int argc, char**argv) {
    opts.seq_len = 625;
    opts.lato = 25;
    opts.n_seq = 2550;
    opts.n_symbols = 2;
    opts.topologia = RETICOLO_2D;
    opts.letto_da = RANDOM;
    opts.simulation_type=MICROCANONICAL;
    opts.beta = 0.45;
    opts.max_link_energy = 10;
    opts.graphics = false;
    opts.verbose = 0;
    opts.write = true;
    opts.distance = true;
    opts.fuzzy = 0;
    opts.threads = 2;
    opts.da_calcolare = 0
            | SHAN | TOP
            | RID | RID_TOP
            // | HAMM
            ;

    int killswitch = 0;

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
            } else if (input == "-sequence") {
                fprintf(stderr, "Analyzing 1d sequences\n");
                opts.topologia = LINEARE;
            } else if (input == "-adj") {
                if (argc - read_argvs < 2)
                    error("Need to specify two vector files\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                strncpy(opts.adj_vec_1, argv[read_argvs++], 255);
                strncpy(opts.adj_vec_2, argv[read_argvs++], 255);
                opts.topologia = FROM_FILE;
            } else if (input == "-length") {
                if (argc - read_argvs < 1)
                    error("Need to specify sequence length\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.seq_len = atoi(argv[read_argvs++]);
                fprintf(stderr, "Sequence length limited to %d\n", opts.seq_len);
            } else if (input == "-square") {
                if (argc - read_argvs < 1)
                    error("Need to specify lattice side length\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.topologia = RETICOLO_2D;
                opts.lato = atoi(argv[read_argvs++]);
                fprintf(stderr, "Using square lattice topology, with side %d\n", opts.lato);
            }
            else if (input == "-num") {
                if (argc - read_argvs < 1)
                    error("Missing max number of sequences to read\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.n_seq = atoi(argv[read_argvs++]);
                fprintf(stderr, "Number of sequences limited to: %d\n", opts.n_seq);
            } else if (input == "-fuzzy") {
                if (argc - read_argvs < 1)
                    error("Missing number of partitioning fuzziness\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.fuzzy = atoi(argv[read_argvs++]);
                if (opts.fuzzy > 0)
                    opts.topologia = FUZZY;
                else
                    opts.topologia = LINEARE;
                fprintf(stderr, "Fuzziness degree set to: %d\n", opts.fuzzy);
            }else if (input == "-beta") {
                if (argc - read_argvs < 1)
                    error("Missing beta parameter\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.beta = atof(argv[read_argvs++]);
                fprintf(stderr, "Beta set to: %.2f\n", opts.beta);
            } else if (input == "-link_energy") {
                if (argc - read_argvs < 1)
                    error("Missing max link energy\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.max_link_energy = atoi(argv[read_argvs++]);
                fprintf(stderr, "Max link energy set to: %d\n", opts.max_link_energy);
            } else if (input == "-microcanonical") {
                fprintf(stderr, "Using microcanonical rule\n");
                opts.simulation_type = MICROCANONICAL;
                opts.letto_da = SIMULATION;
            } else if (input == "-metropolis") {
                fprintf(stderr, "Using Metropolis rule\n");
                opts.simulation_type = METROPOLIS;
                opts.letto_da = SIMULATION;
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

    if (opts.topologia == RETICOLO_2D)
        opts.seq_len = opts.lato * opts.lato;

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
