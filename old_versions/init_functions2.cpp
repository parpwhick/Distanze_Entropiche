#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>


#include "strutture2.h"


void error(const char* message){
    fputs(message,stderr);
    exit(1);
}

void print_help() {
    fprintf(stderr,
            "Usage: distanze [-option1 [arg]] [-option2 [arg]]...\n"
            "\n"
            "Option list for this program:\n"
            "  -random          Turns on random string generation [on]\n"
            "  -file FILENAME   Read sequences from FILENAME [off]\n"
            "  -seqnum N        Limits the number of sequences to N [1000]\n"
            "  -seqlength N     Limits sequence length to N [1000]\n"
            "  -translate       Translate protein sequence into reduced 5 letter code [off]\n"
            "  -symbols N       Generate random strings with N symbols [2]\n"
            "  -write           Write the distance matrices [off]\n"
            "  -seed N          Random number generator seed [37337]\n"
            "  -v [-v -v]       Turns on increasingly verbose messages to stderr [off]\n"
            "  -help            Shows this message\n"
            );
}
void set_program_options(options &opts, int argc, char**argv) {
    opts.seq_len = 1000;
    opts.n_seq = 1000;
    opts.n_symbols = 2;
    opts.random_sequences = true;
    opts.seed = 37337;
    opts.translate = false;
    opts.verbose = 0;
    opts.write=false;
    opts.distance=true;
    int killswitch=0;

    string input;
    if (argc > 1) {
        int read_argvs = 1;
        do {
            input = argv[read_argvs++];
            if (input == "-random") {
                fprintf(stderr, "Specifying random sequence generation\n");
                opts.random_sequences = true;
            } else if (input == "-file") {
                if (argc - read_argvs < 1)
                    error("Missing filename to read!\n");
                if (argv[read_argvs][0] == '-') 
                    error("Expecting argument, not another option\n");
                
                strncpy(opts.filename, argv[read_argvs++], 255);
                fprintf(stderr, "Reading from filename: %s\n", opts.filename);
                opts.random_sequences = false;
            } else if (input == "-seqlength") {
                if (argc - read_argvs < 1)
                    error("Need to specify sequence length\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");
                
                opts.seq_len = atoi(argv[read_argvs++]);
                fprintf(stderr, "Sequence length limited to %d\n", opts.seq_len);
            } else if (input == "-seqnum") {
                if (argc - read_argvs < 1)
                    error("Missing max number of sequences to read\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.n_seq = atoi(argv[read_argvs++]);
                fprintf(stderr, "Number of sequences limited to: %d\n", opts.n_seq);
            } else if (input == "-v") {
                opts.verbose++;
                fprintf(stderr, "Verbosity at %d\n", opts.verbose);
            } else if (input == "-translate") {
                opts.translate = true;
                fprintf(stderr, "Simplifying sequence alphabet\n");
            } else if (input == "-write") {
                opts.write = true;
                fprintf(stderr, "Writing out the distance matrices\n");
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
            } else if (input == "-seed") {
                if (argc - read_argvs < 1)
                    error("Expecting random number generation seed\n");
                if (argv[read_argvs][0] == '-')
                    error("Expecting argument, not another option\n");

                opts.seed = atoi(argv[read_argvs++]);
                fprintf(stderr, "Seed set to: %d\n", opts.seed);
            } else if (input == "-h" || input == "-help" || input == "--help") {
                print_help();
                killswitch=1;
            }
            else fprintf(stderr, "Unknown option: %s\n", input.c_str());
        } while (argc - read_argvs > 0);
    }
    else{
        print_help();
        killswitch=1;
    }
    if(killswitch)
        exit(0);
    
    srand(opts.seed);

}

void fill_entries_randomly(const options opts, entry *entries){
for (int i = 0; i < opts.n_seq; i++) {
        char *name = new char[10];
        char *seq = new char[opts.seq_len + 1];
        snprintf(name, 10, "seq%03d", i + 1);
        entries[i].name = name;
        for (int j = 0; j < opts.seq_len; j++)
            seq[j] = rand() % opts.n_symbols + 'A';
        seq[opts.seq_len + 1] = 0;
        entries[i].seq = seq;
        entries[i].n = opts.seq_len;
    }
}

void fill_entries_from_file(options &opts, entry *entries){
    char *line_buffer = new char[1001];
    char *return_code;
    int max_length = 0, current_length = 0, line_length = 0;
    FILE *in;


    in = fopen(opts.filename, "r");
    if (in == 0) {
        fprintf(stderr, "Can't open file %s\n", opts.filename);
        exit(2);
    }

    int cur_entry = -1;
    while (1) {
        //get every line and check if read properly
        return_code = fgets(line_buffer, 1000, in);
        if (!return_code) {
            if(cur_entry<1){
            fprintf(stderr, "Problem reading line\n");
            exit(3);
            }
            else{
                break;
            }
        }

        //skip whitespace at the end, doesn't check for whitespace inside
        line_length = strnlen(line_buffer, 1000);
        while ((line_buffer[line_length - 1] == '\n') ||
                (line_buffer[line_length - 1] == '\r') ||
                (line_buffer[line_length - 1] == ' ')) {
            line_length -= 1;
            line_buffer[line_length] = 0;
        }

        //skip comments in the sequence file, i.e. lines marked with #
        if (line_buffer[0] == '#')
            continue;

        //a line beginning with ">" is the name of a new sequence
        if (line_buffer[0] == '>') {
            cur_entry += 1;
            if (!(cur_entry < opts.n_seq))
                break;
            entries[cur_entry].name = strdup(line_buffer + 1);
            entries[cur_entry].index=cur_entry+1;
            entries[cur_entry].seq = new char[opts.seq_len + 1];
            for(int i=0;i<opts.seq_len;i++)
                entries[cur_entry].seq[i]=' ';
            entries[cur_entry].seq[opts.seq_len]=0;
            current_length = 0;            
        } else {
            int j = 0;
            //safely copy sequence bytes, never going over the hard-coded limit
            //to avoid buffer-overflows
            for (j = 0; j < line_length && j + current_length < opts.seq_len; j++)
                entries[cur_entry].seq[current_length + j] = line_buffer[j];
            current_length += j;
            entries[cur_entry].seq[current_length] = 0;
            entries[cur_entry].n=current_length;
        }
        max_length = (current_length > max_length) ? current_length : max_length;
    }
    //adding 1 to set the correct number of entries (0..cur_entry -> makes for cur_entry+1 entries)
    cur_entry+=1; 
    if (opts.verbose) 
        fprintf(stderr, "Read %d sequences, the longest was %d bytes\n\n", cur_entry, max_length);

    opts.n_seq = cur_entry;
    opts.seq_len=max_length;
    
    if (opts.verbose > 2)
        for (int i = 0; i < opts.n_seq; i++) {
                printf("Seq[%d].name=%s\n", i, entries[i].name);
                printf("Seq[%d].seq=%s\n", i, entries[i].seq);
        }

    fclose(in);
    delete line_buffer;
}



