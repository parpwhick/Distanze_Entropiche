//  CREATION OF THE RANDOM LATTICE
    if (opts.from == LATTICE) {
        
        //fill_entries_randomly(opts, entries);
        colore = new int[opts.seq_len];
        for (int i = 0; i < opts.seq_len; i++) {
            colore[i] = xrand();
        }
        colore[0] = 0x0F5A3A1F; // blue almost black
    }
        
    //  CREATION OF RANDOM SEQUENCES (as entries)
    //
    if(opts.from==RANDOM){
        fill_entries_randomly(opts,entries);
    }
    //
    //  READ SEQUENCES FROM FILE
    //
    if(opts.from==SEQUENCE)
        //fill_entries_from_file_streams(opts,entries);
        fill_seq_from_file(opts,entries);
    
    if(opts.verbose)
        fprintf(stderr,"Entries filled\n");
    //
    //  ANALISYS OF LOADED SEQUENCES
    switch (opts.from) {
        case(LATTICE):
            for (i = 0; i < opts.n_seq; i++){
                Z[i].from_square_lattice(entries[i].c_str(), opts.lato, 2);                
            }
            da_calcolare &= ~(SHAN | SHAN_TOP | RID | RID_TOP);
            opts.da_calcolare=da_calcolare;
            break;
        default:
        case(SEQUENCE):
            for (i = 0; i < opts.n_seq; i++) {
                X[i].fill(entries[i].c_str(), entries[i].size());
                Z[i].from_linear_sequence(entries[i].c_str(), entries[i].size());
            }
            print_partition_stats(X,"semplici");
            break;
    }
