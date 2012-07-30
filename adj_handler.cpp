/* 
 * File:   adj_handler.cpp
 * Author: fake
 * 
 * Created on May 27, 2012, 8:31 PM
 */

#include "adj_handler.h"
#include "strutture.h"

extern options opts;

#define nnu(i) (i - (i % lato)+ ((i+lato-1)%lato))
#define nnd(i) ((i/lato)*lato + ((i+lato+1)%lato))
#define nnl(i) (i+N-lato)%N
#define nnr(i) (i+N+lato)%N

adj_struct adiacenza_square_lattice(int lato){    
    int N=lato*lato;

    int *adj=new int[4*N+1];
    int *adi=new int[4*N+1];
    int *index=new int[N+1];
    for (int i=0; i<N; i++){
        adj[4*i]   =  nnu(i);
        adj[4*i+1] =  nnl(i);
        adj[4*i+2] =  nnd(i);
        adj[4*i+3] =  nnr(i);
        
        adi[4*i]   =  i;
        adi[4*i+1] =  i;
        adi[4*i+2] =  i;
        adi[4*i+3] =  i;
        
        index[i]=4*i;
    }
    adj[4*N]=LEAST;
    index[N]=4*N;
    
    adj_struct temp(adi,adj,index);
    temp.N=N;
    temp.n_link=4*N;
    temp.zmax=4;
    return(temp);
}

adj_struct adiacenza_simple_line(int N){ 
    int *adj=new int[N+1];
    int *index=new int[N+1];
    
    adj[0]=LEAST;
    index[0]=0;
    for (int i=1; i<N; i++){
        adj[i]=i-1;
        index[i]=i;
    }
    adj[N]=LEAST;
    index[N]=N;

    adj_struct temp(nullptr,adj,index);
    temp.N=N;
    temp.n_link=N;
    temp.zmax=1;
    return(temp);
}

adj_struct adiacenza_fuzzy_line(int N){
    int *adj=new int[(opts.fuzzy+1)*N];
    int *index=new int[N+1];
    int adj_count=0;
    
    adj[0]=LEAST;
    index[0]=0;
    for (int i=1; i<N; i++){
        index[i]=adj_count;
        adj[adj_count++]=-(i-1);
    
        for(int j=2; i-j>=0 && j<=opts.fuzzy+1; j++)
                adj[adj_count++]=i-j;
    }
    adj[adj_count]=-1;
    index[N]=adj_count;

    adj_struct temp(nullptr,adj,index);
    temp.N=N;
    temp.n_link=adj_count;
    temp.zmax=opts.fuzzy+1;
    return(temp);
}

/*
char *colori=new char[total_size];
    colori[0]=1;
    
    for (int i=1; i<total_size; i++){
        int min_col=1;
    
        for(int j=0; j<4;j++)
                if(nn(i,j) != -1 && colori[nn(i,j)] == min_col) min_col++ ;
        colori[i]=min_col;        
    }
    int popcol[4];
    for(int i=0; i<4;i++) popcol[i]=0;
    
    for (int i=0; i<total_size; i++)
        popcol[colori[i]-1]++;
    
    for(int i=0; i<4;i++) 
        printf("Popolazione di %d: %.2f%%\n",i+1,(popcol[i]*100.0)/total_size);
*/  

#include "matrix.h"
adj_struct adiacenza_sierpinski(int GEN, int &total_size){
    int a0, a1, a2;
    int size, old_size;

    total_size = 6;
    for (int g = 2; g <= GEN; g++)
        total_size = 3 * total_size - 3;
    fprintf(stderr, "Generazione %d, size %d, links %d\n", GEN, total_size, total_size * 4 - 6);

    // prima generazione, dimensione 6
    old_size = 6;
    size = 6;

    /*allocazione */
    
    // vettore di coordinazione
    memory < 1 > z(total_size);
    // matrice (n,4) nearest neighbour
    memory < 4 > nn(total_size, -1);
    
    /*riempimento del triangolo generatore */
    z(0) = 2;
    nn(0, 0) = 1;
    nn(0, 1) = 2;

    z(1) = 4;
    nn(1, 0) = 2;
    nn(1, 1) = 0;
    nn(1, 2) = 3;
    nn(1, 3) = 4;


    z(2) = 4;
    nn(2, 0) = 0;
    nn(2, 1) = 1;
    nn(2, 2) = 4;
    nn(2, 3) = 5;

    z(3) = 2;
    nn(3, 0) = 1;
    nn(3, 1) = 4;

    z(4) = 4;
    nn(4, 0) = 3;
    nn(4, 1) = 1;
    nn(4, 2) = 2;
    nn(4, 3) = 5;

    z(5) = 2;
    nn(5, 0) = 2;
    nn(5, 1) = 4;

    /*angoli*/
    a0 = 0;
    a1 = 3;
    a2 = 5;

    /*cicli per la costruzione ed il riempimento delle generazioni successive*/

    for (int g = 2; g <= GEN; g++) {
        size = 3 * old_size - 3;
        /*Costruzione dell gasket della generazione g*/
        
        //triangolo 1 - copia identica
        //nulla da fare
        
        //triangolo 2, dimensione: oldsize-2
        //for (int n=oldsize; n<2*oldsize-2; n++)
        for (int n = 1; n < old_size - 1; n++) {
            z(n + old_size - 1) = z(n);
            for (int m = 0; m < z(n + old_size - 1); m++) {
                if (nn(n, m) == a0)
                    nn(n + old_size - 1, m) = a1;
                else if (nn(n, m) == a2)
                    nn(n + old_size - 1, m) = a1 + 2 * old_size - 3;
                else nn(n + old_size - 1, m) = nn(n, m) + old_size - 1;
            }
        }
        //triangolo 3, dimensione: oldsize-1
        //for (int n=2*oldsize-2; n<3*oldsize-3; i++)
        for (int n = 1; n < old_size; n++) {
            z(n + 2 * old_size - 3) = z(n);
            for (int m = 0; m < z(n + 2 * old_size - 3); m++) {
                if (nn(n, m) == a0)
                    nn(n + 2 * old_size - 3, m) = a2;
                else nn(n + 2 * old_size - 3, m) = nn(n, m) + 2 * old_size - 3;
            }
        }

        /* Vengono creati i link mancanti del gasket della generazione g*/		
        z(a1) = 4;
        nn(a1, 2) = old_size;
        nn(a1, 3) = old_size + 1;
        z(a1 + 2 * old_size - 3) = 4;
        nn(a1 + 2 * old_size - 3, 2) = nn(a2, 0) + old_size - 1;
        nn(a1 + 2 * old_size - 3, 3) = nn(a2, 1) + old_size - 1;
        z(a2) = 4;
        nn(a2, 2) = 2 * old_size - 2;
        nn(a2, 3) = 2 * old_size - 1;
		

        /*Individuo i siti di bordo del gasket generazione g*/
        a1 = a1 + old_size - 1;
        a2 = a2 + 2 * old_size - 3;
        old_size = size;


    }
    int *adi = new int[4 * total_size];
    int *adj = new int[4 * total_size];
    int *index = new int[total_size+1];
    int scritti = 0;
    for (int i = 0; i < total_size; i++) {
        index[i]=scritti;
        for (int j = 0; j < z(i); j++) {
            adi[scritti + j] = i;
            adj[scritti + j] = nn(i, j);
        }
        scritti += z(i);
    }
    index[total_size] = scritti;
    adj_struct temp(adi,adj,index);
    temp.N=total_size;
    temp.n_link=scritti;
    temp.zmax=4;
    
//    // test!
//    for(int i=0; i < temp.N; i++){
//        int z=temp.fetch(i);
//        printf("z(%d)=%d\n",i+1,z);
//    }

    return(temp);
}

#include <iostream>
adj_struct adiacenza_from_file(const char *name_vec1,const char *name_vec2, int & N){
    FILE *vec1=fopen(name_vec1,"rb");
    FILE *vec2=fopen(name_vec2,"rb");
    if(vec1==0 || vec2==0){
        fprintf(stderr,"ADJ READ: Error reading adjacency vectors\n");
        exit(1);
    }
    
    long M,M1;
    (void) fseek(vec1, 0L, SEEK_END);
    M=ftell(vec1);
    rewind(vec1);
    
    (void) fseek(vec2, 0L, SEEK_END);
    M1=ftell(vec2);
    rewind(vec2);
    
    
    if(M != M1){
        fprintf(stderr,"ADJ READ: Number of indexes differs from the number of values!");
        exit(1);
    }
    
    // voglio il numero di interi da leggere, non di byte
    M /= sizeof(int32_t);
    
    int *adj=new int[M+1];
    int *tmp_index = new int [M+1];
    int zmax=0;
    
    
    //i valori sono in vec2
    int T1= fread(adj,sizeof(int32_t),(int)M,vec2);
    //gli indici sono in vec1
    int T2= fread(tmp_index,sizeof(int32_t),(int)M,vec1);
    if(T1<M || T2<M){
            fprintf(stderr,"ADJ READ: Error: read %d links, %d indexes, wanted %d\n",T1,T2, (int)M);
            exit(1);
    }
    
    //try detecting matlab's offset
    int offset=tmp_index[0];
    if(offset != 1){
        fprintf(stderr,"ADJ READ Warning: The index is %d-based, not 1-based\n",offset);
    }
    
    //try detecting number of sites
    N=tmp_index[M-1]-offset+1;
    fprintf(stderr,"ADJ READ Info: Reading vectors for %d elements, %ld nonempty links\n",N,M);
    int *index=new int[N+1];    
       
    if(opts.verbose)
        fprintf(stderr,"ADJ READ Info: Finished allocating\n");
    index[0]=0;
    int count=0;
    adj[0] -= offset;
    for(int i=1; i<M; i++){
        //correcting for possible Matlab's offset
        adj[i] -= offset; 
        if(tmp_index[i]!=tmp_index[i-1]){
            count++;
            if(tmp_index[i]-offset >=N || adj[i] >= N){
                fprintf(stderr,"ADJ READ Error: Links to more than the expected number of sites\n");
                exit(1);
            }
            if(tmp_index[i]<tmp_index[i-1]) {
                fprintf(stderr,"ADJ READ Error: Index file not sorted. Hint: try switching rows and columns\n");
                exit(1);
            }
            if(tmp_index[i]-tmp_index[i-1]>1){
                //fprintf(stderr,"ADJ READ Warning: Sites from %d to %d isolated\n",tmp_index[i-1]+1,tmp_index[i]-1);
                fprintf(stderr,"ADJ READ Warning: Site %d isolated\n",tmp_index[i-1]+1);
                for(; count < tmp_index[i]-offset; count++)
                    index[count]=i;
                
            }                
            index[count]=i;
        }   
    }
    index[N]=M;    
    
    /* find zmax*/
    for (int i=0; i<N; i++){
        int tmp=index[i+1]-index[i];
        zmax = (tmp>zmax) ? tmp: zmax;
    }
    
    /* correct offset in second vector */
    for(int i=0; i<M; i++)
        tmp_index[i] -= offset;    
    
    
    adj_struct temp(tmp_index,adj,index);
    temp.n_link=M;
    temp.N=N;
    temp.zmax=zmax;
    
    fclose(vec1);
    fclose(vec2);
    return(temp);
}

/* Pseudo iterator
void neigh_factory::f3() {
    n = 0;
    int s1;
    int quanti = 1;
    if (adj[adj_counter] == LEAST)
        quanti = 0;

    while (adj[adj_counter + quanti + 1] > 0)
        quanti++;
    for (int i = 0; i < quanti; i++) {
        s1 = adj[adj_counter++];
        if (i == 0)
            s1 = -s1;
        //if(s1 < 0 || s1==site || 
        if (configuration[_site] != configuration[s1])
            continue;
        buffer[n++] = s1;
    }
    _site++;
    if (_site >= N)
        _site = -1;
}
 */
