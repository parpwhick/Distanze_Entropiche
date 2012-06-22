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
    
    adj_struct temp;
    temp.adj=adj;
    temp.adi=adi;
    temp.index=index;
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

    adj_struct temp;
    temp.adj=adj;
    temp.index=index;
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

    adj_struct temp;
    temp.adj=adj;
    temp.index=index;
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
    
    
    adj_struct temp;
    temp.index=index;
    temp.adj=adj;
    temp.adi=tmp_index;
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
