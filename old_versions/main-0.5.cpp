/* 
 * File:   main.cpp
 * Author: fake
 *
 * Created on January 11, 2012, 2:30 PM
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#include "easy_vect.hpp"

using namespace std;

int minimum(int a, int b, int c)
/*Gets the minimum of three values*/ {
    int min = a;
    if (b < min)
        min = b;
    if (c < min)
        min = c;
    return min;
}

class atom {
public:
    int start;
    int n;
    char type;

    atom(int in, int num, char a) {
        start = in;
        n = num;
        type = a;
    }

    atom(int in, int num) {
        start = in;
        n = num;
        type = 0;
    }

    atom(int in) {
        start = 0;
        n = 0;
        type = 0;
    }

    atom() {
        start = 0;
        n = 0;
        type = 0;
    }

    void print() {
        if (type) {
            fprintf(stderr, "(");
            for (int i = 0; i < n; i++)
                fprintf(stderr, "%c", type);
            fprintf(stderr, "),");
        } else {
            fprintf(stderr, "%d (%d),", start, n);
        }
    }
};

class partizione {
public:
    easy_vect <atom> atoms;
    int *simple;
    int *simple_binary;
    bool is_simple;
    int n;
    int N;

    double entropia_shannon();

    int entropia_topologica() {
        return n;
    }

    void print() {
        fprintf(stderr, "{");
        for (int k = 0; k < n; k++) {
            atoms[k].print();
        }
        fprintf(stderr, "}\n");
        if (is_simple) {
            fprintf(stderr, "{");
            for (int k = 0; k < n; k++)
                fprintf(stderr, "%d,", simple[k]);
            fprintf(stderr, "}\n");
            fprintf(stderr, "{");
            for (int k = 0; k < N; k++)
                fprintf(stderr, "%d,", simple_binary[k]);
            fprintf(stderr, "}\n");
        }
    }

    partizione() {
        simple = 0;
        is_simple = false;
        n = 0;
        N = 0;
    }

};


double mylog[50];

double entropy_binary_partition(int *p, int n) {
    int i = 0;
    int mu;
    int begin;
    double H = 0;

    begin = 0;
    for (i = 1; i < n; i++) {
        if (p[i]) {
            mu = i - begin;
            begin = i;
            if (mu>1)
                H += (double)mu * mylog[mu];
        }
    }
    mu = n - begin;
    if (mu>1) H += mu * mylog[mu];

    H = -H / n + log((double)n);
    return (H);
}

double partizione::entropia_shannon() {
    double H = 0;
    int mu = 0, total, i;

    total = this->N;

    if (this->is_simple == false) {
        for (i = 0; i<this->n; i++) {
            mu = this->atoms[i].n;
            if (mu == 0) continue;
            H += mu * log(mu);
        }
        H = -H / total + log(total);
    } 
    
    if(this->is_simple)
        H = entropy_binary_partition(this->simple_binary, this->N);
    
    return (H);
}


typedef char* stringa;

class entry {
public:
    stringa name;
    stringa seq;
    int n;
    partizione Z;

    void partition();

    void print() {
        fprintf(stderr, "%s: %s\n", name, seq);
        Z.print();
        printf("Shannon %f, Topologica %d\n\n", Z.entropia_shannon(),
                Z.entropia_topologica());
    }

};

void entry::partition() {
    int i;
    int start;
    char a, b;

    start = 0;
    a = seq[0];

    //Total length of the partition is equal to the sequence
    Z.N = this->n;
    Z.atoms.resize(Z.N);

    //Divide the sequence in atoms of the same letter
    //When finding a different letter than we started with, we begin a new atom
    for (i = 1; i<this->n; i++) {
        b = seq[i];
        if (b != a) {
            this->Z.atoms.push_back(atom(start + 1, i - start, a));
            a = b;
            start = i;
        }
    }
    //taking care of the last atom, left hanging
    this->Z.atoms.push_back(atom(start + 1, i - start, a));

    //total number of atoms in the partition
    Z.n = Z.atoms.used;

    //simplified partition, only with the starting number of the atom in the sequence
    Z.simple = new int[Z.n];
    Z.simple_binary = new int[Z.N];
    for (i = 0; i < Z.N; i++)
        Z.simple_binary[i] = 0;
    for (i = 0; i < Z.n; i++) {
        Z.simple[i] = Z.atoms[i].start;
        Z.simple_binary[Z.simple[i] - 1] = 1;
    }
    Z.is_simple = true;


}

void print_binary_partition(int*p,int N){
    for(int i=0;i<N;i++)
        printf("%c", (p[i]) ? '|' : '.');
    printf("\n");
}

struct dist_tuple {
    double dist;
    double dist_r;
    int dist_top;
    int dist_top_r;
};

int *match = new int[1000];
int *reduced1 = new int[1000];
int *reduced2 = new int[1000];
int *intersezione = new int[1000];
int *intersezione_ridotta = new int[1000];

dist_tuple intersect(partizione &first, partizione &second) {
    dist_tuple res;
    int N=first.N;
    int dist_top=0;
    int dist_top_r=0;
    
    for (int i = 0; i < N; i++)
        match[i] = first.simple_binary[i] & second.simple_binary[i];
    match[0] = 1;

    for (int i = 0; i < N; i++)
        reduced1[i] = match[i] ^ first.simple_binary[i];
    reduced1[0] = 1;

    for (int i = 0; i < N; i++)
        reduced2[i] = match[i] ^ second.simple_binary[i];
    reduced2[0] = 1;

    for (int i = 0; i < N; i++)
        intersezione[i] =first.simple_binary[i] | second.simple_binary[i];
    intersezione[0] = 1;
    
    for (int i = 0; i < N; i++)
        intersezione_ridotta[i] = reduced1[i] | reduced2[i];
    intersezione_ridotta[0] = 1;
    
    for (int i=0;i<N;i++){
        dist_top+=intersezione[i];
        dist_top_r+=intersezione_ridotta[i];
    }
    
    double h1=entropy_binary_partition(first.simple_binary,N),
        h2=entropy_binary_partition(second.simple_binary,N),
        hr1=entropy_binary_partition(reduced1,N),
        hr2=entropy_binary_partition(reduced2,N),
        h12=entropy_binary_partition(intersezione,N),
        hr12=entropy_binary_partition(intersezione_ridotta,N);
    
    res.dist=2*h12-h1-h2;
    res.dist_r=2*hr12-hr1-hr2;
    res.dist_top=dist_top;
    res.dist_top_r=dist_top_r;
    
    
#ifdef DEBUG
    printf("H(1):\t\t%f\n",h1);
    print_binary_partition(first.simple_binary,N);
    printf("H(2):\t\t%f\n",h2);
    print_binary_partition(second.simple_binary,N);
//    printf("Hr(1):\t\t%f\n",hr1);
//    print_binary_partition(reduced1,N);
//    printf("Hr(2):\t\t%f\n",hr2);
//    print_binary_partition(reduced2,N);
    printf("H(12):\t\t%f\n",h12);
    print_binary_partition(intersezione,N);
    printf("Hr(12):\t\t%f\n",hr12);
    print_binary_partition(intersezione_ridotta,N);
    printf("\n");
    
    printf("d(1,2)=\t\t%f\n",res.dist);
    printf("dr(1,2)=\t%f\n\n\n",res.dist_r);
#endif
    
//    delete match; 
//    delete intersezione;
//    delete intersezione_ridotta;
//    delete reduced1;
//    delete reduced2;
    
    return(res);
}

entry new_entry(void) {
    entry temp;
    temp.name = 0;
    temp.seq = 0;
    temp.n = 0;
    return temp;
}




/*
 * 
 */
int main(int argc, char** argv) {

    int i,j;
    int seq_len=500;
    int n_seq=1000;
    double *adj;
    double *adj_r;
    FILE *out;
    dist_tuple res;
    easy_vect <entry> entries(n_seq, new_entry);
    
    
    srand(37337);
    for(i=0;i<50;i++)
        mylog[i]=log(i);
    
    
    for (i = 0; i < n_seq; i++) {
        char *name = new char[50];
        snprintf(name, 50, "seq%03d", i + 1);
        entries[i].name = name;
        char *seq = new char[seq_len+1];
        for (j = 0; j < seq_len; j++)
            seq[j] = rand() % 2 + 'A';
        seq[seq_len+1] = 0;
        entries[i].seq = seq;
        entries[i].n = seq_len;
    }


    for (i = 0; i < n_seq; i++) {
        entries[i].partition();
       // entries[i].print();   
    }
    
    adj=new double[n_seq*n_seq];
    adj_r=new double[n_seq*n_seq];
    
    for(i=0;i<n_seq*n_seq;i++){
        adj[i]=0;
        adj_r[i]=0;
    }
        
    for (i=0;i<n_seq;i++)
//       for(j=0;j<i+1;j++){
//            adj[i*n_seq+j]=0;
//            adj_r[i*n_seq+j]=1;
//       } 
        for(j=i+1;j<n_seq;j++){
            res=intersect(entries[i].Z,entries[j].Z);
            adj[i*n_seq+j]=res.dist;
            adj_r[i*n_seq+j]=res.dist_r;
        }
    
//    for (i=0;i<n_seq;i++){
//        for(j=i+1;j<n_seq;j++){
//            printf("%f ",adj[i*n_seq+j]);
//        }
//        printf("\n");
//    }
    
    out=fopen("output-adj.bin","w");
    i=fwrite(adj,sizeof(double),n_seq*n_seq,out);
    printf("Written %d records",i);
    
    out=fopen("output-adjr.bin","w");
    i=fwrite(adj_r,sizeof(double),n_seq*n_seq,out);
    printf("Written %d records",i);
    
    return 0;
}

