/*
  Date: 06/03/06 11:33
  Description: Generazione matrice di adiacenza per gasket di sierpinski a generazione g = GEN. 

  Versione 2, 31/05/2013, Dawid Crivelli
  -update in place
  -utilizzo memoria ottimale
*/


#include <cstdio>
#include <cstdlib>
#include <cstring>

template <int cols>
class memory {
    int *buffer;
    int nrows;
    int total;
public:

    memory(int rows = 0, int fill=0) {
        total = rows*cols;
        nrows = rows;
        if (nrows) {
            fprintf(stderr, "Requested allocation space of %u kbytes\n", ((sizeof (int)) * total) >> 10);
            buffer = new int[total];
            if(fill) 
                for (int i=0; i<total ;i++)
                        buffer[i]=fill;
        }
    }

    int & operator()(int row, int col) {
        return buffer[cols * row + col];
    }

    int & operator[](int row) {
        return buffer[row];
    }

    int & operator()(int row) {
        return buffer[cols * row];
    }

    ~memory() {
        if(nrows)
                delete buffer;
    }
};

int GEN = 2;

int main(int argc, char **argv) {
    int a0, a1, a2;
    int size, old_size, total_size;

    if (argc > 1) {
        GEN = atoi(argv[1]);
        if (GEN < 1) {
            fprintf(stderr, "Input scorretto, deve essere un intero\n");
            return (2);
        }
        if (GEN > 20) {
            fprintf(stderr, "Specificata una generazione troppo grossa\n");
            return (1);
        }
    }

    total_size = 6;
    for (int g = 2; g <= GEN; g++)
        total_size = 3 * total_size - 3;
    fprintf(stderr, "Generazione %d, size %d\n", GEN, total_size);

    // prima generazione, dimensione 6
    old_size = 6;
    size = 6;

    /*allocazione */
    
    // vettore di coordinazione
    memory < 1 > z  = *new memory < 1 > (total_size);
    // matrice (n,4) nearest neighbour
    memory < 4 > nn = *new memory < 4 > (total_size, -1);
    
    

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


    } //chiudo ciclo g

    /* Creazione vettori di adiacenza per matrice sparse in stile Matlab */
    /* La funzione per caricare i dati cosi creati e':
    -------------------------------
    function adiacenza=load_sierpinski()
	indici_riga=fread(fopen('vector1.bin','r'),inf,'int32');
	indici_colonna=fread(fopen('vector2.bin','r'),inf,'int32');
	N=max(max(indici_riga),max(indici_colonna));
	adiacenza=sparse(indici_riga,indici_colonna,1,N,N);
    end
    ******************************/

    int write = 1;
    if (write) {
        // scrittura in blocchi di lunghezza chunk, per alleggerire le operazioni
        int chunk = 500000;

        FILE *vec1 = fopen("vector1.bin", "wb");
        FILE *vec2 = fopen("vector2.bin", "wb");

        int *vectemp1 = new int[11 * chunk / 10];
        int *vectemp2 = new int[11 * chunk / 10];
        int scritti = 0;
        int totale = 0;
        for (int i = 0; i < size; i++) {

            for (int j = 0; j < z(i); j++) {
                vectemp1[scritti + j] = i + 1;
                vectemp2[scritti + j] = nn(i, j) + 1;
            }
            scritti += z(i);
            totale += z(i);
            if (scritti > chunk) {
                fwrite(vectemp1, sizeof (int), scritti, vec1);
                fwrite(vectemp2, sizeof (int), scritti, vec2);
                scritti = 0;
            }
        }
        fwrite(vectemp1, sizeof (int), scritti, vec1);
        fwrite(vectemp2, sizeof (int), scritti, vec2);
        fprintf(stderr, "%d elementi nonnulli della matrice di adiacenza\n", totale);
    }


}




