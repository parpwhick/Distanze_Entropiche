#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cassert>
#include <map>

template <int max>
class memory {
public:
	int quanti;
	int elementi[max];
        
        memory(){
            quanti=0;
        }
        
        void add(int i){
            if (quanti>=max)
                return;
            elementi[quanti]=i;
            quanti++;
        }
};

memory<3> *vicini;


void addneigh(int sito, int vicino) {
    vicini[sito].add(vicino);
    vicini[vicino].add(sito);
}

int main(int argc, char** argv) {
    int potenza=12;
    if(argc>1){
        potenza=atoi(argv[1]);
        if(potenza<1){
            fprintf(stderr,"Input scorretto, deve essere un intero\n");
            return(2);
        }
        if(potenza>16){
            fprintf(stderr,"Specificata una generazione troppo grossa\n");
            return(1);
        }
    }
    int N = 1 << potenza; // lunghezza del lato del triangolo piu grosso = 2^potenza 
    int *line1 = new int[N];
    int *labels1 = new int[N];
    int *line2 = new int[N];
    int *labels2 = new int[N];
    int sitecount = 0;
    
    int maxsite=1;
    for (int i=0;i<potenza;i++)
        maxsite *=3;    
    vicini = new memory<3>[maxsite];
    
    fprintf(stderr,"\nTriangolo Sierpinski, generazione %d: %d siti nonnulli\n", potenza,maxsite);

    //    FILE* out=fopen("sierpinski.pbm","w");    
    //    fprintf(out,"P1\n");
    //    fprintf(out,"%d %d\n",N,N);
    //    fprintf(out,"\n");
    
    FILE *vec1=fopen("vector1.bin","w");
    FILE *vec2=fopen("vector2.bin","w");

    for (int i = 0; i < N; i++)
        line1[i] = line2[i] = 0;
    line1[0] = 1;

    for (int i = 0; i < N; i++) {
        // riempimento dei siti della riga corrente
        line2[0] = 1;
        for (int j = 1; j < i + 1; j++)
            line2[j] = (line1[j] + line1[j - 1]) % 2;

        // labelling dei siti nonnulli
        for (int j = 0; j < i + 1; j++) {
            if (line2[j]) {
                labels2[j] = sitecount;
                sitecount++;

                if(i==0)
                    continue;
                //aggiungi all'elenco dei vicini
                if (line1[j])
                    addneigh(labels2[j], labels1[j]);
                if (j > 0) {
                    if (line1[j - 1])
                        addneigh(labels2[j], labels1[j - 1]);
                    if (line2[j - 1])
                        addneigh(labels2[j], labels2[j - 1]);
                }
            } else {
                labels2[j] = 0;
            }
            //fprintf(out,"%d",line2[j]);
        }
        //fprintf(out,"\n");
        std::swap(line1, line2);
        std::swap(labels1,labels2);
    }
    
    fprintf(stderr,"Scrittura vettori sparse\n");
    int chunk=500000;
    
    int *vectemp1=new int[11*chunk/10];
    int *vectemp2=new int[11*chunk/10];
    int scritti=0;
    int totale=0;
    for (int i=0; i<sitecount; i++){        
                
        for(int j=0; j<vicini[i].quanti; j++){
            vectemp1[scritti+j]=i+1;
            vectemp2[scritti+j]=vicini[i].elementi[j]+1;
        }        
        scritti+=vicini[i].quanti;
        totale+=vicini[i].quanti;
        if(scritti > chunk){            
               fwrite(vectemp1,sizeof(int),scritti,vec1);
               fwrite(vectemp2,sizeof(int),scritti,vec2);
               scritti=0;
        }
    }
    fprintf(stderr,"%d elementi nonnulli della matrice di adiacenza\n",totale);
    
}
