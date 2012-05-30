#include <cstdlib>
#include <cstdio>
#include <vector>
int isSierpinskiCarpetPixelFilled2(int x, int y)
{
    while(x>0 || y>0) // when either of these reaches zero the pixel is determined to be on the edge 
                               // at that square level and must be filled
    {
        if(x %3==1 && y%3==1) //checks if the pixel is in the center for the current square level
            return 0;
        x /= 3; //x and y are decremented to check the next larger square level
        y /= 3;
    }
    return 1; // if all possible square levels are checked and the pixel is not determined 
                   // to be open it must be filled
}

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
int isSierpinskiCarpetPixelFilled3(int x, int y)
{
    while(x>0 || y>0) // when either of these reaches zero the pixel is determined to be on the edge 
                               // at that square level and must be filled
    {
        if(x % 2 ==1 && y % 2 == 1) //checks if the pixel is in the center for the current square level
            return 0;
        x /=2; //x and y are decremented to check the next larger square level
        y /=2;
    }
    return 1; // if all possible square levels are checked and the pixel is not determined 
                   // to be open it must be filled
}

memory<4> *vicini;


void addneigh(int sito, int vicino) {
    vicini[sito].add(vicino);
    vicini[vicino].add(sito);
}

inline int isSierpinskiCarpetPixelFilled(int x, int y)
{
    return((x & y) == 0);
}

int disegno=1;

int main(int argc, char** argv) {
    int potenza=12;
    if(argc>1){
        potenza=atoi(argv[1]);
        if(potenza<1){
            fprintf(stderr,"Input scorretto, deve essere un intero\n");
            return(2);
        }
        if(potenza>20){
            fprintf(stderr,"Specificata una generazione troppo grossa\n");
            return(1);
        }
    }
    long maxsite=1;
    long N=1;
    for (int i=0;i<potenza;i++){
        //triangolo
        maxsite *= 3;
        N *= 2;
        //carpet
        //maxsite *= 8;
        //N *= 3;
    }
    long filled=0;     
    int *line1 = new int[N];
    int *labels1 = new int[N];
    int *line2 = new int[N];
    int *labels2 = new int[N];
    int sitecount = 0;
       
    vicini = new memory<4>[maxsite];
    
    fprintf(stderr,"\nTriangolo Sierpinski, generazione %d: lato %ld, %ld siti nonnulli\n", potenza,N,maxsite);

    FILE* out;
    if(disegno){
    out=fopen("sierpinski.pbm","w");    
    fprintf(out,"P1\n");
    fprintf(out,"%ld %ld\n",N,N);
    fprintf(out,"\n");
    }
    FILE *vec1=fopen("vector1.bin","wb");
    FILE *vec2=fopen("vector2.bin","wb");
	
	if(!vec1 || !vec2)
		printf("Cant open files for output\n");

    for (int i = 0; i < N; i++) {
        // riempimento dei siti della riga corrente        
        line2[0] = 1;
        for (int j = 1; j < i + 1; j++)
            line2[j] = (line1[j] + line1[j - 1]) % 2;

        if(disegno)
            for (int j = 0; j< N; j++)
                fprintf(out,"%d",line2[j]);

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
        }
        if(disegno) fprintf(out,"\n");
        std::swap(line1, line2);
        std::swap(labels1,labels2);
    }
    if(disegno) fclose(out);
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
    fwrite(vectemp1,sizeof(int),scritti,vec1);
    fwrite(vectemp2,sizeof(int),scritti,vec2);
    fprintf(stderr,"%d elementi nonnulli della matrice di adiacenza\n",totale);
}
