#include <cstdlib>
#include <cstdio>
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
    fprintf(stderr,"Sierpinski carpet, generazione %d, lato %ld, siti totali %ld\n",potenza,N,maxsite);
    int *line1 = new int[N];
    
    //int *cache=cachesearch(N);
        
    FILE* out=0;
    if(disegno){
    out=fopen("carpet.pbm","wb");
    fprintf(out,"P1\n");
    fprintf(out,"%ld %ld\n",N,N);
    fprintf(out,"\n");
    }

    // righe
    for (int i = 0; i < N; i++) {
        for (int j=0; j < N-i; j++){
            line1[j] =  isSierpinskiCarpetPixelFilled(i,j);
	    filled+=line1[j];
	}
	if(disegno) {
		for (int j=0; j < N; j++)
			fprintf(out,"%d ",line1[j]);
		fprintf(out,"\n");
	}
    }
    if(disegno) fclose(out); 
    
    fprintf(stderr,"Filled %ld/%ld (%g%%)\n",filled,N*N,filled/((N+0.0)*N)*100);
    return(0);
    
}
