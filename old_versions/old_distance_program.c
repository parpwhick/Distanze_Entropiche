#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

int levenshtein(char *s,char*t);
int minimum(int a,int b,int c);

#define min(a,b)  {(a<b)?a:b}

/****************************************/
/*Implementation of Levenshtein distance*/
/****************************************/

int levenshtein(char *s,char*t)
/*Compute levenshtein distance between s and t*/
{
  //Step 1
  int k,i,j,n,m,cost,*d,distance;
  n=strlen(s); 
  m=strlen(t);
  if(n!=0&&m!=0)
  {
    d=malloc((sizeof(int))*(m+1)*(n+1));
    m++;
    n++;
    //Step 2	
    for(k=0;k<n;k++)
	d[k]=k;
    for(k=0;k<m;k++)
      d[k*n]=k;
    //Step 3 and 4	
    for(i=1;i<n;i++)
      for(j=1;j<m;j++)
	{
        //Step 5
        if(s[i-1]==t[j-1])
          cost=0;
        else
          cost=1;
        //Step 6			 
        d[j*n+i]=minimum(d[(j-1)*n+i]+1,d[j*n+i-1]+1,d[(j-1)*n+i-1]+cost);
      }
    distance=d[n*m-1];
    free(d);
    return distance;
  }
  else 
    return -1; //a negative return value means that one or both strings are empty.
}

int minimum(int a,int b,int c)
/*Gets the minimum of three values*/
{
  int min=a;
  if(b<min)
    min=b;
  if(c<min)
    min=c;
  return min;
}

typedef struct {
	char *name;
	int n;
	char *seq;
} entry;

typedef struct {
	int *punti;
	int n;
} atom;

typedef struct {
	atom *atomi;
	int n;
	int cumsum;
} partizione;

atom a_temp; 
atom intersect(atom a,atom b) {
	int shorter,longer;
	int i,j,k;
	int *primo,*secondo;
	atom def;
	if(a.n < b.n) {
		shorter=a.n;
		longer=b.n;
		primo=a.punti;
		secondo=b.punti;
	} else {
		shorter=b.n;
		longer=a.n;
		primo=b.punti;
		secondo=a.punti;
	}
	if(shorter>a_temp.n){
		free(a_temp.punti);
		a_temp.n=shorter;
	 	a_temp.punti=malloc(sizeof(int)*shorter);
	}
	// check every member of the shorter
	i=0;
	for(k=0;k<shorter;k++){
		//if it matches any member of the longer
		for(j=0;j<longer;j++){
			if (a.punti[k]==b.punti[j]){
				a_temp.punti[i++]=a.punti[k];
				break;
			}
		}
	}
	def.n=i;
	def.punti=malloc(sizeof(int)*i);
	memcpy(def.punti,a_temp.punti,i*sizeof(int));
	return def;
}

int printatom(atom a) {
	int i;
	printf("atomo {");
	for(i=0;i<a.n;i++)
		printf("%d, ",a.punti[i]);
	printf("}\n");
	return i;
}


double entropia(partizione a){
	int i;
	int n=0;
	double h=0;
	for(i=0;i<a.n;i++){	
		n+=a.atomi->n;
		h+=(a.atomi->n)*log(a.atomi->n);
	}
	h=-h/n+log(n);
}

int* findmatches(entry *elenco,int max){
	int k,i,j;
	int *distanze=calloc(sizeof(int),max*max);
	k=1;
	fprintf(stderr,"Computing distance\n");
	for(i=0;i<max;i++) {
		for(j=i+1;j<max;j++) 
		//	printf("(%d,%d) ",i,j);
			distanze[j*max+i]=
				levenshtein(elenco[i].seq,elenco[j].seq);
//				strlen(elenco[j].seq);
		if((i) / max > k) {
			fprintf(stderr,".");
			fprintf(stderr," %d%% ",k);
			k++;
		}
		for(j=0;j<max;j++) 
			printf("%3d ",distanze[j*max+i]);
		printf("\n");

	}
	return distanze;
}

int main(int argc,char **argv) {
	char buf[32768];
	int i=0,i_atom=0,i_punto;
	int max;
	atom input;
	if (argc > 1) {
		max=atoi(argv[1]);
	}
	else {
		fprintf(stderr,"Setting 512 max input lines\n");
		max=512;
	}

//	a_temp.n=50;
//	a_temp.punti=malloc(sizeof(int)*a_temp.n);
//	partizione.n=50;
//	partizione.atomi=malloc(sizeof(atom)*partizione.n);
	
	/* READ INPUT STRINGS WITH A NAME ATTACHED
	 *
	 *
	 */
	char namebuf[256];
	char seqbuf[32768];
	entry *elenco;
	entry temp;
	elenco=malloc(sizeof(entry)*max);
	while(!(feof(stdin)) && scanf("%s %s",namebuf,seqbuf) && i < max){
		temp.seq=strdup(seqbuf);
		temp.name=strdup(namebuf);
		temp.n=i;
		elenco[i]=temp;
		i++;
	}
	max=i-1;
	fprintf(stderr,"read %d strings\n",max);
	findmatches(elenco,max);
/*	for(i=0;i<max;i++)
		printf("%d: %s %s\n",elenco[i].n,
				elenco[i].name,elenco[i].seq);
*/

		return 0;
}

