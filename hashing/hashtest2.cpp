#include <cstdio>
#include <cstdlib>


int main(){
	const int MAX=600;
	unsigned long int supertotal=0;
    for(int i=400; i <= MAX; i++){
	 unsigned long ihash=5381;
	 u_int8_t *word= (u_int8_t*) &i;
	 for(int m=0;m<4;m++)
		 ihash = (ihash << 5) + ihash + word[m];
	 for(int j=i+1; j <= MAX; j++){
		 u_int8_t *word= (u_int8_t*) &j;
		 unsigned long jhash=ihash;
		 for(int m=0;m<4;m++)
			 jhash = (jhash << 5) + jhash + word[m];
		 for(int k=j+1; k <= MAX; k++){
		 u_int8_t *word= (u_int8_t*) &k;
			 unsigned long khash=jhash;
			 for(int m=0;m<4;m++)
				 khash = (khash << 5) + khash + word[m];
			 printf("%lu\n",khash);
			 supertotal+=khash;
			 if(khash==193388366)
				 printf("multiple! %d %d %d\n",i,j,k);
		 }
	 }
	}
	printf("%lx supertotal\n",supertotal);
	return(0);
}
