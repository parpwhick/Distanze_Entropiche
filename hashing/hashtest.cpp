#include <stdio.h>
#include <cstdlib>
#include <stdint.h>

#include "GeneralHashFunctions.cpp"
#include "MurmurHash3.cpp"

const int MAX = 2500;

//unsigned int *arch=new unsigned int[MAX * (MAX-1) * (MAX-2) / 6];
int index=0;

int compare (const void * a, const void * b) {
  return ( *(unsigned int *)a - *(unsigned int*) b );
}

/* STUPID MURMURHASH !!!!!!!!!!!!!!!!1
 */
/*

FORCE_INLINE void fmix ( u_int32_t &h ) {
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
}
inline u_int32_t ROTL32 ( u_int32_t x, int8_t r )
{
  return (x << r) | (x >> (32 - r));
}
FORCE_INLINE void HASHSTEP(u_int32_t &h1, u_int32_t &k1){
    k1 *= 0xcc9e2d51;
    k1 = ROTL32(k1,15);
    k1 *= 0x1b873593;
    h1 ^= k1;
    h1 = ROTL32(h1,13); 
    h1 = h1*5+0xe6546b64; 
}
FORCE_INLINE void FINALIZEHASH(u_int32_t &h1, int len){
     h1 ^= len;
     fmix(h1);
}
u_int32_t MurmurHash3_x86_32 ( u_int32_t * blocks, int len, u_int32_t seed) {
    u_int32_t h1 = seed;
  for(int i = -len; i; i++)  {
    u_int32_t k1 = blocks[i];
    HASHSTEP(h1,k1);
  } 
  FINALIZEHASH(h1,len);  
  return(h1);
} 
*/

inline void HASHSTEP(u_int32_t &hash, u_int32_t value) {
    hash = (hash << 5) + hash + (0xff000000 & value);
    hash = (hash << 5) + hash + (0x00ff0000 & value);
    hash = (hash << 5) + hash + (0x0000ff00 & value);
    hash = (hash << 5) + hash + (0x000000ff & value);
}

int main() {
    
    unsigned long supertotal=0;
    unsigned int triplet[3];
    for (int i = 1; i <= MAX; i++) {
        
        
            for (int j = i + 1; j <= MAX; j++) {
            
            
    
            for (int k = j + 1; k <= MAX; k++) {
                unsigned int khash;
                triplet[0]=i;
                triplet[1]=j;
                triplet[2]=k;
                
                MurmurHash3_x86_32 ((unsigned char*)triplet,sizeof(int)*3,0,&khash);
                //khash = DJBHash ((unsigned char*)triplet,sizeof(int)*3);
//                khash=5381;
//                HASHSTEP(khash,i);
//                HASHSTEP(khash,j);
//                HASHSTEP(khash,k);
                           
                //arch[index++]=khash;
                supertotal+=khash;
                
                
            }
        }
    }
//    printf("index %d\n",index);
//    qsort(arch,index,sizeof(unsigned int),compare);
//    for(int i=1; i<index; i++){
//        if(arch[i]==arch[i-1])
//                printf("repeat! index %d, %x\n",i,arch[i]);
//        
//    }
    
    printf("%lx supertotal\n",supertotal);
    return (0);
}
