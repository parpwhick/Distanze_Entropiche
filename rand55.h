typedef struct Lnk_List *Lnk_List_Ptr;

struct Lnk_List {
	unsigned long Y;
	Lnk_List_Ptr next;
} ;




class rand55 {
    Lnk_List_Ptr Ran, n1, n2;
    void rand_init(long idum);
    
    int have_next_normal;
    double next_normal;
    
public:
    double rand() ;
    double semi_norm() ;
    
    rand55(long idum=-1) { rand_init(idum); have_next_normal=0;}
    ~rand55() { delete[]Ran;}
    
    double operator()(){
        return rand();
    }
      
    
};



#ifdef LONG31 /* x^31 + x^3 + 1 */
#define SIZE 31
#define SIZE1 30
#define P1 3
#define P2 0
#else /* LONG63: x^63 + x + 1 */
#define SIZE 63
#define SIZE1 62
#define P1 1
#define P2 0
#endif


#define LONG_MAX 0x7fffffff

long xrand();
void xrandinit(long seed);

