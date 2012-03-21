#define Norm_Factor   2.328306437080797e-10

#ifdef MACRO
#define Ranf()  (((n1=n1->next)->Y += (n2=n2->next)->Y), (double)n1->Y * Norm_Factor)
#else
#define Ranf() rand55()
double rand55();
#endif

void rand_init(long);

typedef struct Lnk_List *Lnk_List_Ptr;
struct Lnk_List {
	unsigned long Y;
	Lnk_List_Ptr next;
} ;


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

