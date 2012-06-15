#ifndef RAND55_H
#define	RAND55_H


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
    unsigned long rand_long();
    unsigned long get_int(){ return rand_long();}
    double get_float(){ return rand();}
    double semi_norm() ;
    
    rand55(long idum=-1) { rand_init(idum); have_next_normal=0;}
    ~rand55() { delete[]Ran;}
    
    double operator()(){
        return rand();
    }
      
    
};

#endif