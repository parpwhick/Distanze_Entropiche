#include <stdio.h>
#include <stdlib.h>
#include "rand55.h"

long Seed = 161803398L;

Lnk_List_Ptr Ran, n1, n2;

//inline 
double rand55()
{
    n1 = n1->next;
    n2 = n2->next;
    n1->Y += n2->Y;

    return (double) n1->Y * Norm_Factor;
}

void rand_init(long idum)
{
    int i, j;
    long tmp, aux;
    double null;

    if ((Ran = (Lnk_List_Ptr) calloc(55, sizeof(struct Lnk_List))) == NULL)
	exit(1);
    for (i = 0; i < 54; ++i) {
	(Ran + i)->next = Ran + i + 1;
    }
    (Ran + 54)->next = Ran;
    tmp = Seed + idum;
    Ran[54].Y = tmp;
    aux = 1;
    j = 20;
    for (i = 0; i < 54; ++i, j += 21) {
	j %= 55;
	Ran[j].Y = aux;
	aux += tmp;
	tmp = Ran[j].Y;
    }

    n1 = Ran;
    n2 = Ran + 31;

    for (j = 0; j < 55 * 4; ++j) {
	null += rand55();
    }

}

/* LONG RAND GENERATOR */


int p1 = P1, p2 = P2;
long table[SIZE];

long xrand()
{
    int r;


    table[p1] = table[p1] + table[p2];	/* add two table elements */
    r = (table[p1] >> 1) & LONG_MAX;	/* throw least significant bit away */


    if (p1 == SIZE1) {		/* increment the table indexes */
	p1 = 0;
	p2 = p2 + 1;
    } else if (p2 == SIZE1) {
	p1 = p1 + 1;
	p2 = 0;
    } else {
	p1 = p1 + 1;
	p2 = p2 + 1;
    }


    return (r);
}




void xrandinit(long seed)
{
    int i;


    table[0] = seed;
    for (i = 1; i < SIZE; ++i)
	table[i] = (table[i - 1] * 1103515145) + Seed;	/* lousy */


    for (i = 0; i < 10 * SIZE; ++i)
	xrand();
}




/*** a small test program ***/
