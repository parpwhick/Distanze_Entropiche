#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rand55.h"

#ifdef __i386
#define Norm_Factor   2.328306437080797e-10
#elif __amd64__
#define Norm_Factor   5.421010862427522170037264e-20 // 1/__UINT_MAX__ 
#endif

#define std_correction 1.414213562373095 // sqrt(2)

long Seed = 161803398L;

double rand55::rand() {
    n1 = n1->next;
    n2 = n2->next;
    n1->Y += n2->Y;
    return ((double) n1->Y * Norm_Factor);
}

unsigned long rand55::rand_long() {
    n1 = n1->next;
    n2 = n2->next;
    n1->Y += n2->Y;
    return (n1->Y);
}
//generatore di numeri gaussiani positivi

double rand55::semi_norm() {
    if (have_next_normal) {
        have_next_normal = 0;
        return (next_normal);
    }
    double u1 = 0, u2 = 0, s = 0;

    // s = R^2 = u1^2 + u2^2
    while (s >= 1 || s < 1e-15) {
        u1 = rand();
        u2 = rand();
        s = u1 * u1 + u2*u2;
    }

    double logs = std_correction * sqrt(-log(s) / s);

    u1 *= logs;
    u2 *= logs;

    have_next_normal = 1;
    next_normal = u2;
    return (u1);
}

void rand55::rand_init(long idum) {
    int i, j;
    long tmp, aux;
    double null;

    if ((Ran = new struct Lnk_List[55]) == NULL) { //(Lnk_List_Ptr) calloc(55, sizeof(struct Lnk_List))) == NULL)
        printf("Failed allocating memory for rand55\n");
        exit(1);
    }
    for (i = 0; i < 54; ++i) {
        (Ran + i)->next = Ran + i + 1;
    }
    (Ran + 54)->next = Ran;

    if (idum == -1) {
        FILE *out = fopen("/dev/urandom", "r");
        int bytes_read = fread(&idum, sizeof (long), 1, out);
        if (!bytes_read) {
            printf("Failed reading the seed, which wasn't provided\n");
            exit(1);
        }
        fclose(out);
    }

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
        null += rand();
    }
}

/* LONG RAND GENERATOR */


int p1 = P1, p2 = P2;
long table[SIZE];

long xrand() {
    int r;


    table[p1] = table[p1] + table[p2]; /* add two table elements */
    r = (table[p1] >> 1) & LONG_MAX; /* throw least significant bit away */


    if (p1 == SIZE1) { /* increment the table indexes */
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

void xrandinit(long seed) {
    int i;

    table[0] = seed;
    for (i = 1; i < SIZE; ++i)
        table[i] = (table[i - 1] * 1103515145) + Seed; /* lousy */

    for (i = 0; i < 10 * SIZE; ++i)
        xrand();
}
