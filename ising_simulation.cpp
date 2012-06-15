#include "adj_handler.h"
#include "strutture.h"
#include "rand_marsenne.h"
#include "rand55.h"
#include <cmath>


//Attenzione, generatore numeri casuali globale, bisogna prestare attenzione a 
//race conditions e rendere il tutto thread-safe innanzitutto!
//rand55 random55;
RandMT random;

void metropolis_step(int * &config, adj_struct NN, double beta) {
    static double *myexp = 0;
    int s,z,somma_vicini,dH;

    if (config == 0) {
        config = new int [NN.N];
        for (int i = 0; i < NN.N; i++)
            config[i] = random.get_int() % NN.N;
    }
    
    if (myexp == 0) {
        myexp = new double[NN.zmax + 2];
        for (int i = 0; i <= NN.zmax + 1; i++)
            myexp[i] = exp(-2 * beta * i);
    }

    /* Dinamica di Metropolis, 1 passo temporale */
    for (int j = 0; j < NN.N; j++) {
        // il generatore di numeri casuali influisce per un 5% sulla performance
        // totale, in realta' e' l'accesso disordinato alla memoria che 
        // uccide in confronto con s=j!
        s = random.get_int() % NN.N;
        z = NN.fetch(s);
        somma_vicini = 0;
        for (int m = 0; m < z; m++)
            somma_vicini += config[NN.vicini[m]];

        dH = config[s] * somma_vicini;
        if (dH <= 0 || random.get_float() < myexp[dH]) // exp(-2 * beta * dH))
            config[s] = -config[s];
    }
}

int max_link_energy=4;
void microcanonical_step(int * &config, int * &link_energies, adj_struct NN) {
    int dH;
    int somma_vicini;
    int z;
    
    //generazione dei buffer necessari, se vuoti
    if (config == 0) {
        config = new int [NN.N];
        for (int i = 0; i < NN.N; i++)
            config[i] = random.get_int() % NN.N;
    }    
    if (link_energies == 0) {
        link_energies = new int [NN.n_link];
        for (int i = 0; i < NN.n_link; i++)
            link_energies[i] = random.get_int() % max_link_energy;
    }
    
    /* Dinamica Microcanonica, 1 passo temporale */
    for (int j = 0; j < NN.N/2; j++) {
        // il generatore di numeri casuali influisce per un 5% sulla performance
        // totale, in realta' e' l'accesso disordinato alla memoria che 
        // uccide in confronto con s=j!
        //int s=generatore.rand_long() % N;
        uint32_t randnum = random.get_int();
        int link = randnum % NN.n_link;
        int accept1 = randnum & (1 << 29);
        int accept2 = randnum & (1 << 30);
        int accept12 = accept1 && accept2;

        //se nessuna mossa e' accettata, salta
        if (!(accept1 || accept2))
            continue;

        //creazione variabili locali, sperando di aumentare 
        //il pre-caricamento delle aree di memoria opportune
        //un 15% di miglioramento solo con la prima
        //   25% totale, cambiando le posizioni peggiora la performance!
        int &linkenergy = link_energies[link];
        int &s1 = NN.adi[link];
        int &chain1 = config[s1];
        int &s2 = NN.adj[link];
        int &chain2 = config[s2];

        //energia dai vicini di s1, se cambia                
        if (accept1) {
            z = NN.fetch(s1);
            somma_vicini = 0;
            for (int m = 0; m < z; m++)
                somma_vicini += config[NN.vicini[m]];
            dH = 2 * chain1 * somma_vicini;
        } else
            dH = 0;

        //energia dai vicini di s2, se cambia 
        if (accept2) {
            z = NN.fetch(s2);
            somma_vicini = 0;
            for (int m = 0; m < z; m++)
                somma_vicini += config[NN.vicini[m]];
            dH += 2 * chain2 * somma_vicini;
        }

        //energia dal flip simultaneo di s1 e s2
        if (accept12)
            dH -= 4 * chain1 * chain2;

        //caso dH==0 a parte, 5% di performance in piu
        if (dH == 0) {
            //flip di s1 o s2
            if (accept1) chain1 = -chain1;
            if (accept2) chain2 = -chain2;
        } else if (dH < 0 || linkenergy >= dH) {
            //dH < 0: l'energia del link e' incrementata
            //dH > 0: il link cede energia ai siti
            linkenergy -= dH;

            //flip di s1 o s2
            if (accept1) chain1 = -chain1;
            if (accept2) chain2 = -chain2;
        }

    }
}

void ising_simulation(adj_struct NN, int T=1000) { //, general_partition *partitions) {
    int *config=0, *link_energies=0;
    int N = NN.N;
        
/* inizializza i primi due terzi a +1
   l'ultimo terzo a -1  */    
    config = new int[N];
    int size=N/3+1;
    for (int i=0; i<size; i++)
        config[i]=+1;
    for (int i=size; i<2*size-2; i++)
        config[i]=+1;    
    for (int i=2*size-2; i<3*size-3; i++)
        config[i]=-1;
    double M1,M2,M3;

    for (int t=0; t< T; t++){
        
        M1=M2=M3=0.0;        
        for (int i=0; i<size; i++)
                M1+=(double)config[i];
        for (int i=size; i<2*size-2; i++)
                M2+=(double)config[i];  
        for (int i=2*size-2; i<3*size-3; i++)
                M3+=(double)config[i];
        M1/=size;
        M2/=size-2;
        M3/=size-1;
        
        printf("%d %f %f %f\n",t,M1,M2,M3);
        //printf("%2d/%d done, M: %g\n", t, T, (M_medio + 0.0) / (N + 0.0));
        //        partitions[i].from_square_lattice(chain, lato, 2);
        
        microcanonical_step(config,link_energies,NN);
        //si scartano i run con indice negativo, per termalizzare il sistema
        if (t < 0)
            continue;
    }
    delete[]config;

}

int main(int argc, char** argv) {
    int N = 0;
    int T=1000;
    adj_struct da_file = adiacenza_from_file("vector1.bin", "vector2.bin", N);
    //adj_struct square=adiacenza_square_lattice(200);
    //N=square.N;
    
    fprintf(stderr,"Loaded ADJ matrix\n\n");
    
    if(argc>1)
        T=atoi(argv[1]);
    fprintf(stderr,"Simulation %d steps long\n", T);
    if(argc>2)
        max_link_energy=atoi(argv[2]);
    fprintf(stderr,"Max link energy: %d\n", max_link_energy);
    fprintf(stderr,"\n");
    ising_simulation(da_file,T);    
}


