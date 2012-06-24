#include "adj_handler.h"
#include "strutture.h"
#include "rand_mersenne.h"
#include "rand55.h"
#include <cmath>

options opts;

class ising_simulation {
    simulation_type update_rule;

    adj_struct NN;
    int N;
    int T;
    int *config;
    int *link_energies;
    int max_link_energy;
    bool running;
    double beta;
    RandMT random;

    void metropolis_step();
    void microcanonical_step();
public:
    ising_simulation(adj_struct NN1, simulation_type TT);
    void measure();
    void step(int steps = 1);
    void init_config();
    int *copy();
};

int* ising_simulation::copy() {
    int *second = new int[N];
    for (int i = 0; i < N; i++)
        second[i] = config[i];
    return second;
}

void ising_simulation::step(int steps){
    if(update_rule==METROPOLIS)
        for(int i=0; i<steps; i++)
            metropolis_step();
    else if(update_rule==MICROCANONICAL)
        for(int i=0; i<steps; i++)
            microcanonical_step();
}

void ising_simulation::metropolis_step() {
    static double *myexp = 0;
    int s, z, somma_vicini, dH;

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

    int s1 = random.get_int() % NN.N;    
    /* Dinamica di Metropolis, 1 passo temporale */
    for (int j = 0; j < NN.N; j++) {
        // il generatore di numeri casuali influisce per un 5% sulla performance
        // totale, in realta' e' l'accesso disordinato alla memoria che 
        // uccide in confronto con s=j!
        s=s1;
        s1 = random.get_int() % NN.N;
        __builtin_prefetch(config+s1,0,0);
        __builtin_prefetch(NN.index+s1,0,0);        
        
        z = NN.fetch(s);
        somma_vicini = 0;
        for (int m = 0; m < z; m++)
            somma_vicini += config[NN.vicini[m]];

        dH = config[s] * somma_vicini;
        if (dH <= 0 || random.get_double() < myexp[dH]) // exp(-2 * beta * dH))
            config[s] = -config[s];
    }
}

void ising_simulation::microcanonical_step() {
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

    uint32_t randnum = random.get_int();
    int link = randnum % NN.n_link;
    int s1,s2;
    /* Dinamica Microcanonica, 1 passo temporale */
    for (int j = 0; j < NN.N / 2; j++) {
        // il generatore di numeri casuali influisce per un 5% sulla performance
        // totale, in realta' e' l'accesso disordinato alla memoria che uccide 
        
        int accept1 = randnum & (1 << 29);
        int accept2 = randnum & (1 << 30);
        int accept12 = accept1 && accept2;

        s1 = NN.adi[link];
        s2 = NN.adj[link];
        __builtin_prefetch(NN.index+s1,0,0);
        __builtin_prefetch(NN.index+s2,0,0);
        __builtin_prefetch(config+s1,0,0);
        __builtin_prefetch(config+s2,0,0);
        //creazione variabili locali, sperando di aumentare 
        //il pre-caricamento delle aree di memoria opportune
        //un 15% di miglioramento solo con la prima
        //   25% totale, cambiando le posizioni peggiora la performance!
        int &linkenergy = link_energies[link];        
        
        randnum = random.get_int();
        link = randnum % NN.n_link;
        __builtin_prefetch(link_energies+link,1,0);
        __builtin_prefetch(NN.adi+link,0,0);
        __builtin_prefetch(NN.adj+link,0,0);
        
        if (!(accept1 || accept2))
            continue;

        //energia dai vicini di s1, se cambia                
        if (accept1) {
            z = NN.fetch(s1);
            somma_vicini = 0;
            for (int m = 0; m < z; m++)
                somma_vicini += config[NN.vicini[m]];
            dH = 2 * config[s1] * somma_vicini;
        } else
            dH = 0;

        //energia dai vicini di s2, se cambia 
        if (accept2) {
            z = NN.fetch(s2);
            somma_vicini = 0;
            for (int m = 0; m < z; m++)
                somma_vicini += config[NN.vicini[m]];
            dH += 2 * config[s2] * somma_vicini;
        }
        
        //energia dal flip simultaneo di s1 e s2
        if (accept12)
            dH -= 4 * config[s1] * config[s2];

        //caso dH==0 a parte, 5% di performance in piu
        if (dH == 0) {
            //flip di s1 o s2
            if (accept1) config[s1] = -config[s1];
            if (accept2) config[s2] = -config[s2];
        } else if (dH < 0 || linkenergy >= dH) {
            //dH < 0: l'energia del link e' incrementata
            //dH > 0: il link cede energia ai siti
            linkenergy -= dH;

            //flip di s1 o s2
            if (accept1) config[s1] = -config[s1];
            if (accept2) config[s2] = -config[s2];
        }

    }
}

void ising_simulation::measure() {
    static int t = 0;
    double M1 = 0, M2 = 0, M3 = 0;
    int size = N / 3 + 1;

    for (int i = 0; i < size; i++)
        M1 += (double) config[i];
    for (int i = size; i < 2 * size - 2; i++)
        M2 += (double) config[i];
    for (int i = 2 * size - 2; i < 3 * size - 3; i++)
        M3 += (double) config[i];
    M1 /= size;
    M2 /= size - 2;
    M3 /= size - 1;

    printf("%d %f %f %f\n", t, M1, M2, M3);
    //printf("%2d/%d done, M: %g\n", t, T, (M_medio + 0.0) / (N + 0.0));
    t++;
}

void ising_simulation::init_config() {
    config = new int[N];
    int size = N / 3 + 1;
    for (int i = 0; i < size; i++)
        config[i] = +1;
    for (int i = size; i < 2 * size - 2; i++)
        config[i] = +1;
    for (int i = 2 * size - 2; i < 3 * size - 3; i++)
        config[i] = -1;
}

ising_simulation::ising_simulation(adj_struct NN1, simulation_type TT) {
    config = 0;
    link_energies = 0;
    max_link_energy = 8;
    running = false;
    T=1000;

    NN = NN1;
    N = NN.N;

    update_rule = TT;

    init_config();
    measure();

    for (int t = 0; t < T; t++) {
        //microcanonical_step();
        metropolis_step();
        //si scartano i run con indice negativo, per termalizzare il sistema
        if (t < 0)
            continue;
        //measure();
    }
    measure();
    delete[]config;

}

int main(int argc, char** argv) {
    int N = 0;
    int T = 1000;
    int max_link_energy=4;
    adj_struct da_file = adiacenza_from_file("vector1.bin", "vector2.bin", N);
    //adj_struct square=adiacenza_square_lattice(200);
    //N=square.N;

    fprintf(stderr, "Loaded ADJ matrix\n\n");

    if (argc > 1)
        T = atoi(argv[1]);
    fprintf(stderr, "Simulation %d steps long\n", T);
    if (argc > 2)
        max_link_energy = atoi(argv[2]);
    fprintf(stderr, "Max link energy: %d\n", max_link_energy);
    fprintf(stderr, "\n");
    ising_simulation run(da_file, MICROCANONICAL);
}


