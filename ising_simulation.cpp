#include "adj_handler.h"
#include "strutture.h"
#include "rand_mersenne.h"
#include "ising_simulation.h"
#include <cmath>

#ifndef STANDALONE
extern
#endif
options opts;

double *myexp=0;

ising_simulation::config_t* ising_simulation::copy() {
    config_t *second = new config_t[N];
    for (int i = 0; i < N; i++)
        second[i] = config[i];
    return second;
}

int generate_sierpinski_borders(int N, int * &bordo_sinistro, int * &bordo_destro, int * &bordo_sotto){
    int gen=1;
    int size=6;
    int bordersize=2;

    // determine generation and proper dimensions
    while(size<N){
	gen+=1;
	size = 3*size-3;
	bordersize = 2*bordersize;
    }
    
    bordo_sinistro = new int[bordersize];
    bordo_destro = new int[bordersize];
    bordo_sotto = new int[bordersize];
    bordo_sinistro[0]=2;
    bordo_sinistro[1]=4;
    bordo_destro[0]=3;
    bordo_destro[1]=6;
    
    // we start over and populate
    size=6;
    bordersize=2;
    // iterate again over generations, to copy from the previous
    for(int g=2; g <= gen; g++){
	//copy left border, adding size of left triangle
	for(int i=0; i < bordersize; i++)
	    bordo_sinistro[bordersize+i]=bordo_sinistro[i]+size-1;
	//copy right border, adding size of right triangle
	for(int i=0; i < bordersize; i++)
	    bordo_destro[bordersize+i]=bordo_destro[i]+2*size-3;
	
	size = 3*size -3;
	bordersize = 2*bordersize;
    }
    for(int i=0; i< bordersize; i++)
        bordo_sotto[i] = N-1-i;
    return bordersize;
}

int generate_square_border(int lato, int * &bsx){
    bsx = new int[lato];    
    
    for(int i=0; i < lato; i++)
        bsx[i]=i;
    
    return(lato);
}

void ising_simulation::step(int steps){
    if(skip){
        steps+=skip;
        skip=0;
    }
    steps*=steps_per_time;
    
    if(update_rule==METROPOLIS)
        for(int i=0; i<steps; i++)
            metropolis_step();
    
    else if(update_rule==MICROCANONICAL)
        for(int i=0; i<steps; i++){
            microcanonical_step();
            if(border1)
                metropolis_subset(border1,border_size);
            if(border2)
                metropolis_subset(border2,border_size);
            if(border3)
                metropolis_subset(border2,border_size);
        }
}

int ising_simulation::energia_cinetica(){
    int totale = 0;
    for (int i = 0; i < NN.n_link; i++)
        totale += link_energies[i];
    return (totale);
}

double ising_simulation::magnetizzazione(){
    int totale = 0;
    for (int i = 0; i < N; i++)
        totale += config[i];
    double media = totale;
    media /= N;
    return (media);
}


int ising_simulation::energia_magnetica() {
    int dH=0;
    // somma su i
    for (int i = 0; i < NN.n_link; i++)
        dH += - config[NN.adi[i]] * config[NN.adj[i]];
    dH /= 2;
    return(dH);
}

void ising_simulation::metropolis_step() {
    int s, z, somma_vicini, dH;

    if (config == 0) {
        config = new config_t [NN.N];
        for (int i = 0; i < NN.N; i++)
            config[i] = 2 * (random.get_double() < 0.5) -1 ;
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

void ising_simulation::metropolis_subset(int *subset, int length) {
    int s, z, somma_vicini, dH;
      
    if (myexp == 0) {
        myexp = new double[NN.zmax + 2];
        for (int i = 0; i <= NN.zmax + 1; i++)
            myexp[i] = exp(-2 * beta * i);
    }
    
    /* Dinamica di Metropolis, 1 passo temporale */
    
    //update red
    for (int j = 0; j < length; j+=2) {        
        s = subset[j];
        z = NN.fetch(subset[j]);
        __builtin_prefetch(config+subset[j+2],0,0);
        __builtin_prefetch(NN.index+subset[j+2],0,0);
        somma_vicini = 0;
        for (int m = 0; m < z; m++)
            somma_vicini += config[NN.vicini[m]];

        dH = config[s] * somma_vicini;
        if (dH <= 0 || random.get_double() < myexp[dH])
            config[s] = -config[s];
    }
    
    //update black
    for (int j = 1; j < length; j+=2) {        
        s = subset[j];
        z = NN.fetch(subset[j]);
        __builtin_prefetch(config+subset[j+2],0,0);
        __builtin_prefetch(NN.index+subset[j+2],0,0);
        somma_vicini = 0;
        for (int m = 0; m < z; m++)
            somma_vicini += config[NN.vicini[m]];

        dH = config[s] * somma_vicini;
        if (dH <= 0 || random.get_double() < myexp[dH])
            config[s] = -config[s];
    }
}

void ising_simulation::microcanonical_step() {
    int dH;
    int somma_vicini;
    int z;

    //generazione dei buffer necessari, se vuoti
    if (config == 0) {
        config = new config_t [NN.N];
        for (int i = 0; i < NN.N; i++)
            config[i] = 2 * (random.get_double() < 0.5) - 1 ;
    }
    if (link_energies == 0) {
        link_energies = new config_t [NN.n_link];
        for (int i = 0; i < NN.n_link; i++)
            link_energies[i] = random.get_int() % max_link_energy;
    }

    uint32_t randnum = random.get_int();
    int link = randnum % NN.n_link;
    int s1,s2;
    /* Dinamica Microcanonica, 1 passo temporale */
    for (int j = 0; j < NN.N;) {
        // il generatore di numeri casuali influisce per un 5% sulla performance
        // totale, in realta' e' l'accesso disordinato alla memoria che uccide 
        
        int accept1 = (randnum & (1 << 29)) != 0;
        int accept2 = (randnum & (1 << 30)) != 0;
        int accept12 = accept1 && accept2;

        // adi e adj forniscono i due siti collegati dal 'link'
        s1 = NN.adi[link];
        s2 = NN.adj[link];
        
        /*
         L'accesso in memoria e' stato attentamente ottimizzato,
         precaricando i dati necessari all'iterazione corrente, poi generando 
         il link per l'iterazione successiva e caricando quello che dipende
         solo dal link prescelto.
         
         La differenza e' un 50% di tempo in piu' se si rimuovono i prefetch!
         */
        __builtin_prefetch(NN.index+s1,0,0);
        __builtin_prefetch(NN.index+s2,0,0);
        __builtin_prefetch(config+s1,0,0);
        __builtin_prefetch(config+s2,0,0);
        config_t &linkenergy = link_energies[link];        
        
        randnum = random.get_int();
        link = randnum % NN.n_link;
        __builtin_prefetch(link_energies+link,1,0);
        __builtin_prefetch(NN.adi+link,0,0);
        __builtin_prefetch(NN.adj+link,0,0);
        
        if (!(accept1 || accept2))
            continue;
	//se la mossa e' stata accettata, allora contala
	j++;

        //energia dai vicini di s1, se cambia                
        if (accept1) {
            // z numero dei vicini di s1
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
    config = new config_t[N];
    int size = N / 3 + 1;    
    for (int i = 0; i < size; i++)
        config[i] = +1;
    for (int i = size; i < 2 * size - 2; i++)
        config[i] = +1;
    for (int i = 2 * size - 2; i < 3 * size - 3; i++)
        config[i] = -1;
}

ising_simulation::ising_simulation(const adj_struct & NN1, simulation_t TT,
        int time_length,int initial_time_skip) : NN(NN1) {
    config = 0;
    link_energies = 0;
    max_link_energy = 8;
    beta = 0.45;
    running = false;
    steps_per_time=time_length;
    skip=initial_time_skip;
    border_size=0;
    border1=0;
    border2=0;
    border3=0;

    //NN = NN1;
    N = NN.N;

    update_rule = TT;
}

void ising_simulation::test_run(int T){
    init_config();
    measure();

    for (int t = 0; t < T; t++) {
        step();
        //si scartano i run con indice negativo, per termalizzare il sistema
        if (t < 0)
            continue;
        measure();
    }
    measure();
}

template <typename data_t> void write_binary_array(data_t *array, int N, const char *filename){
    FILE *out;

    if(filename==0)
        return;
    
    out = fopen(filename, "ab");
    if (out == 0) {
        fprintf(stderr, "Error opening file for writing\n");
        exit(1);
    }
    
    fwrite(array,sizeof(data_t),N,out);
}

#ifndef STANDALONE
void time_series(const adj_struct &adj){
    general_partition Z1, Z2;
    distance d(adj.N);
    int E_kin=0, E_mag=0;
    double mag;
    
    ising_simulation sim(adj,opts.simulation_type,3,0);
    sim.set_beta(opts.beta);
    sim.set_max_energy(opts.max_link_energy);
    if(opts.topologia == SIERPINSKI){
        sim.border_size=generate_sierpinski_borders(adj.N,sim.border1,sim.border2,sim.border3);
    }
        
    //sim.init_config();
    sim.step();
       
    Z1.from_configuration(sim.config_reference(),adj);
    E_kin = sim.energia_cinetica();
    E_mag = sim.energia_magnetica();
    mag = sim.magnetizzazione();
    
    printf("%%t\tatomi\tentropia\te kin\te mag\tdist\t\tdist_ridotta\tmagnetizz\n");
    printf("%d\t%d\t%.6f\t%d\t%d\t%.6f\t%.6f\t%.6f\n",1,Z1.n,Z1.entropia_shannon,E_kin, E_mag, 0.0, 0.0,mag);
    for(int i=1; i<opts.n_seq; i++){
        sim.step();
        if(opts.graphics)
            write_binary_array(sim.config_reference(),adj.N, "configurazioni.bin");
        if(i%2){
            Z2.from_configuration(sim.config_reference(),adj);
            if(opts.graphics)
            write_binary_array(&Z2.labels[0],adj.N,"partizioni.bin");
        }
        else
            Z1.from_configuration(sim.config_reference(),adj);
        d.fill(Z1,Z2);
        E_kin = sim.energia_cinetica();
        E_mag = sim.energia_magnetica();
        mag = sim.magnetizzazione();
        printf("%d\t%d\t%.6f\t%d\t%d\t%.6f\t%.6f\t%.6f\n",i+1,Z2.n,Z2.entropia_shannon,E_kin,E_mag,
                                d.dist_shan,d.dist_shan_r,mag);
    }
}
#endif


#ifdef STANDALONE
int main(int argc, char** argv) {
//    int N = 0;
    int T = 1000;
    int max_link_energy=4;
    adj_struct da_file = adiacenza_square_lattice(50);

    if (argc > 1)
        T = atoi(argv[1]);
    fprintf(stderr, "Simulation %d steps long\n", T);
    if (argc > 2)
        max_link_energy = atoi(argv[2]);
    fprintf(stderr, "Max link energy: %d\n", max_link_energy);
    fprintf(stderr, "\n");
    ising_simulation sim(da_file, MICROCANONICAL,1,0);
    sim.set_max_energy(max_link_energy);
    sim.test_run(T);
}
#endif

