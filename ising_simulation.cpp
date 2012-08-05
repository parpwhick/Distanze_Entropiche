#include <vector>
#include <cmath>
#include "adj_handler.h"
#include "strutture.h"
#include "rand_mersenne.h"
#include "ising_simulation.h"
#include "distance.h"
#include "smart_data_types.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using std::cout;
using std::endl;

#ifndef SOLO_SIMULAZIONE
extern
#endif
options opts;

using std::vector;
double *myexp=0;
#ifndef _GLIBCXX_DEBUG
template <typename T> inline void prefetch(T& element) {
    __builtin_prefetch(&element,0,0);
}
#else
#define prefetch(something); ;
#endif

std::vector<ising_simulation::config_t> ising_simulation::copy() {
    return config;
}

int generate_sierpinski_borders(int N, vector<int> &bordo_sinistro, vector<int> &bordo_destro, vector<int> &bordo_sotto){
    int gen=1;
    int size=6;
    int bordersize=2;

    // determine generation and proper dimensions
    while(size<N){
	gen+=1;
	size = 3*size-3;
	bordersize = 2*bordersize;
    }
    
    bordo_sinistro.resize(bordersize);
    bordo_destro.resize(bordersize);
    bordo_sotto.resize(bordersize);
    bordo_sinistro[0]=1;
    bordo_sinistro[1]=3;
    bordo_destro[0]=2;
    bordo_destro[1]=5;
    
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

int generate_square_border(int lato, vector<int> &bsx){
    bsx.resize(lato);
    
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
            if(!border1.empty())
                metropolis_subset(border1);
            if(!border2.empty())
                metropolis_subset(border2);
            if(!border3.empty())
                metropolis_subset(border3);
        }
}

double ising_simulation::energia_cinetica(){
    int totale = 0;
    for (int i = 0; i < NN.n_link; i++)
        totale += link_energies[i];
    return (totale+0.0)/NN.n_link;
}

double ising_simulation::magnetizzazione(){
    int totale = 0;
    for (int i = 0; i < N; i++)
        totale += config[i];
    double media = totale;
    media /= N;
    return (std::abs(media));
}


double ising_simulation::energia_magnetica() {
    int dH=0;
    // somma su i
    for (int i = 0; i < NN.n_link; i++)
        dH += - config[NN.adi[i]] * config[NN.adj[i]];
    dH /= 2;
    return (dH+0.0)/NN.n_link;
}

void ising_simulation::metropolis_step() {
    int s, z, somma_vicini, dH;

    if (config.empty()) {
        config.resize(NN.N);
        for(auto &x : config)
            x = 2 * (random.get_double() < 0.2) - 1 ;
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
        prefetch(config[s1]);
        prefetch(NN.index[s1]);
        
        z = NN.fetch(s);
        somma_vicini = 0;
        for (int m = 0; m < z; m++)
            somma_vicini += config[NN.vicini[m]];

        dH = config[s] * somma_vicini;
        if (dH <= 0 || random.get_double() < myexp[dH]) // exp(-2 * beta * dH))
            config[s] = -config[s];
    }
}

void ising_simulation::metropolis_subset(vector<int> subset) {
    int s, z, somma_vicini, dH;
      
    if (myexp == 0) {
        myexp = new double[NN.zmax + 2];
        for (int i = 0; i <= NN.zmax + 1; i++)
            myexp[i] = exp(-2 * beta * i);
    }
    
    /* Dinamica di Metropolis, 1 passo temporale */
    
    //update pari
    for (size_t j = 0; j < subset.size(); j+=2) {
        s = subset[j];
        z = NN.fetch(subset[j]);
        prefetch(config[subset[j+2]]);
        prefetch(NN.index[subset[j+2]]);
        somma_vicini = 0;
        for (int m = 0; m < z; m++)
            somma_vicini += config[NN.vicini[m]];

        dH = config[s] * somma_vicini;
        if (dH <= 0 || random.get_double() < myexp[dH])
            config[s] = -config[s];
    }
    
    //update dispari
    for (size_t j = 1; j < subset.size(); j+=2) {
        s = subset[j];
        z = NN.fetch(subset[j]);
        prefetch(config[subset[j+2]]);
        prefetch(NN.index[subset[j+2]]);
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
    if (config.empty()) {
        config.resize(NN.N);
        for (int i = 0; i < NN.N; i++)
            config[i] = 2 * (random.get_double() < 0.2) - 1 ;
    }
    if (link_energies.empty()) {
        link_energies.resize(NN.n_link);
        for (int i = 0; i < NN.n_link; i++)
            //distribuizione uniforme
            //link_energies[i] = random.get_int() % max_link_energy;
            
            //distribuzione esponenziale inversa
            //N = -1/lambda * ln(1-r) -1   poi da ceil'are.
            link_energies[i] = std::ceil(-1/opts.beta * std::log(1-random.get_double()) - 1);
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
        prefetch(NN.index[s1]);
        prefetch(NN.index[s2]);
        prefetch(config[s1]);
        prefetch(config[s2]);
        config_t &linkenergy = link_energies[link];        
        
        randnum = random.get_int();
        link = randnum % NN.n_link;
        __builtin_prefetch(&link_energies[link],1,0);
        prefetch(NN.adi[link]);
        prefetch(NN.adj[link]);
        
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
    config.resize(N);
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
    max_link_energy = 8;
    beta = 0.45;
    running = false;
    steps_per_time=time_length;
    skip=initial_time_skip;
    
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

template <typename data_t> void write_binary_array(const data_t *array, int N, const char *filename){
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

#ifndef SOLO_SIMULAZIONE
void time_series(const adj_struct &adj){
    general_partition Z1, Z2;
    distance dist(adj.N);
    auto_stats<double> d_top("Distanza_Topol"),d_top_r("Distanza_Top_r");
    auto_stats<double> d_shan("Distanza_Rohlin"),d_shan_r("Distanza_Rid");
    auto_stats<int> n_atomi("Numero_atomi");
    auto_stats<double> H("Entropia");
    auto_stats<double> mag("Magnetizzazione"), beta_est("Beta_stimato");
    auto_stats<double> E_kin("Energia_kin"), E_mag("Energia_mag");

    ising_simulation sim(adj,opts.simulation_type,opts.sweeps,opts.skip);
    sim.set_beta(opts.beta);
    sim.set_max_energy(opts.max_link_energy);
    if(opts.topologia == SIERPINSKI){
        generate_sierpinski_borders(adj.N,sim.border1,sim.border2,sim.border3);
    }        
    //sim.init_config();
    sim.step();
       
    Z1.from_configuration(sim.config_reference(),adj);
        
    /** Piccola nota sulla distribuzione delle energie:
     * se f(En) = exp(-beta * En) a meno di costante
     * allora <En> = 1 / (exp(beta) - 1) ====> beta = log[(1 + 1 / <En>]
     * delta(beta) = -1/[(<En>+1)*<En>)]  * sigma(En)/sqrt(N)
     */
    printf("%%t,\tatomi,\tH,\te_kin,\te_mag,\tdist,\tdist_r,\tM,\tbeta_est\n");
    
    for (int i = 1; i < opts.n_seq + 1; i++) {
        //calcolo quantita' da stampare
        E_kin = sim.energia_cinetica();
        E_mag = sim.energia_magnetica();
        mag = sim.magnetizzazione();
        beta_est = std::log(1. + 1. / E_kin);
        d_shan = dist.dist_shan;
        d_shan_r = dist.dist_shan_r;
        d_top = dist.dist_top;
        d_top_r = dist.dist_top_r;
        n_atomi = Z1.n;
        H = Z1.entropia_shannon;
        //stampa
        printf("%d\t", i);
        printf("%d\t", (int) n_atomi);
        printf("%.4f\t",(double) H);
        printf("%.4f\t",(double) E_kin);
        printf("%.4f\t",(double) E_mag);
        printf("%.4f\t",(double) d_shan);
        printf("%.4f\t",(double) d_shan_r);
        printf("%.4f\t",(double) mag);
        printf("%.4f\t",(double) beta_est);
        printf("\n");

        //simulation step
        sim.step();
        //aggiornamento partizione
        Z2.from_configuration(sim.config_reference(), adj);
        if (opts.graphics){
            write_binary_array(sim.config_reference(), adj.N, "configurazioni.bin");
            write_binary_array(Z2.show_labels(), adj.N, "partizioni.bin");
        }

        //calcolo distanze
        dist(Z1, Z2);
        //adesso Z1 conterra' la vecchia partizione, Z2 la prossima
        std::swap(Z1, Z2);
    }
    //stampa medie
    std::ofstream out1("medie.txt", std::ios::app);
    out1 << std::fixed << std::setprecision(4)
            << "%options: " << opts.command_line << endl
            << "%beta,\tatomi,\tH,\te_kin,\te_mag,\tdist,\tdist_r,\tdist_t,\tdist_tr\tM" << endl;
    out1 << beta_est.mean() << " \t"
            << n_atomi.mean() << "\t"
            << H.mean() << " \t"
            << E_kin.mean() << " \t"
            << E_mag.mean() << " \t"
            << d_shan.mean() << " \t"
            << d_shan_r.mean() << " \t"
            << d_top.mean() << " \t"
            << d_top_r.mean() << " \t"
            << mag.mean() << endl;
    //stampa varianze
    std::ofstream out2("varianze.txt", std::ios::app);
    out2 << std::fixed << std::setprecision(4)
            << "%options: " << opts.command_line << endl
            << "%beta,\tatomi,\tH,\te_kin,\te_mag,\tdist,\tdist_r,\tdist_t,\tdist_tr\tM" << endl;
    out2 << beta_est.var() << " \t"
            << n_atomi.var() << "\t"
            << H.var() << " \t"
            << E_kin.var() << " \t"
            << E_mag.var() << " \t"
            << d_shan.var() << " \t"
            << d_shan_r.var() << " \t"
            << d_top.var() << " \t"
            << d_top_r.var() << " \t"
            << mag.var() << endl;
    //le variabili auto_stat qui stampano le loro medie, lasciamo una riga di spazio
    printf("\n");
}
#endif


#ifdef SOLO_SIMULAZIONE
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

