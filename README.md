Pagina principale        {#mainpage}
============

File importanti:
------------------
In general_distance.cpp tutti i file sono caricati, analizzati ed e' delegato al resto.
Il file piu' interessante e' partitions.cpp, in cui la classe general_partition e' contenuta - tutte le operazioni tra partizioni, altamente ottimizzate, sono li' definite. 
In ising_simulation.cpp sono definite le classi per le simulazioni con generiche adiacenze. 
In adj_handler.cpp i file per caricare e gestire le possibili topologie.


Come utilizzarlo?
------------------
Un paio di note:

I file in input sono in formato binario, intendo un file di interi a 32 bit, in cui ogni intero rappresenta un dato della configurazione, 
senza nessun spazio ne "a capo" tra le diverse configurazioni. In Matlab si ottiene molto semplicemente:
> fid = fopen('configurazioni.bin','wb')\n
> fwrite(fid,matricione,'int32')

L'output, se sono richieste le matrici delle distanze, sono matrici binarie di double:
-# output-shan.bin: distanza di Rohlin, usando entropia di Shannon, non ridotta
-# output-shan_r.bin: come sopra, ma con la riduzione specificata (partizione comune o eliminazione diretta, epsilon 0 o no)
-# output-top.bin: distanza di Rohlin, usando entropia topologica, non ridotta
-# output-top_r.bin: come sopra, con la riduzione specificata, anche se non ha senso per epsilon > 0
> fid = fopen('output-shan_r.bin','rb')\n
> dist=fread(fid,[dimensione,dimensione],'double');

Esempi:
> distanze_generiche -torus 100 -file input.bin -num 1000
Legge dal file "input.bin" 1000 configurazioni, crea le partizioni sapendo che rappresentano delle configurazioni su un toro 100x100, calcola le complete matrici delle distanze
Se il numero di configurazioni da leggere non e' specificato, controlla che il file "input.bin" abbia il numero giusto di dati, e ne carica tante quanto il file contiene.

> distanze_generiche -sequence 566 -file input.bin -num 1400 -epsilon 3
Legge stavolta 1400 configurazioni, sempre da un file binario, in cui ogni stato e' rappresentato da un intero, considerate come catene lunghe 566. Le distanze sono calcolate con riduzione "a meno di epsilon uguale o minore a 3".

> distanze_generiche -random -num 50 -symbols 2 -sierpinski 5 -common
Crea 50 configurazioni random, usando 2 simboli (equivalente ad un Ising), su un sierpinski di generazione 5, utilizzando come metodo di riduzione, quello del confronto con una partizione comune

> distanze_generiche -adj colonne.bin righe.bin -file configurazioni.bin -num 100
Carica le informazioni sull'adiacenza dai due file "colonne.bin" e "righe.bin", contenenti gli elementi nonnulli della matrice di adiacenza, legge 100 configurazioni e le trasforma in partizioni utilizzando le informazioni dell'adiacenza, e ne calcola tutte le distanze, usando come riduzione il metodo di default, la riduzione diretta tra atomi identici nelle due partizioni

> distanze_generiche -sierpinski 8 -microcanonical -num 5000 -sweeps 10 -skip 1000 -beta 0.44
Su un Sierpinski gasket di generazione 8, esegue 5000 iterazioni di Ising con la dinamica di tipo microcanonico, con 10 full sweep del reticolo per istante di tempo, saltando 1000 istanti di tempo all'inizio, calcolando le distanze tra una configurazione e la successiva, ad un beta (approssimativamente) di 0.44

> distanze_generiche -sierpinski 8 -microcanonical -num 5000 -sweeps 10 -skip 1000 -beta 0.44 -v -v 
Come prima, ma generando molto output, nei file "states.bin" e "partitions.bin" ci saranno 5000 configurazioni e partizioni rispettivamente, codificate come int32. Il file "output.txt" conterra' lo stato delle variabili iterazione per iterazione, "medie.txt" e "varianze.txt" le medie e varianze delle osservabili calcolate durante tutto il run del programma.

> distanze_generiche -sierpinski 8 -microcanonical -num 5000 -sweeps 10 -skip 1000 -beta 0.3,1.0 -v -v 
Come prima, solo che un lato del gasket ha beta pari a 0.3, l'altro 1.0

> distanze_generiche -sierpinski 8 -file states.bin -num 200
Calcola la matrici complete delle distanze tra le prime 200 configurazioni in output dal run precedente (ci mette davvero un attimo)

> distanze_generiche -demo -square 1000
Mostra statistiche sulla velocita' per operazioni eseguite su un reticolo con periodicita' cilindrica, 1000x1000



Opzioni generali:
-------------------------------------
opzione           |  cosa vuol dire
----------------- |  ---------------------------
  -random         |  Genera configurazioni a random, usando "symbols" simboli diversi
  -file FILENAME  |  Legge le configurazioni da un file, per poi calcolarne le matrici delle distanze
  -simulation     |  Genera una time series di configurazioni con la dinamica scelta
  -num N          |  Legge al massimo N / genera al massimo N configurazioni
  -nodistance     |  Non calcola le distanze di Rohlin
  -symbols N      |  Numero di simboli usati per la generazione random
  -nowrite        |  Non scrive le matrici distanze [non-default]
  -threads N      |  Usa N threads per il calcolo delle matrici distanze
  -v [-v -v]      |  -v: scrive i file output.txt, -v -v scrive gli stati e le partizioni, -v -v -v scrive anche le energie microcanoniche di ogni iterazione
  -help           |  Mostra questo messaggio
  -demo           |  Qualche test e benchmark del programma
  -graphics       |  Crea disegni in formato .ppm per tutte le operazioni svelte, se si usa la topologia quadrata

Opzioni per la riduzione:
-----------------
opzione           |  cosa vuol dire
----------------- |  ---------------------------
  -epsilon N      |  Usa epsilon riduzione con epsilon pari a N [default 0]
  -common         |  Riduce tramite confronto con partizione comune, invece che eliminando atomi "simili"

Scelta di una topologia:
--------------------------------
opzione           |  cosa vuol dire
----------------- |  ---------------------------
  -sequence N     | Sequenza lineare con 2 primi vicini lunga N
  -fuzzy NN       | In congiunzione con -sequence, i primi vicini diventano NN
  -torus L        | Reticolo bidimensionale periodico di lato L
  -square L       | Reticolo cilindrico, lato L
  -sierpinski G   | Sierpinski gasket generazione G
  -adj file1 file2| Topologia specificata a run time, tramite i due file, ad es. @c-adj rows.bin cols.bin


Opzioni per la simulazione:
-------------------
opzione           |  cosa vuol dire
----------------- |  ---------------------------
  -microcanonical | Dinamica icrocanonica
  -metropolis     | Dinamica di Metropolis
  -beta B[,B2,..] | Specifica la beta per ogni lato, oppure beta globale per Metropolis
  -sweeps N       | Numero di full lattice sweeps per ogni istante di tempo [default 1]
  -skip N         | Istanti di tempo da saltare inizialmente [default 50000]
  -suffix_out X   | Scrive file di output con suffisso X (eg. statesX.bin)
