/* 
 * File:   matrix.h
 * Author: fake
 *
 * Created on June 27, 2012, 3:55 PM
 */

#ifndef MATRIX_H
#define	MATRIX_H

#include <cmath>
#include <vector>
#include <cstdio>
#include <stdexcept>
#include <sstream>

///Funziona come una normale variabile, ma con nome, statistiche interne e storia completa
/**
 * Per utilizzarla, basta creare una variabile con nome. Ogni volta che viene assegnata, aggiunge alle statistiche il proprio valore:
 * \code
 * auto_stats<double> E("Energia cinetica");
 * E=1;
 * E=3;
 * E=0.5;
 * E.print();  //Energia cinetica: 1.5 +- 1.08
 * \endcode
 */
template <typename T>
class auto_stats{
    private:
        ///Quanti valori diversi sono stati assegnati
        int counter;
        ///Somma, per la media
        double sum;
        ///Somma dei quadrati, per la varianza
        double sum2;
        ///Ultimo valore assegnato, per poter leggere
        T value;
        ///Nome della variabile
        std::string name;
        ///Storia completa della variabile
        std::vector<T> history;
    public:
        ///Operatore di assegnazione(generico), che implementa il calcolo delle statistiche ogni volta che un valore viene assegnato
        void operator=(T val){
            value=val;
            counter++;
            sum+=static_cast<double>(val);
            sum2+=static_cast<double>(val*val);
            history.push_back(val);
        }
        ///Lettura della variabile come un normale valore di tipo T (solo rvalue)
        operator T(){
            return value;
        }
        ///Media
        double mean(){
            return sum/counter;
        }
        ///Variazione standard con attenzione a troncamenti numerici per evitare radici di numeri negativi
        double std(){
            if (var()<1e-13)
                return 0;
            else
                return std::sqrt(var());
        }
        ///Varianza
        double var(){
            return (sum2/counter - sum/counter*sum/counter);
        }
        ///Stampa delle statistiche
        void print(){
            printf("%s\t %.3g \t %.3g\n",name.c_str(),mean(),std());
        }
        ///Quanti valori contiene
        int size(){
            return counter;
        }
        ///Restituisce storia read-only
        const std::vector<T> & read_history(){
            return history;
        }
        ///Costruttore con valore di default zero e battesimo della variabile
        auto_stats(std::string _name="") : counter(0), sum(0.0), sum2(0.0), name(_name) {}
        ///Distruttore che stampa automaticamente le statistiche accumulate
        ~auto_stats(){
            print();
        }
};

///Matrice con dimensioni definite a runtime e protezione accessi oltre i limiti
template <typename data_t>
class matrix {
    ///Buffer per i dati, acceduti poi tramite righe e colonne
    std::vector<data_t> buffer;
    ///Numero di righe
    int nrows;
    ///Numero di colonne
    int ncols;
    ///Dimensioni totali
    int total;
public:

    ///Costruisce una matrice con numero di righe e colonne fornite a runtime, riempita con \e fill
    matrix(int rows, int cols, data_t fill = 0) {
        total = rows*cols;
        nrows = rows;
        ncols = cols;
        if (nrows) {
#ifdef DEBUG
            fprintf(stderr, "Requested allocation space of %d kbytes\n", (((int) sizeof (data_t)) * total) >> 10);
#endif
            buffer.reserve(total);
            buffer.assign(total, fill);
        }
    }

    ///Accesso tramite riga e colonna. Ad es: matrice(5,3)=1. Se gli indici sono troppo grossi, stampa un errore
    inline data_t & operator()(int row, int col) {
        if (row >= nrows || col >= ncols) {
            std::stringstream msg;
            msg << "Requesting matrix element out of bounds at (" << row << "," << col << ")!";
            throw(std::out_of_range(msg.str()));
        }
        return buffer[ncols * row + col];
    }
};

///Matrice con numero di colonne definito a run-time, permette grandi ottimizzazioni
template <int cols>
class memory {
    ///Buffer per i dati, acceduti poi tramite righe e colonne
    std::vector<int> buffer;
    ///Numero di righe
    int nrows;
    ///Dimensione totale
    int total;
public:

    ///Crea una struttura con \e rows righe, riempita con valore \e fill
    memory(int rows = 0, int fill = 0) {
        total = rows*cols;
        nrows = rows;
        if (nrows) {
#ifdef DEBUG
            fprintf(stderr, "Requested allocation space of %d kbytes\n", (((int) sizeof (int)) * total) >> 10);
#endif
            buffer.reserve(total);
            buffer.assign(total, fill);
        }
    }
    ///Operatore per leggere e scrivere tramite righe e colonne, es: nn(1,4)=5
    inline int & operator()(int row, int col) {
        return buffer[cols * row + col];
    }

    ///Operatore per accedere al contenitore in modo sequenziale
    inline int & operator[](int row) {
        return buffer[row];
    }

    ///Sinonimo di []
    inline int & operator()(int row) {
        return buffer[row];
    }
};

#endif	/* MATRIX_H */

