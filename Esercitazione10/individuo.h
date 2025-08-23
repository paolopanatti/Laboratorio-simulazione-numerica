#ifndef __Individuo__
#define __Individuo__

#include <iostream>
#include <cmath>
#include <armadillo>
#include <math.h>
#include "random.h"

using namespace std;
using namespace arma;

class Individuo{

 private:
    Random* _point_rnd;       // Puntatore al generatore di numeri casuali
    int _ncities = 110;        // Numero totale di città
    Col<int> _path;           // Percorso dell'individuo
    double _length;           // Lunghezza totale del percorso
    double _pm = 0.3;         // Probabilità di mutazione

 public:
    Individuo& operator=(const Individuo& other); // Operatore di assegnazione

    // Getters e setters
    int getcities(){return _ncities;}
    Col<int> getpath(){return _path;}
    void setpath(Col<int> path){_path = path;}
    int getstop(int pos){return _path(pos);}
    void setstop(int pos, int city){_path(pos) = city;}
    double getlength() const{return _length;}

    // Inizializza percorso casuale
    void initialize(const mat* distanze, Random& rnd);

    // Controlla validità percorso
    bool checkbonds();

    // Boundary condition cicliche
    int pbc(int pos);

    // Calcola lunghezza percorso
    void calc_length(const mat* distanze);

    // Operatori di mutazione
    void swap(int i, int j);   // Scambia due città specifiche
    void swap();               // Scambia due città casuali
    void shift();              // Sposta un blocco di città
    void swap_block();         // Scambia due blocchi di città
    void invert_block();       // Inverte un blocco di città
    void mutate();             // Applica mutazione casuale con probabilità _pm
    
};

#endif // __Individuo__