#ifndef __Popolazione__
#define __Popolazione__

#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>
#include "random.h"
#include "individuo.h"

using namespace std;
using namespace arma;

class Popolazione{

 private:
    Random _rnd;                // Generatore di numeri casuali
    mat _distanze;              // Matrice delle distanze tra le città
    int _nindividuals;          // Numero di individui nella popolazione
    field <Individuo> _pop;     // Vettore di individui
    const double _selexp = 2.4; // Esponente per la selezione degli individui
    const double _pcross = 0.8; // Probabilità di crossover
    int _ngen;                   // Numero di generazioni

 public:
    void initialize_pop(int nindividuals, int ngen, mat* distanze); // Inizializza la popolazione
    int getdim(){return _nindividuals;}                             // Restituisce il numero di individui
    Individuo get_individual(int i);                                 // Restituisce l'i-esimo individuo
    void sort_by_length();                                           // Ordina la popolazione in base alla lunghezza del percorso
    Individuo select();                                              // Seleziona un individuo secondo probabilità
    Col<int> sort_by_reference(Col<int> a, Col<int> ref);           // Ordina un vettore in base a un riferimento
    field<Individuo> crossover(Individuo a, Individuo b);           // Applica crossover tra due individui
    void evolve();                                                   // Evolve la popolazione per _ngen generazioni

};

#endif