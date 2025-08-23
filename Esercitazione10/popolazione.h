#ifndef __Popolazione__
#define __Popolazione__

#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>
#include "random.h"
#include "individuo.h"
#include "mpi.h"

using namespace std;
using namespace arma;

class Popolazione{

 private:
   Random _rnd;               // Generatore di numeri casuali
   mat _distanze;             // Matrice delle distanze tra le città
   int _nindividuals;         // Numero di individui nella popolazione
   field <Individuo> _pop;    // Popolazione di individui
   const double _selexp = 2.4; // Esponente per selezione probabilistica
   const double _pcross = 0.8; // Probabilità di crossover
   int _ngen;                 // Numero di generazioni

 public:
   // Inizializza la popolazione con individui casuali e RNG differente per ogni rank MPI
   void initialize_pop(int nindividuals, int ngen, mat* distanze, int rank);

   int getdim(){return _nindividuals;} // Restituisce la dimensione della popolazione
   Individuo get_individual(int i);    // Restituisce l'individuo i-esimo

   void sort_by_length();              // Ordina la popolazione in base alla lunghezza del percorso
   Individuo select();                 // Seleziona un individuo secondo probabilità ponderata
   Col<int> sort_by_reference(Col<int> a, Col<int> ref); // Ordina vettore secondo referenza

   field<Individuo> crossover(Individuo a, Individuo b); // Applica crossover tra due individui
   void evolve(int size, int rank, MPI_Comm comm, int Nmigr); // Evoluzione della popolazione con migrazione
   void migrate(int size, int rank, MPI_Comm comm, Random& rnd); // Migrazione di individui tra processi MPI

};

#endif