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
   Random _rnd;
   mat _distanze;
   int _nindividuals;
   field <Individuo> _pop;
   const double _selexp = 2.4;
   const double _pcross = 0.8;
   int _ngen;

 public:
   void initialize_pop(int nindividuals, int ngen, mat* distanze, int rank);
   int getdim(){return _nindividuals;}
   Individuo get_percorso(int i);
   void sort_by_length();
   Individuo select();
   Col<int> sort_by_reference(Col<int> a, Col<int> ref);
   field<Individuo> crossover(Individuo a, Individuo b);
   void evolve(int size, int rank, MPI_Comm comm, int Nmigr);
   void migrate(int size, int rank, MPI_Comm comm, Random& rnd);

};

#endif