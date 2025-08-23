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
   Random* _point_rnd;
   int _ncities = 110;
   Col<int> _path;
   double _length;
   double _pm = 0.3;

 public:
   Individuo& operator=(const Individuo& other);

   int getcities(){return _ncities;}
   void initialize(const mat* distanze, Random& rnd);
   bool checkbonds();
   int pbc(int pos);
   Col<int> getpath(){return _path;}
   void setpath(Col<int> path){_path = path;}
   int getstop(int pos){return _path(pos);}
   void setstop(int pos, int city){_path(pos) = city;}
   double getlength() const{return _length;}
   void calc_length(const mat* distanze);

   // Operatori di mutazione
   void swap(int i, int j);
   void swap();
   void shift();
   void swap_block();
   void invert_block();
   void mutate();
    
};

#endif // __Individuo__