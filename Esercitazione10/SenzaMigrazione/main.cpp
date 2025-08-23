#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "individuo.h"
#include "popolazione.h"
#include "mpi.h"

using namespace std;
using namespace arma;

void generate_cities(int choice, int n_cities, mat& distanze);
void load_cities(mat& distanze);

int main(int argc, char* argv[]){

    int size, rank;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int n_cities = 110;
    int n_individuals = 1000*4;
    int n_gen = 1000;
    int Nmigr = 10000;

    mat distanze(n_cities,n_cities);
    load_cities(distanze);

    Popolazione pop;
    pop.initialize_pop(n_individuals,n_gen,&distanze,rank);
    pop.evolve(size,rank,MPI_COMM_WORLD,Nmigr);

    MPI_Finalize();
    return 0;
}

void generate_cities(int choice, int n_cities, mat& distanze){
    Random rnd_cities;
    int p1, p2; // Read from ../Primes a pair of numbers to be used to initialize the RNG
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    int seed[4]; // Read the seed of the RNG
    ifstream Seed("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd_cities.SetRandom(seed,p1,p2);
    Seed.close();

    ofstream outc("cities.dat");
    outc << "#    CITY:             x:             y:" << endl;
    vec x(n_cities);
    vec y(n_cities);
    if(choice == 0){
        for(int i = 0; i < n_cities; i++){
            double theta = rnd_cities.Rannyu(0,2*M_PI);
            x[i] = cos(theta);
            y[i] = sin(theta);
            outc << setw(10) << i+1 << setw(15) << x[i] << setw(15) << y[i] << endl;
        }
    }
    else{
        for(int i = 0; i < n_cities; i++){
            x[i] = rnd_cities.Rannyu();
            y[i] = rnd_cities.Rannyu();
            outc << setw(10) << i+1 << setw(15) << x[i] << setw(15) << y[i] << endl;
        }
    }
    for(int i = 0; i < n_cities; i++){
        for(int j = 0; j < i+1; j++){
            distanze(i,j) = pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
            distanze(j,i) = distanze(i,j);
        }
    }
    outc.close();
}

void load_cities(mat& distanze){
    int n_cities = 110;
    vec x(n_cities);
    vec y(n_cities);

    ifstream filein("cap_prov_ita.dat");
    if(!filein){
        cerr << "Errore: file cap_prov_ita.dat non trovato" << endl;
        exit(1);
    }
    for(int i = 0; i < n_cities; i++){
        filein >> x[i] >> y[i];
    }
    filein.close();

    for(int i = 0; i < n_cities; i++){
        for(int j = 0; j < i+1; j++){
            distanze(i,j) = pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
            distanze(j,i) = distanze(i,j);
        }
    }
}