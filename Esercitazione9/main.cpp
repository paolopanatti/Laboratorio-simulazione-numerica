#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "individuo.h"
#include "popolazione.h"

using namespace std;
using namespace arma;

// Funzione per generare le coordinate delle città e calcolare la matrice delle distanze
void generate_cities(int choice, int n_cities, mat& distanze);

int main(){

    int choice;

    cout << "Inserisci 0 per città disposte su circonferenza, 1 per città disposte all'interno di un quadrato" << endl;
    cin >> choice;
    if(choice != 0 && choice != 1){
        cout << "Scelta non valida! Per favore inserire 0 oppure 1" << endl;
        return -1;
    }

    int n_cities = 34;         // Numero di città
    int n_individuals = 200;   // Numero di individui nella popolazione
    int n_gen = 500;           // Numero di generazioni

    mat distanze(n_cities,n_cities);
    generate_cities(choice,n_cities,distanze); // Genera città e calcola distanze

    Popolazione pop;
    pop.initialize_pop(n_individuals,n_gen,&distanze); // Inizializza popolazione
    pop.evolve();                                     // Evolve la popolazione

    return 0;
}

void generate_cities(int choice, int n_cities, mat& distanze){
    Random rnd_cities;
    int p1, p2; // Numeri primi per inizializzare RNG
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    int seed[4]; // Seed del RNG
    ifstream Seed("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd_cities.SetRandom(seed,p1,p2);

    ofstream outc("cities.dat");
    outc << "#    CITY:             x:             y:" << endl;
    vec x(n_cities);
    vec y(n_cities);

    if(choice == 0){ // città disposte su circonferenza
        for(int i = 0; i < n_cities; i++){
            double theta = rnd_cities.Rannyu(0,2*M_PI);
            x[i] = cos(theta);
            y[i] = sin(theta);
            outc << setw(10) << i+1 << setw(15) << x[i] << setw(15) << y[i] << endl;
        }
    }
    else{ // città disposte all'interno di un quadrato
        for(int i = 0; i < n_cities; i++){
            x[i] = rnd_cities.Rannyu();
            y[i] = rnd_cities.Rannyu();
            outc << setw(10) << i+1 << setw(15) << x[i] << setw(15) << y[i] << endl;
        }
    }

    // Calcolo della matrice delle distanze euclidee
    for(int i = 0; i < n_cities; i++){
        for(int j = 0; j < i+1; j++){
            distanze(i,j) = pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
            distanze(j,i) = distanze(i,j);
        }
    }
    outc.close();
}