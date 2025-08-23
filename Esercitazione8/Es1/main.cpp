#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <armadillo>
#include <algorithm>
#include "random.h"

using namespace std;

// Funzione d'onda di prova: somma di due gaussiane centrate in ±mu
double psi(double x, double mu, double sigma);

// Energia locale: (Hψ/ψ)(x)
double Hpsi(double x, double mu, double sigma);

// Algoritmo di Metropolis: decide se accettare nuova configurazione
bool metro(double x, double y, double mu, double sigma, Random& rnd);

// Funzione per calcolare l'errore statistico con blocking method
double error(double acc, double acc2, int blk);

// Funzione che regola automaticamente delta per ottenere ~50% di accettanza
double set_delta(Random& rnd, double mu, double sigma, int steps);

int main (int argc, char *argv[]){

    // Inizializzazione generatore di numeri casuali
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
 
    ifstream input("seed.in");
    string property;
    if (input.is_open()){
       while ( !input.eof() ){
          input >> property;
          if( property == "RANDOMSEED" ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             rnd.SetRandom(seed,p1,p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    // Parametri della simulazione
    double x = 0.; // configurazione iniziale
    double mu = 0.817861;   // parametro mu della wave function
    double sigma = 0.605328; // parametro sigma della wave function

    int n_steps = 10000; // numero di passi per blocco
    int n_blocks = 100;  // numero di blocchi
    int tune_steps = 1000; // passi per tuning del delta

    // Calibrazione di delta per ottenere acceptance ~0.5
    double delta = set_delta(rnd,mu,sigma,tune_steps);

    // File di output
    ofstream couta("acceptance.dat"); // accettanza per blocco
    couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;

    ofstream coute("energy.dat"); // energia stimata
    coute << "#     BLOCK:   ACTUAL_E:     E_AVE:       ERROR:" << endl;

    ofstream coutc("configurations.dat"); // configurazioni campionate

    int n_attemps;
    int n_accepted;

    double global_av = 0.;  // media cumulativa energia
    double global_av2 = 0.; // media cumulativa energia^2

    // Ciclo sui blocchi
    for(int i = 0; i < n_blocks; i++){
        double measure = 0.;
        n_attemps = 0;
        n_accepted = 0;
        // Ciclo sui passi all'interno del blocco
        for(int j = 0; j < n_steps; j++){
            n_attemps++;
            double y = x + rnd.Rannyu(-1.0,1.0) * delta; // proposta nuova configurazione
            if(metro(x,y,mu,sigma,rnd)){ // Metropolis
                x = y;
                n_accepted++;
            }
            measure += Hpsi(x,mu,sigma); // accumulo energia locale
            coutc << x << endl; // salvo configurazione
        }
        // Salvo accettanza per il blocco
        couta << setw(12) << i+1 << setw(12) << double(n_accepted)/double(n_attemps) << endl;
        // Aggiorno medie globali
        global_av += measure/double(n_steps);
        global_av2 += pow(measure/double(n_steps),2);
        // Stampo su file energia istantanea, cumulativa ed errore
        coute << setw(12) << i+1
              << setw(12) << measure/double(n_steps)
              << setw(12) << global_av/double(i+1)
              << setw(12) << error(global_av, global_av2, i+1) << endl;
    }

    // Chiudo i file
    couta.close();
    coute.close();
    coutc.close();

    return 0;
}

// Funzione d’onda trial: somma di due gaussiane centrate in ±mu
double psi(double x, double mu, double sigma){
    return exp(-pow(x-mu,2)/(2*pow(sigma,2)))+exp(-pow(x+mu,2)/(2*pow(sigma,2)));
}

// Energia locale Hψ/ψ: calcola derivata seconda analitica + potenziale
double Hpsi(double x, double mu, double sigma){
    return -0.5 * (-1./pow(sigma,2) + (pow(x-mu,2)/pow(sigma,4) * exp(-pow(x-mu,2)/(2*pow(sigma,2))) + pow(x+mu,2)/pow(sigma,4) * exp(-pow(x+mu,2)/(2*pow(sigma,2))))/psi(x,mu,sigma)) + pow(x,4) - 5./2.*pow(x,2);
}

// Algoritmo Metropolis: accetta con probabilità min(1, |ψ(y)|² / |ψ(x)|²)
bool metro(double x, double y, double mu, double sigma, Random& rnd){
    bool decision = false;
    double acceptance = min(1.0,pow(psi(y,mu,sigma),2)/pow(psi(x,mu,sigma),2));
    if(rnd.Rannyu() < acceptance){
        decision = true;
    }
    return decision;
}

// Calcolo errore statistico con metodo dei blocchi
double error(double acc, double acc2, int blk){
    if(blk <= 1) return 0.0;
    else return sqrt(fabs(acc2/double(blk)-pow(acc/double(blk),2))/double(blk-1));
}

// Regolazione automatica di delta per ottenere acceptance ~0.5
double set_delta(Random& rnd, double mu, double sigma, int steps){
    double delta = 0.5;
    double acceptance;
    double target_acceptance = 0.5; // valore target
    double tol = 0.1; // tolleranza
    double x = 0.; // punto di partenza
    int max_iter = 1000;
    int iter = 0;
    while(iter < max_iter){
        int accepted = 0;
        for(int i = 0; i < steps; i++){
            double y = x + rnd.Rannyu(-1.0,1.0) * delta;
            if(metro(x,y,mu,sigma,rnd)){
                x = y;
                accepted++;
            }
        }
        acceptance = double(accepted)/steps;
        if(fabs(acceptance-target_acceptance) < tol){
            break;
        }
        if(acceptance < target_acceptance)
            delta *= 0.9; // riduci se acceptance troppo bassa
        else
            delta *= 1.1; // aumenta se troppo alta

        iter++;
    }
    return delta;
}