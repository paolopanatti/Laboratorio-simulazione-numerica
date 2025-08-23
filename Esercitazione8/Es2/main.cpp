#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <armadillo>
#include <algorithm>
#include "random.h"

using namespace std;

// Funzione d’onda di prova ψ(x; μ,σ): somma di due gaussiane centrate in ±μ
double psi(double x, double mu, double sigma);

// Energia locale Hψ/ψ calcolata analiticamente
double Hpsi(double x, double mu, double sigma);

// Algoritmo Metropolis per campionare |ψ(x)|²
bool metro_H(double x, double y, double mu, double sigma, Random& rnd);

// Calcolo dell’errore statistico (blocking method)
double error(double acc, double acc2, int blk);

// Regolazione automatica del passo Metropolis (delta)
double set_delta(Random& rnd, double mu, double sigma, int steps);

// Calcolo di energia e errore per dati μ e σ
vector<double> compute_H (Random& rnd, double mu, double sigma);

// Metropolis per l’aggiornamento dei parametri (μ,σ) nel simulated annealing
bool metro_par(double H_old, double H_new, double beta, Random& rnd);

int main(int argc, char *argv[]){

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

    // Parametri iniziali per μ e σ
    double mu = 1.;
    double sigma = 1.;
    // Calcolo energia iniziale con questi parametri
    vector<double> H = compute_H(rnd,mu,sigma);

    // File di output parametri + energia
    ofstream coutp("parameters.dat");
    coutp << "#         T:            mu:         sigma:             H:           err:" << endl;

    // Ciclo di simulated annealing sulla temperatura
    for(double T = 2.; T >= 0.01; T *= 0.99){
        double beta = 1./T; // β = 1/T
        double delta_mu = 0.5 * T;    // ampiezza proposte per μ
        double delta_sigma = 0.5 * T; // ampiezza proposte per σ
        // Per ogni temperatura faccio 100 proposte di nuovi parametri
        for(int n = 0; n < 100; n++){
            double mu_new = fabs(mu + rnd.Rannyu(-1.0,1.0) * delta_mu);
            double sigma_new;
            // σ deve rimanere positivo
            do{
                sigma_new = fabs(sigma + rnd.Rannyu(-1.0,1.0) * delta_sigma);
            } while(sigma_new <= 0.1);
            // Calcolo energia con i nuovi parametri
            vector<double> H_new = compute_H(rnd,mu_new,sigma_new);
            // Metropolis sui parametri (accept/reject)
            if(metro_par(H[0],H_new[0],beta,rnd)){
                mu = mu_new;
                sigma = sigma_new;
                H = H_new;
            }
        }
        // Stampo temperatura, parametri e energia media con errore
        coutp << setw(12) << T << setw(15) << mu << setw(15) << sigma << setw(15) << H[0] << setw(15) << H[1] << endl;
    }

    return 0;
}

vector<double> compute_H (Random& rnd, double mu, double sigma){
    
    // Parametri della simulazione Monte Carlo
    double x = 0.; // configurazione iniziale (posizione della particella)

    int n_steps = 1000;  // passi per blocco
    int n_blocks = 100;  // numero di blocchi
    int tune_steps = 1000; // passi per tarare delta

    // Regolo delta per avere acceptance ~0.5
    double delta = set_delta(rnd,mu,sigma,tune_steps);

    double global_av = 0.;
    double global_av2 = 0.;

    vector<double> H; // conterrà energia media ed errore

    // Ciclo a blocchi
    for(int i = 0; i < n_blocks; i++){
        double measure = 0.;
        // Campionamento dentro il blocco
        for(int j = 0; j < n_steps; j++){
            double y = x + rnd.Rannyu(-1.0,1.0) * delta; // nuova configurazione
            if(metro_H(x,y,mu,sigma,rnd)){ // accettazione Metropolis
                x = y;
            }
            measure += Hpsi(x,mu,sigma); // accumulo energia locale
        }
        // Aggiorno medie globali
        global_av += measure/double(n_steps);
        global_av2 += pow(measure/double(n_steps),2);
    }

    // Energia media e errore stimato
    H.push_back(global_av/double(n_blocks));
    H.push_back(error(global_av, global_av2, n_blocks));

    return H;
}

// Funzione d’onda trial ψ(x)
double psi(double x, double mu, double sigma){
    return exp(-pow(x-mu,2)/(2*pow(sigma,2)))+exp(-pow(x+mu,2)/(2*pow(sigma,2)));
}

// Energia locale Hψ/ψ
double Hpsi(double x, double mu, double sigma){
    return -0.5 * (-1./pow(sigma,2) + (pow(x-mu,2)/pow(sigma,4) * exp(-pow(x-mu,2)/(2*pow(sigma,2))) + pow(x+mu,2)/pow(sigma,4) * exp(-pow(x+mu,2)/(2*pow(sigma,2))))/psi(x,mu,sigma)) + pow(x,4) - 5./2.*pow(x,2);
}

// Algoritmo Metropolis per la configurazione x → y
bool metro_H(double x, double y, double mu, double sigma, Random& rnd){
    bool decision = false;
    double acceptance = min(1.0,pow(psi(y,mu,sigma),2)/pow(psi(x,mu,sigma),2));
    if(rnd.Rannyu() < acceptance){
        decision = true;
    }
    return decision;
}

// Calcolo errore statistico (blocking method)
double error(double acc, double acc2, int blk){
    if(blk <= 1) return 0.0;
    else return sqrt(fabs(acc2/double(blk)-pow(acc/double(blk),2))/double(blk-1));
}

// Regolazione automatica di delta per acceptance ~0.5
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
            if(metro_H(x,y,mu,sigma,rnd)){
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

// Metropolis per l’aggiornamento dei parametri (μ,σ) nel simulated annealing
bool metro_par(double H_old, double H_new, double beta, Random& rnd){
    bool decision = false;
    double acceptance = exp(-beta * (H_new - H_old)); // criterio Boltzmann per nuovi parametri
    if(rnd.Rannyu() < acceptance){
        decision = true; // accetta i nuovi parametri
    }
    return decision;
}