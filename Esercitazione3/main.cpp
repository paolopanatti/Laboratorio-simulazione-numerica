#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

// Funzione per la stima dell'incertezza statistica
double error(vector<double>, vector<double>, int);

int main (int argc, char *argv[]){

    Random rnd;
    int seed[4];
    int p1, p2;

    // Lettura dei numeri primi da "Primes"
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    // Lettura del seed da "seed.in"
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

    // Parametri
    int M = 100000; // Numero totale di simulazioni
    int N = 100;    // Numero di blocchi
    int L = M/N;    // Numero di simulazioni per blocco

    double S0 = 100.;   // Prezzo iniziale
    double T = 1.;      // Tempo a scadenza
    double K = 100.;    // Strike price
    double r = 0.1;     // Tasso risk-free
    double sigma = 0.25;// Volatilità

    vector<double> av;
    vector<double> av2;
    vector<double> sum_prog(N);
    vector<double> sum2_prog(N);
    vector<double> err_prog(N);

    ofstream fileout;

    // ============================
    // Punto 1A: Call diretta
    // ============================
    fileout.open("data1.dat");
    fileout << M << " " << N << endl;

    for (int i = 0; i < N; i++){
        double C = 0;
        for (int j = 0; j < L; j++){
            double z = rnd.Gauss(0,1);
            // Simulazione diretta alla scadenza
            double S = S0*exp((r - pow(sigma,2)/2)*T + sigma*z*sqrt(T));
            C += exp(-r*T)*fmax(0,S-K);
        }
        av.push_back(C/L);            // media nel blocco
        av2.push_back(pow(av[i],2));  // quadrato della media
    }

    // Calcolo progressivo e errore
    for (int i = 0; i < N; i++){
        for (int j = 0; j < i+1; j++){
            sum_prog[i] += av[j];
            sum2_prog[i] += av2[j];
        }
        sum_prog[i] /= (i+1);  // media cumulativa
        sum2_prog[i] /= (i+1); // media quadratica cumulativa
        err_prog[i] = error(sum_prog, sum2_prog, i);
        fileout << sum_prog[i] << " " << err_prog[i] << endl;
    }
    fileout.close();

    // ============================
    // Punto 1B: Put diretta
    // ============================
    fileout.open("data2.dat");

    for (int i = 0; i < N; i++){
        double P = 0;
        for (int j = 0; j < L; j++){
            double z = rnd.Gauss(0,1);
            double S = S0*exp((r - pow(sigma,2)/2)*T + sigma*z*sqrt(T));
            P += exp(-r*T)*fmax(0,K-S);
        }
        av[i] = P/L;
        av2[i] = pow(av[i],2);
    }

    for (int i = 0; i < N; i++){
        sum_prog[i] = 0;
        sum2_prog[i] = 0;
        err_prog[i] = 0;
        for (int j = 0; j < i+1; j++){
            sum_prog[i] += av[j];
            sum2_prog[i] += av2[j];
        }
        sum_prog[i] /= (i+1);
        sum2_prog[i] /= (i+1);
        err_prog[i] = error(sum_prog, sum2_prog, i);
        fileout << sum_prog[i] << " " << err_prog[i] << endl;
    }
    fileout.close();

    // ============================
    // Punto 2A: Call tramite simulazione del cammino (100 step)
    // ============================
    fileout.open("data3.dat");

    for (int i = 0; i < N; i++){
        double C = 0;
        for (int j = 0; j < L; j++){
            double S1 = S0;
            double S2;
            for (int k = 0; k < 100; k++){ // evoluzione passo passo
                double z = rnd.Gauss(0,1);
                S2 = S1*exp((r - pow(sigma,2)/2)*T/100 + sigma*z*sqrt(T/100));
                S1 = S2;
            }
            C += exp(-r*T)*fmax(0,S2-K);
        }
        av[i] = C/L;
        av2[i] = pow(av[i],2);
    }

    for (int i = 0; i < N; i++){
        sum_prog[i] = 0;
        sum2_prog[i] = 0;
        err_prog[i] = 0;
        for (int j = 0; j < i+1; j++){
            sum_prog[i] += av[j];
            sum2_prog[i] += av2[j];
        }
        sum_prog[i] /= (i+1);
        sum2_prog[i] /= (i+1);
        err_prog[i] = error(sum_prog, sum2_prog, i);
        fileout << sum_prog[i] << " " << err_prog[i] << endl;
    }
    fileout.close();

    // ============================
    // Punto 2B: Put tramite simulazione del cammino (100 step)
    // ============================
    fileout.open("data4.dat");

    for (int i = 0; i < N; i++){
        double P = 0;
        for (int j = 0; j < L; j++){
            double S1 = S0;
            double S2;
            for (int k = 0; k < 100; k++){
                double z = rnd.Gauss(0,1);
                S2 = S1*exp((r - pow(sigma,2)/2)*T/100 + sigma*z*sqrt(T/100));
                S1 = S2;
            }
            P += exp(-r*T)*fmax(0,K-S2);
        }
        av[i] = P/L;
        av2[i] = pow(av[i],2);
    }

    for (int i = 0; i < N; i++){
        sum_prog[i] = 0;
        sum2_prog[i] = 0;
        err_prog[i] = 0;
        for (int j = 0; j < i+1; j++){
            sum_prog[i] += av[j];
            sum2_prog[i] += av2[j];
        }
        sum_prog[i] /= (i+1);
        sum2_prog[i] /= (i+1);
        err_prog[i] = error(sum_prog, sum2_prog, i);
        fileout << sum_prog[i] << " " << err_prog[i] << endl;
    }
    fileout.close();

    rnd.SaveSeed();
    return 0;
}

// Funzione errore statistico
double error(vector<double> AV, vector<double> AV2, int n){
    if (n==0){
        return 0;
    } else{
        return sqrt((AV2[n]-pow(AV[n],2))/n);
    }
}