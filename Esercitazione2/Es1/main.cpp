#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

//Funzione per la stima dell'incertezza statistica
double error(vector<double>, vector<double>, int);
//Funzioni da integrare
double funzione1(double);
double funzione2(double);

int main (int argc, char *argv[]){

    Random rnd;
    int seed[4];
    int p1, p2;

    //Lettura dei numeri primi dal file "Primes"
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    //Lettura del seed dal file "seed.in"
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

    //Parametri simulazione
    int M = 100000; //Numero totale di lanci
    int N = 100;    //Numero di blocchi
    int L = M/N;    //Numero di lanci per blocco

    double x;

    vector<double> av;          //Valori medi dei blocchi
    vector<double> av2;         //Quadrati dei valori medi
    vector<double> sum_prog(N); //Medie progressive
    vector<double> sum2_prog(N);//Quadrati medie progressive
    vector<double> err_prog(N); //Errori statistici

    ofstream fileout;
    fileout.open("data1.dat");
    fileout << M << " " << N << endl;

    //Punto 1: integrazione con metodo Monte Carlo (uniforme)
    for (int i = 0; i < N; i++){
        double f = 0;
        for (int j = 0; j < L; j++){
            x = rnd.Rannyu();          //Genera x in [0,1)
            f += funzione1(x);         //Valuta la funzione integranda
        }
        av.push_back(f/L);             //Stima del blocco
        av2.push_back(pow(av[i],2));   //Quadrato della stima
    }

    //Medie progressive e incertezze
    for (int i = 0; i < N; i++){
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

    //Punto 2: integrazione con Importance Sampling
    fileout.open("data2.dat");

    for (int i = 0; i < N; i++){
        double f = 0;
        for (int j = 0; j < L; j++){
            x = rnd.Distr();           //Genera x secondo pdf scelta
            f += funzione2(x);         //Valuta funzione pesata
        }
        av[i] = f/L;                   //Stima del blocco
        av2[i] = pow(av[i],2);         //Quadrato della stima
    }

    //Medie progressive e incertezze
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

//Funzione che calcola l'errore statistico
double error(vector<double> AV, vector<double> AV2, int n){
    if (n==0){
        return 0;
    } else{
        return sqrt((AV2[n]-pow(AV[n],2))/n);
    }
}

//Funzione integranda con campionamento uniforme
double funzione1(double x){
    return M_PI/2*cos((M_PI*x)/2);
}

//Funzione integranda con Importance Sampling
double funzione2(double x){
    return M_PI/4*cos((M_PI*x)/2)/(1-x);
}