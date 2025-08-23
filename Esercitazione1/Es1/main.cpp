#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

//Funzione per la stima dell'incertezza statistica (data la media e la media dei quadrati)
double error(vector<double>, vector<double>, int);
 
int main (int argc, char *argv[]){

    Random rnd;
    int seed[4];
    int p1, p2;

    //Legge due numeri primi dal file "Primes" per inizializzare il generatore
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    //Legge il seed da "seed.in" e inizializza il generatore
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

    int M = 100000; //Numero totale di lanci
    int N = 100;    //Numero di blocchi
    int L = M/N;    //Numero di lanci in ciascun blocco

    //Punto 1: stima del valor medio e relativa incertezza

    vector<double> av;
    vector<double> av2;
    vector<double> sum_prog(N);
    vector<double> sum2_prog(N);
    vector<double> err_prog(N);

    ofstream fileout;
    fileout.open("data_1.dat");
    fileout << M << " " << N << endl;

    //Calcolo delle medie per ciascun blocco
    for (int i = 0; i < N; i++){
        double sum1 = 0;
        for (int j = 0; j < L; j++){
            sum1 += rnd.Rannyu(); //Genera numero uniforme [0,1)
        }
        av.push_back(sum1/L); // A_i 
        av2.push_back(pow(av[i],2)); // (A_i)^2
    }

    //Calcolo delle medie progressive e dell’incertezza statistica
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

    //Punto 2: stima del valor medio di (x-0.5)^2 con incertezza

    fileout.open("data_2.dat");

    for (int i = 0; i < N; i++){
        double sum2 = 0;
        for (int j = 0; j < L; j++){
            sum2 += pow((rnd.Rannyu()-0.5),2);
        }
        av[i] = sum2/L;
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

    //Punto 3: test del Chi-quadro per verificare la distribuzione uniforme

    M = 100;        //Numero di intervalli
    N = 10000;      //Numero di lanci per ogni test
    int NTest = 100;//Numero di test indipendenti

    vector<double> chi2(NTest);
    vector<int> obs(M);
    double r;

    fileout.open("chi2.dat");

    //Ciclo sui test
    for (int i = 0; i < NTest; i++){
        //Distribuzione degli N numeri generati in M intervalli
        for (int j = 0; j < N; j++){
            r = rnd.Rannyu();
            for (int k = 0; k < M; k++){
                if (r >= double(k)/M and r < double(k+1)/M){
                    obs[k]++;
                    break;
                }
            }
        }
        //Calcolo del chi-quadro
        for (int k = 0; k < M; k++){
            chi2[i] += pow((obs[k]-N/M),2)/(N/M);
            obs[k] = 0; //reset conteggi
        }
        fileout << chi2[i] << endl;
    }

    fileout.close();

    rnd.SaveSeed();
    return 0;
}

//Funzione che calcola l’errore statistico sulla media progressiva
double error(vector<double> AV, vector<double> AV2, int n){
    if (n==0){
        return 0;
    } else{
        return sqrt((AV2[n]-pow(AV[n],2))/n);
    }
}