#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){

    Random rnd;
    int seed[4];
    int p1, p2;

    //Legge i numeri primi dal file "Primes" per inizializzare il generatore
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    //Legge il seed dal file "seed.in" e inizializza il generatore
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

    //Punto 2: generazione di medie di variabili casuali con diverse distribuzioni

    int N[4] = {1,2,10,100};  //Numero di variabili da mediare
    int M = 10000;            //Numero di esperimenti
    vector<double> s1(4);     //Medie per distribuzione uniforme
    vector<double> s2(4);     //Medie per distribuzione esponenziale
    vector<double> s3(4);     //Medie per distribuzione di Cauchy

    ofstream fileout1;
    ofstream fileout2;
    ofstream fileout3;
    fileout1.open("standard.dat");
    fileout2.open("exp.dat");
    fileout3.open("cauchy.dat");

    //Ciclo sugli M esperimenti
    for (int i = 0; i < M; i++){
        for (int j = 0; j < 4; j++){
            //Somma di N[j] variabili per ogni distribuzione
            for (int k = 0; k < N[j]; k++){
                s1[j] += rnd.Rannyu();     //Uniforme [0,1)
                s2[j] += rnd.Exp(1);       //Esponenziale con lambda=1
                s3[j] += rnd.Cauchy(0,1);  //Cauchy centrata in 0 con gamma=1
            }
            //Media delle variabili
            s1[j] /= N[j];
            s2[j] /= N[j];
            s3[j] /= N[j];
        }
        //Scrittura delle medie nei file
        fileout1 << s1[0] << " " << s1[1] << " " << s1[2] << " " << s1[3] << endl;
        fileout2 << s2[0] << " " << s2[1] << " " << s2[2] << " " << s2[3] << endl;
        fileout3 << s3[0] << " " << s3[1] << " " << s3[2] << " " << s3[3] << endl;
        //Reset dei vettori per la prossima iterazione
        for (int j = 0; j < 4; j++){
            s1[j] = 0;
            s2[j] = 0;
            s3[j] = 0;
        }
    }

    fileout1.close();
    fileout2.close();
    fileout3.close();

    rnd.SaveSeed();
    return 0;
}