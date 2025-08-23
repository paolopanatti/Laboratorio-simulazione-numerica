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

    //Lettura dei numeri primi da "Primes"
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    //Lettura del seed da "seed.in"
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

    //Punto 1: random walk su reticolo cubico

    int M = 100000; //Numero totale di cammini
    int N = 100;    //Numero di blocchi
    int L = M/N;    //Numero di cammini in ciascun blocco

    double x;
    vector<int> pos(6);      //Conta passi nelle 6 direzioni ±x, ±y, ±z
    vector<double> dist(100);//Distanza media per ciascun passo

    vector<double> av(100);  //Medie dei blocchi
    vector<double> av2(100); //Quadrati delle medie

    ofstream fileout;
    fileout.open("data1.dat");

    for (int i = 0; i < N; i++){
        for (int j = 0; j < L; j++){
            for (int k = 0; k < 100; k++){
                x = rnd.Rannyu(); //Numero uniforme [0,1)
                //Decide la direzione del passo tra le 6 possibili
                for (int m = 0; m < 6; m++){
                    if (x >= m/6. and x < (m+1)/6.){
                        pos[m]++;
                        break;
                    }
                }
                //Distanza quadratica = somma dei quadrati degli spostamenti
                dist[k] += pow(pos[0]-pos[1],2)+pow(pos[2]-pos[3],2)+pow(pos[4]-pos[5],2);
            }
            //Reset conteggi direzioni per il cammino successivo
            for (int m = 0; m < 6; m++){
                pos[m] = 0;
            }
        }
        //Media e accumulo dei risultati
        for (int k = 0; k < 100; k++){
            dist[k] /= L;
            dist[k] = sqrt(dist[k]);
            av[k] += dist[k];
            av2[k] += pow(dist[k],2);
            dist[k] = 0;
        }
    }
    //Scrittura dei risultati con errori statistici
    for (int k = 0; k < 100; k++){
        av[k] /= N;
        av2[k] /= N;
        fileout << av[k] << " " << sqrt((av2[k]-pow(av[k],2))/(N-1)) << endl;
        av[k] = 0;
        av2[k] = 0;
    }

    fileout.close();

    //Punto 2: random walk nel continuo (passi isotropi in 3D)

    double a = 1.;  //Lunghezza del passo
    double theta;
    double phi;
    double y = 0, z = 0;
    x = 0;

    fileout.open("data2.dat");

    for (int i = 0; i < N; i++){
        for (int j = 0; j < L; j++){
            for (int k = 0; k < 100; k++){
                //Genera direzione isotropa
                theta = rnd.Theta();
                phi = rnd.Rannyu(0,2*M_PI);
                //Aggiorna posizione
                x += a*sin(theta)*cos(phi);
                y += a*sin(theta)*sin(phi);
                z += a*cos(theta);
                //Accumula distanza quadratica
                dist[k] += pow(x,2)+pow(y,2)+pow(z,2);
            }
            //Reset posizione per nuovo cammino
            x = 0;
            y = 0;
            z = 0;
        }
        //Media e accumulo dei risultati
        for (int k = 0; k < 100; k++){
            dist[k] /= L;
            dist[k] = sqrt(dist[k]);
            av[k] += dist[k];
            av2[k] += pow(dist[k],2);
            dist[k] = 0;
        }
    }
    //Scrittura dei risultati con errori statistici
    for (int k = 0; k < 100; k++){
        av[k] /= N;
        av2[k] /= N;
        fileout << av[k] << " " << sqrt((av2[k]-pow(av[k],2))/(N-1)) << endl;
        av[k] = 0;
        av2[k] = 0;
    }

    fileout.close();

    rnd.SaveSeed();
    return 0;
}