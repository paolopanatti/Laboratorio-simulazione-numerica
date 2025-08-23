#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

//Funzione per la stima dell'incertezza statistica
double error(vector<double>, vector<double>, int);

int main (int argc, char *argv[]){

    Random rnd;
    int seed[4];
    int p1, p2;

    //Lettura dei numeri primi da "Primes" per inizializzare il generatore
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

    int M = 100000; //Numero totale di lanci
    int N = 100;    //Numero di blocchi
    int L = M/N;    //Numero di lanci in ciascun blocco

    int d = 3;      //Distanza tra le righe
    int l = 1;      //Lunghezza dell'ago
    int NHit = 0;   //Numero di aghi che intersecano una linea
    double x;       //Centro dell'ago
    double xp, yp, sin; //Coordinate per generare l'orientazione casuale

    vector<double> av;          //Valori medi blocco
    vector<double> av2;         //Quadrati dei valori medi blocco
    vector<double> sum_prog(N); //Medie progressive
    vector<double> sum2_prog(N);//Quadrati delle medie progressive
    vector<double> err_prog(N); //Errori statistici

    ofstream fileout;
    fileout.open("data.dat");
    fileout << M << " " << N << endl;

    //Simulazione dell’esperimento di Buffon
    for (int i = 0; i < N; i++){
        for (int j = 0; j < L; j++){
            //Posizione casuale del centro dell'ago
            x = rnd.Rannyu(0,5*d);

            //Generazione casuale dell'orientazione (uniforme su una semicirconferenza)
            do{
                xp = rnd.Rannyu(-1.,1.);
                yp = rnd.Rannyu(-1.,1.);
            } while(pow(xp,2)+pow(yp,2) >= 1);

            sin = fabs(yp/sqrt(pow(xp,2)+pow(yp,2))); //Seno dell’angolo

            //Controllo se l’ago interseca una linea
            for (int m = 0; m < 6; m++){
                if (m*d >= (x-l/2.*sin) and m*d <= (x+l/2.*sin)){
                    NHit++;
                    break;
                }
            }
        }
        //Stima di π per il blocco
        av.push_back((2.*l*L)/(NHit*d)); 
        av2.push_back(pow(av[i],2));     
        NHit = 0; //Reset contatore per il blocco successivo
    }

    //Calcolo delle medie progressive e degli errori
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

    rnd.SaveSeed();
    return 0;
}

//Funzione che calcola l’errore statistico
double error(vector<double> AV, vector<double> AV2, int n){
    if (n==0){
        return 0;
    } else{
        return sqrt((AV2[n]-pow(AV[n],2))/n);
    }
}