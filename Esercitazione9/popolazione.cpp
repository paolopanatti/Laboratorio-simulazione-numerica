#include "popolazione.h"
#include <algorithm>
#include <fstream>
#include <iomanip>

// Inizializza la popolazione con nindividuals e numero di generazioni ngen
void Popolazione:: initialize_pop(int nindividuals, int ngen, mat* distanze){
    int p1, p2; // Read from ../Primes a pair of numbers to be used to initialize the RNG
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    int seed[4]; // Read the seed of the RNG
    ifstream Seed("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    _rnd.SetRandom(seed,p1,p2);

    _nindividuals = nindividuals;
    _ngen = ngen;
    _distanze = *distanze;

    _pop.set_size(_nindividuals);
    for(int i = 0; i < _nindividuals; i++){
        _pop(i).initialize(&_distanze, _rnd); // Inizializza ogni individuo
    }
}

// Restituisce l'individuo i-esimo
Individuo Popolazione:: get_individual(int i){
    return _pop(i);
}

// Ordina la popolazione in base alla lunghezza del percorso
void Popolazione:: sort_by_length(){
    vector<Individuo> temp_vec(_nindividuals);
    for(int i = 0; i < _nindividuals; i++){
        temp_vec[i] = _pop(i);
    }
    std::sort(temp_vec.begin(), temp_vec.end(), [&](const Individuo& a, Individuo& b){
        return a.getlength() < b.getlength();
    });
    for(int i = 0; i < _nindividuals; i++){
        _pop(i) = temp_vec[i];
    }
}

// Seleziona un individuo con probabilità ponderata
Individuo Popolazione:: select(){ 
    double r = _rnd.Rannyu(); 
    int j = int(_nindividuals * pow(r, _selexp));
    return _pop(j);
}

// Ordina un vettore a in base all'ordine dei valori in ref
Col<int> Popolazione:: sort_by_reference(Col<int> a, Col<int> ref){
    Col<int> ref_pos(ref.size()+1);
    for(int i = 0; i < ref.size(); i++){
        ref_pos[ref[i]] = i;
    }

    vector<int> sorted(a.begin(),a.end());
    sort(sorted.begin(), sorted.end(), [&ref_pos](int j, int k){
        return ref_pos[j] < ref_pos[k];
    });

    return Col<int>(sorted);
}

// Crossover tra due individui a e b con probabilità _pcross
field<Individuo> Popolazione:: crossover(Individuo a, Individuo b){
    if(_rnd.Rannyu() < _pcross){
        field<Individuo> figli(2);
        figli[0].initialize(&_distanze,_rnd);
        figli[1].initialize(&_distanze,_rnd);

        Col<int> genitore1 = a.getpath();
        Col<int> genitore2 = b.getpath();

        int pos = int(_rnd.Rannyu(0,a.getcities()));
        Col<int> sub1 = genitore1.subvec(0,pos);
        Col<int> sub2 = genitore2.subvec(0,pos);

        Col<int> figlio1_temp = sub1;
        Col<int> figlio2_temp = sub2;

        if(pos < a.getcities()-1){
            Col<int> sub3 = genitore1.subvec(pos+1,a.getcities()-1);
            Col<int> sub4 = genitore2.subvec(pos+1,b.getcities()-1);
            Col<int> sub3_ord = sort_by_reference(sub3,genitore2);
            Col<int> sub4_ord = sort_by_reference(sub4,genitore1);
            figlio1_temp = join_vert(figlio1_temp,sub3_ord);
            figlio2_temp = join_vert(figlio2_temp,sub4_ord);
        }

        figli[0].setpath(figlio1_temp);
        figli[1].setpath(figlio2_temp);

        // Controlla validità percorsi
        if(!figli[0].checkbonds() or !figli[1].checkbonds()){
            cerr << "Error: path not valid after crossover" << endl;
            exit(0);
        }

        return figli;
    }
    else{
        return {a,b};
    }
}

// Evoluzione della popolazione per _ngen generazioni
void Popolazione:: evolve(){
    ofstream coutf("results.dat");
    coutf << "#     GEN:            L2:          <L2>:" << endl;
    sort_by_length(); // Ordina inizialmente
    for(int i = 0; i < _ngen; i++){
        field<Individuo> figli(_nindividuals);
        for(int j = 0; j < _nindividuals/2; j++){
            Individuo a = select(); // Seleziona genitori
            Individuo b = select();
            field<Individuo> appo = crossover(a,b); // Applica crossover
            figli[j].initialize(&_distanze,_rnd);
            figli[j].setpath(appo[0].getpath());
            figli[j+_nindividuals/2].initialize(&_distanze,_rnd);
            figli[j+_nindividuals/2].setpath(appo[1].getpath());
        }
        for(int j = 0; j < _nindividuals; j++){
            figli[j].mutate(); // Mutazione casuale
            figli[j].calc_length(&_distanze); // Ricalcola lunghezza percorso
        }
        _pop = figli;
        sort_by_length(); // Riordina popolazione
        double sum = 0;
        for (int j = 0; j < _nindividuals/2; j++){
            sum += _pop[j].getlength();
        }
        coutf << setw(10) << i+1 << setw(15) << _pop[0].getlength() << setw(15) << sum/(_nindividuals/2) << endl;
    }
    coutf.close();
    
    // Salva il percorso migliore
    ofstream out("best_path.dat");
    out << "# BEST PATH:" << endl;
    for(int i = 0; i < _pop[0].getcities(); i++){
        out << _pop[0].getpath()[i] << endl;
    }
    out.close();
}