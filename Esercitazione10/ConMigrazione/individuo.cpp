#include "individuo.h"

Individuo& Individuo:: operator=(const Individuo& other){
    if(this != &other){
        this->_point_rnd = other._point_rnd;
        this->_ncities = other._ncities;
        this->_path = other._path;
        this->_length = other._length;
    }
    return *this;
}

void Individuo:: initialize(const mat* distanze, Random& rnd){
    _point_rnd = &rnd;
    _path.resize(_ncities);
    for(int i = 0; i < _ncities; i++){
        _path(i) = i+1;
    }
    int n_swaps = _ncities;
    for(int n = 0; n < n_swaps; n++){
        int i = int(rnd.Rannyu(1,_ncities));
        int j = int(rnd.Rannyu(1,_ncities));
        swap(i,j);
    }
    calc_length(distanze);
}

bool Individuo:: checkbonds(){
    if(_path(0) != 1) return false;
    for(int i = 0; i < _ncities; i++){
        for(int j = i+1; j < _ncities; j++){
            if(_path(i) == _path(j)){
                return false;
            }
        }
    }
    return true;
}

int Individuo:: pbc(int pos){
    if(pos >= _ncities) pos = pos - _ncities + 1;
    return pos;
}

void Individuo:: calc_length(const mat* distanze){
    double length = 0;
    for(int i = 0; i < _ncities - 1; i++){
        length += distanze->at(_path(i)-1,_path(i+1)-1);
    }
    length += distanze->at(_path(_ncities-1)-1,_path(0)-1);
    _length = length;
}

void Individuo:: swap(int i, int j){
    int temp = _path(i);
    _path(i) = _path(j);
    _path(j) = temp;
}

void Individuo:: swap(){
    int i = int(_point_rnd->Rannyu(1,_ncities));
    int j = int(_point_rnd->Rannyu(1,_ncities));
    swap(i,j);
    if(!checkbonds()){
        cerr << "Error: path not valid after swap" << endl;
        exit(0);
    }
}

void Individuo:: shift(){
    int m = int(_point_rnd->Rannyu(1,_ncities-1));
    int start = int(_point_rnd->Rannyu(1,_ncities-m+1));
    int n = int(_point_rnd->Rannyu(1,_ncities-1));
    for(int i = 0; i < m; i++){
        swap(start+i,pbc(i+n));
    }
    if(!checkbonds()){
        cerr << "Error: path not valid after shift" << endl;
        exit(0);
    }
}

void Individuo:: swap_block(){
    int m = int(_point_rnd->Rannyu(1,_ncities/2));
    int start1 = int(_point_rnd->Rannyu(1,_ncities-2*m+1));
    int start2 = int(_point_rnd->Rannyu(start1+m,_ncities-m+1));
    for(int i = 0; i < m; i++){
        swap(start1+i,start2+i);
    }
    if(!checkbonds()){
        cerr << "Error: path not valid after swap_block" << endl;
        exit(0);
    }
}

void Individuo:: invert_block(){
    int m = int(_point_rnd->Rannyu(1,_ncities));
    int start = int(_point_rnd->Rannyu(1,_ncities-m+1));
    for(int i = 0; i < m/2; i++){
        swap(start+i,start+m-1-i);
    }
    if(!checkbonds()){
        cerr << "Error: path not valid after invert_block" << endl;
        exit(0);
    }
}

void Individuo:: mutate(){
    double p = _point_rnd->Rannyu();
    if(p < _pm){
        double sel = _point_rnd->Rannyu();
        if(sel < 0.25){
            swap();
        } else if(sel < 0.5){
            shift();
        } else if(sel < 0.75){
            swap_block();
        } else{
            invert_block();
        }
    }
}