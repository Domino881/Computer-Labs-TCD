/*
* Author: Dominik Kuczynski
*
* This code simulates an Ising Model using the Markov Chain Monte
* Carlo method. It is my rewrite of the original code, written
* during the TPSA workshop, from Python to C++ to improve performance.
*
* It outputs the result to the file "example.txt" to be plotted by
* "plotter.py".
*/

#include<cstdio>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<chrono>
using namespace std::chrono;
using namespace std;
const int xdim = 10;
const int ydim = 10;


class Lattice {
    public:
        short matrix[xdim][ydim];

        Lattice(short matrix[xdim][ydim]){
            for(int i=0;i<xdim;i++){
                for(int j=0;j<ydim;j++){
                    this->matrix[i][j] = matrix[i][j];
                }
            }
        }

        int sum_nn(int ix, int iy){
            int r = 0;
            r += this->matrix[(ix-1+xdim)%xdim][iy];
            r += this->matrix[(ix+1)%xdim][iy];
            r += this->matrix[ix][(iy-1+xdim)%ydim];
            r += this->matrix[ix][(iy+1)%ydim];
            return r;
        }

        int H(){
            int E = 0;
            for(int i=0;i<xdim;i++){
                for(int j=0;j<ydim;j++){
                    //printf("%d ", this->sum_nn(i,j));
                    E += this->matrix[i][j] * this->sum_nn(i,j);
                }
                //printf("\n");
            }
            return -E/2;
        }

        int M(){
            int M = 0;
            for (int i = 0; i < xdim; i++) {
                for (int j = 0; j < ydim; j++) {
                    M += this->matrix[i][j];
                }
            }
            return M;
            
        }

        void print(){
            for (int i = 0; i < xdim; i++) {
                for (int j = 0; j < ydim; j++) {
                    printf("%+4d ", this->matrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
            
        }
};

const int N = 2000;
const int burnin = 100;
const int ntemps = 50;

int energies[N+3];
int magnets[N+3];
double avg_energies[ntemps];
double avg_magnets[ntemps];
double temps[ntemps];

int main(){

    srand(time(NULL));

    short m[xdim][ydim];
    for (int  i = 0; i < xdim; i++) {
        for (int j = 0; j < ydim; j++) {
            m[i][j] = rand()>RAND_MAX/2? -1 : 1;
        }
        
    }

    Lattice l(m);

    printf("%d", l.H());
    l.print();

    int H=0; int M=0;

    for (int it=0;it<ntemps;it++){
        double t = 0.1 + (7-0.8)/double(ntemps-1) * it;
        double beta = 1.0/t;

        // auto start = high_resolution_clock::now();
        for (int n = 0; n < N; n++) {

            for (int ix=0;ix<xdim;ix++){
                for (int iy=0;iy<ydim;iy++){
                    int deltaH = 2 * l.matrix[ix][iy] * l.sum_nn(ix,iy);

                    if(deltaH <= 0){
                        H += deltaH;
                        l.matrix[ix][iy] = -l.matrix[ix][iy];
                        M += 2*l.matrix[ix][iy];
                    }
                    else{
                        double u = double(rand())/RAND_MAX;
                        if(u <= exp(-beta * deltaH)){
                            H += deltaH;
                            l.matrix[ix][iy] = -l.matrix[ix][iy];
                            M += 2 * l.matrix[ix][iy];
                        }
                    }
                }
            }
            energies[n]=H;
            magnets[n]=M;

            // l.print();
            // printf("%d\n", H);
        }
        // auto duration = duration_cast<microseconds>(high_resolution_clock::now()-start);
        // printf("time per sweep: %.0fus\n", (double)duration.count()/(N+burnin));
        double avg_E = 0.0;
        double avg_M = 0.0;
        for(int i=burnin;i<N;i++){
            avg_E += energies[i];
            avg_M += magnets[i];
        }
        avg_energies[it] = avg_E/double(N-burnin);
        avg_magnets[it] = avg_M/double(N-burnin);
        // printf("%.2f ", avg_energies[it]);
        temps[it] = t;
    }

    FILE *fp;

    fp = fopen("example.txt", "w");

    for(int j=0;j<ntemps;j++){
        fprintf(fp, "%f ", temps[j]);
    }

    fprintf(fp, "\n");

    for (int i = 0; i < ntemps; i++) {
        fprintf(fp, "%f ", avg_energies[i]);
    }

    fprintf(fp, "\n");

    for (int i = 0; i < ntemps; i++) {
        fprintf(fp, "%f ", abs(avg_magnets[i]));
    }


    fclose(fp);

}