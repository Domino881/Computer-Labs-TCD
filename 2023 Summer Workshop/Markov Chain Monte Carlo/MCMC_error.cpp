/*
* Author: Dominik Kuczynski
*
* This code (based heavily on MCMC.cpp) calculates the 
* error of the MCMC method for different bin sizes of the
* final data, in order to approximate the actual (uncorrelated)
* error.
*/

#include<cstdio>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<chrono>
#include<vector>
#include<numeric>
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

const double logN = 20;
const int N = 1<<int(logN);
const int burnin = 128;

int energies[N+3];
int magnets[N+3];

int main(){

    srand(time(NULL));

    short m[xdim][ydim];
    for (int i = 0; i < xdim; i++) {
        for (int j = 0; j < ydim; j++) {
            m[i][j] = rand() > RAND_MAX / 2 ? -1 : 1;
        }

    }

    Lattice l(m);

    printf("%d", l.H());
    l.print();

    int H=0; int M=0;

    double beta = 1.0/5.0;

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

    const int nbins = 30;
    vector<double> binned[nbins];
    double errors[nbins];
    for(int b=0; b<nbins; b++){

        int binsize = (int)pow(2.0, (logN-1.0)/double(nbins-1)*b);
        double cur_avg = 0.0;
        for(int i=burnin+1; i<N; i++){
            cur_avg += (double)energies[i]/binsize;
            if((i-burnin) % binsize == 0){
                binned[b].push_back(cur_avg);
                cur_avg=0.0;
            }
        }

        double mean = reduce(binned[b].begin(), binned[b].end()) / double(binned[b].size());
        double var = 0.0;
        for (auto x : binned[b]){
            var += (x-mean)*(x-mean)/double(binned[b].size());
        }
        int m = (N-burnin) / binsize;
        errors[b] = var/double(m);
    }


    FILE *fp;

    fp = fopen("example.txt", "a");

    // for(int j=0;j<nbins;j++){
    //     fprintf(fp, "%d ", (int)pow(2.0, (logN-1.0)/double(nbins-1)*j));
    // }

    fprintf(fp, "\n");

    for (int i = 0; i < nbins; i++) {
        fprintf(fp, "%f ", errors[i]);
    }

    fclose(fp);

}