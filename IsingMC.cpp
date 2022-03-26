#include<cstdlib>
#include<cmath>
#include<fstream>
#include "IsingMC.h"

// uniform random number in [0,1)
#define frand() (rand()/(RAND_MAX+1.))

IsingMC::IsingMC(long L1) :
// L1 = number of spins in one direction
// temperature is set to zero
// spins are all initialized to align upward
// N = number of spins is set to L1*L1
L(L1), N(L*L), K(1e99), E(-2*N),
s(new char*[L]), ip(new long[L]), im(new long[L])
{
    long i;
    s[0] = new char[N];
    for(i=1; i<L; i++) s[i] = s[i-1] + L;
    for(i=0; i<N; i++) s[0][i] = 1;
    // periodic boundary
    for(i=0; i<L; i++) ip[i] = (i+1)%L;
    for(i=0; i<L; i++) im[i] = (i-1+L)%L;
}

IsingMC::~IsingMC() {
    delete[] s[0];
    delete[] s;
    delete[] ip;
    delete[] im;
}

void IsingMC::set_K(double K1, long n)
// set K = J/(kT), T = temperature, J = spin interaction
// n = number of Monte-Carlo steps to equilibrate / N
{
    long i;
    K = K1;
    for(i=0; i<3; i++) p[i] = exp(-4*i*K);// Boltzmann factor
    run(n*N);// run until equilibrium
}

void IsingMC::get_UC(double& U, double& C, long m, long n)
// U = mean longernal energy per spin
// C = specific heat per spin (variance of energy)
// m = sample size to compute mean and variance
// n = number of Monte-Carlo steps between samples / N
{
    long k,nN(n*N);
    U=C=0;
    for(k=0; k<m; k++) {
        run(nN);
        U += E;
        C += pow(E,2);
    }
    U /= m;
    C = (C/m - U*U)*K*K;
    U /= N;
    C /= N;
}

void IsingMC::run(long n)
// run simulation by Metropolis algorithm
// n = number of Monte-Carlo steps
{
    long i,j,k,e;
    for(k=0; k<n; k++) {
        i = (long)(frand()*N);// choose a spin at random
        j = i%L;
        i /= L;
        e = s[ip[i]][j] + s[im[i]][j];
        e+= s[i][ip[j]] + s[i][im[j]];
        e*= s[i][j];// change in energy / 2
        if(e <= 0 || frand() < p[e/2]) {// Metropolis
            s[i][j] = -s[i][j];// flip spin
            E += 2*e;
        }
    }
}

double IsingMC::correlation(long l)
// l = distance between two spins
// return <s(0)s(l)>
{
    long i,j;
    double c(0);
    for(i=0; i<L; i++)
        for(j=0; j<L; j++)
            c += s[i][j] * (s[(i+l)%L][j] + s[i][(j+l)%L]);
    return c /= 2*N;
}

void IsingMC::zoom_out()// Kadanoff transformation
// lattice is shrinked (zoomed out) by a factor of sqrt(2),
//   rotated by +45 deg, and half of spins
//   are removed by majority rule in neiboring spins
// assume L is even
{
    long i,j,k,l,L2(L/2);
    char t[L][L],u;
    for(i=0; i<L; i++)
        for(j=0; j<L; j++) t[i][j] = s[i][j];
    for(i=0; i<L; i++) {
        for(j=i%2; j<L; j+=2) {
            u = t[ip[i]][j] + t[im[i]][j];
            u+= t[i][ip[j]] + t[i][im[j]] + t[i][j];
            u = (u<0 ? -1:1);// majority rule
            k = (i+j)/2;// rotate by +45 deg
            l = (i-j+L)/2;
            s[k][l] = u;
            if(k<L2) k+=L2; else k-=L2;
            if(l<L2) l+=L2; else l-=L2;
            s[k][l] = u;// fill by periodicity
        }
    }
}

void IsingMC::write(const char *fname) {
    long i,j;
    std::ofstream f(fname);
    for(i=0; i<L; i++) {
        f << (int)s[i][0];
        for(j=1; j<L; j++) f << ' ' << (int)s[i][j];
        f << std::endl;
    }
}