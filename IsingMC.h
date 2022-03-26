#ifndef __IsingMC_h__
#define __IsingMC_h__

// Monte Carlo simulation of Ising model in two dimension
struct IsingMC {
    long L;// number of spins in one direction
    long N;// L*L = total number of spins
    char **s;// matrix of spin values (1 or -1)
    double K;// J/(kT), T = temperature, J = spin interaction
    long E;// total energy / J
    double p[3];// transition probability
    long *ip,*im;// indices of nearest neighbor
    IsingMC(long);
    ~IsingMC();
    void set_K(double, long=1);
    void get_UC(double&, double&, long=8, long=1);
    void run(long);
    double correlation(long);
    void zoom_out();
    void write(const char*);
};

#endif // __IsingMC_h__