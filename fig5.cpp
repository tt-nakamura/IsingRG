#include "IsingMC.h"
#include<iostream>
#include<fstream>

void fig5(const char *fname, long L,
          double K1, double K2, long n) {
    std::ofstream f(fname);
    double dK((K2-K1)/n),K,U,C;
    long k;
    IsingMC m(L);
    for(k=0; k<=n; k++) {
        K = K1 + k*dK;
        m.set_K(K,1<<10);
        m.get_UC(U,C,1<<16);
        f << K << ' ' << U << ' ' << C << std::endl;
        std::cout << K << ' ' << U << ' ' << C << std::endl;
    }
}

main() {
    fig5("fig5_L32.txt", 32, 0.4, 0.5, 40);
    fig5("fig5_L128.txt", 128, 0.4, 0.5, 40);
}