#include "IsingMC.h"
#include<cmath>

main() {
    long L(128);
    double Kc(log(1+sqrt(2))/2);
    IsingMC m(L);

    m.set_K(0.3, 1<<14);
    m.write("fig4_K3.txt");
    m.zoom_out();
    m.write("fig4_K3z.txt");

    m.set_K(Kc, 1<<14);
    m.write("fig4_Kc.txt");
    m.zoom_out();
    m.write("fig4_Kcz.txt");
    
    m.set_K(0.5, 1<<14);
    m.write("fig4_K5.txt");
    m.zoom_out();
    m.write("fig4_K5z.txt");
}