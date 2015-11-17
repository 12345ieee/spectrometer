#include <iostream>
#include <fstream>
#include <cmath>

//~ #include "TH1.h"
//~ #include "TH2.h"
//~ #include "TH3.h"
//~ #include "TCanvas.h"
//~ #include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
//~ #include "TLorentzVector.h"

/* Simulated experiment is a K_L->PiPi decay (CP violating)*/

const double Pi_mass          /*GeV*/ =  0.139570; // charged pi

namespace MC // Contains things not known to experimenters
{
    const double K_energy_average /*GeV*/ =  2; // Around the value of Cronin-Fitch experiment
    const double K_energy_sigma   /*GeV*/ =  0.5;
    const double K_path_average   /* m */ = 15.34;
    const double K_mass           /*GeV*/ =  0.497611;
    
    const double Pi_energy_cm = K_mass/2;
    const double Pi_momentum_mod_cm = sqrt(Pi_energy_cm*Pi_energy_cm - Pi_mass*Pi_mass);

    /* Detector smearing */
    const double dx /* m */ =  0.01;
    const double dy /* m */ =  0.01;
}

/* Detector coordinates */
const double z1 /* m */ = 200;
const double z2 /* m */ = 201;
const double z3 /* m */ = 250;
const double z4 /* m */ = 251;

/* dp kick - applied at the center of the field region */
const double z_drift /* m */ = z3-z2;
const double z_kick /* m */ = (z2+z3)/2;
const double p_kick /* GeV */ = 1e-5*z_drift;

/* Plots constants */
const double hits_bound /* m */ = 100;

using namespace std;

ostream& operator<<(ostream& ost, TVector2& vec)
{
    return ost << vec.X() << " " << vec.Y();
}

ostream& operator<<(ostream& ost, TVector3& vec)
{
    return ost << vec.X() << " " << vec.Y() << " " << vec.Z();
}

istream& operator>>(istream& ist, TVector2& vec)
{
    double x, y;
    ist >> x >> y;
    vec = TVector2(x, y);
    return ist;
}

istream& operator>>(istream& ist, TVector3& vec)
{
    double x, y, z;
    ist >> x >> y >> z;
    vec = TVector3(x, y, z);
    return ist;
}

inline TVector2 line_propagation(TVector3 Pi_momentum, double z_end, double z_start, double x_start=0, double y_start=0)
{
    double z_travelled = z_end - z_start;
    double scale = z_travelled/Pi_momentum.Pz();
    return TVector2(Pi_momentum.Px()*scale + x_start, Pi_momentum.Py()*scale + y_start);
}

inline TVector3 momentum_versor_from_2hits(double deltaz, TVector2 hit1, TVector2 hit2)
{
    double px = (hit2.X() - hit1.X())/deltaz;
    double py = (hit2.Y() - hit1.Y())/deltaz;
    return TVector3(px, py, 1);
}

int analysis()
{
    // Histograms declaration
    
    // In file
    ifstream infile;
    infile.open("experiment.txt");
    
    int N_events = 0;
    
    // Main loop
    while (true) {
        
        // Check file
        int nev;
        infile >> nev;
        if (!infile.good()) break;
        N_events++;
        
        // Declare vars
        TVector2 hit_p1;
        TVector2 hit_n1;
        TVector2 hit_p2;
        TVector2 hit_n2;
        TVector2 hit_p3;
        TVector2 hit_n3;
        TVector2 hit_p4;
        TVector2 hit_n4;
        double K_energy;
        double path;
        TVector3 Pi_pos;
        TVector3 Pi_neg;
        
        // Read from file
        infile >> hit_p1;
        infile >> hit_n1;
        infile >> hit_p2;
        infile >> hit_n2;
        infile >> hit_p3;
        infile >> hit_n3;
        infile >> hit_p4;
        infile >> hit_n4;
        infile >> K_energy;
        infile >> path;
        infile >> Pi_pos;
        infile >> Pi_neg;
    }
    
    // Histogram drawing
    
    return 0;
}

int main(int, char**)
{
    return analysis();
}
