#include <iostream>
#include <fstream>
#include <cmath>

#include "TH1.h"
//~ #include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TVector3.h"
//~ #include "TLorentzVector.h"

#include "detector.hpp"

/* Simulated experiment is a K_L->PiPi decay (CP violating)*/

const double Pi_mass        /*GeV*/ = 0.139570; // charged pi

const double momentum_z_cut /*GeV*/ = 1;

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
    TCanvas* canv_Pi_pos_reco = new TCanvas("canv_Pi_pos_reco", "Pi pos", 1000, 600);
    TH3D* histo_Pi_pos_reco = new TH3D("histo_Pi_pos_reco", "Pi pos reco; Px (GeV); Py (GeV); Pz (GeV)",
                                       100, -0.3, +0.3, 100, -0.3, +0.3, 100,  0,   +4);
    
    TCanvas* canv_delta_Pi_pos_x = new TCanvas("canv_delta_Pi_pos_x", "Delta Pi pos x", 1000, 600);
    TH1D* histo_delta_Pi_pos_x = new TH1D("histo_delta_Pi_pos_x", "Delta Pi pos x; Dx (GeV); N", 100, -0.05, +0.05);
    
    TCanvas* canv_delta_Pi_pos_y = new TCanvas("canv_delta_Pi_pos_y", "Delta Pi pos y", 1000, 600);
    TH1D* histo_delta_Pi_pos_y = new TH1D("histo_delta_Pi_pos_y", "Delta Pi pos y; Dy (GeV); N", 100, -0.05, +0.05);
    
    TCanvas* canv_delta_Pi_pos_z = new TCanvas("canv_delta_Pi_pos_z", "Delta Pi pos z", 1000, 600);
    TH1D* histo_delta_Pi_pos_z = new TH1D("histo_delta_Pi_pos_z", "Delta Pi pos z; Dz (GeV); N", 100, -0.3,  +0.3);
    
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
        
        // Derivation of initial momentum
        // Uses analytical formulas, information contained in y3,y4 is not used,
        // could be used to over-constrain Py
        const double z12 = z2 - z1;
        const double z34 = z4 - z3;
        
        TVector3 momentum_versor_p12 = momentum_versor_from_2hits(z12, hit_p1, hit_p2);
        TVector3 momentum_versor_p34 = momentum_versor_from_2hits(z34, hit_p3, hit_p4);
        double momentum_z_p = p_kick / (momentum_versor_p34.Px() - momentum_versor_p12.Px());
        TVector3 Pi_pos_reco = momentum_z_p * momentum_versor_p12;
        
        TVector3 momentum_versor_n12 = momentum_versor_from_2hits(z12, hit_n1, hit_n2);
        TVector3 momentum_versor_n34 = momentum_versor_from_2hits(z34, hit_n3, hit_n4);
        double momentum_z_n = - p_kick / (momentum_versor_n34.Px() - momentum_versor_n12.Px());
        TVector3 Pi_neg_reco = momentum_z_n * momentum_versor_n12;
        
        // Cut out events with pz > 1, which are poorly reconstructed due to low delta{x,y}
        if (momentum_z_p > momentum_z_cut || momentum_z_n > momentum_z_cut) continue;
        histo_Pi_pos_reco->Fill(Pi_pos_reco.Px(), Pi_pos_reco.Py(), Pi_pos_reco.Pz());
        histo_delta_Pi_pos_x->Fill(Pi_pos_reco.Px() - Pi_pos.Px());
        histo_delta_Pi_pos_y->Fill(Pi_pos_reco.Py() - Pi_pos.Py());
        histo_delta_Pi_pos_z->Fill(Pi_pos_reco.Pz() - Pi_pos.Pz());

    }
    
    // Histogram drawing
    
    canv_Pi_pos_reco->cd();
    histo_Pi_pos_reco->Draw();
    
    canv_delta_Pi_pos_x->cd();
    histo_delta_Pi_pos_x->Draw();
    
    canv_delta_Pi_pos_y->cd();
    histo_delta_Pi_pos_y->Draw();
    
    canv_delta_Pi_pos_z->cd();
    histo_delta_Pi_pos_z->Draw();
    
    return 0;
}

int main(int, char**)
{
    return analysis();
}
