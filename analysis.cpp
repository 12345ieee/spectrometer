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

inline TVector3 Pi_versor_reco_from_2hits(double deltaz, TVector2 hit1, TVector2 hit2)
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
    TH1D* histo_delta_Pi_pos_x = new TH1D("histo_delta_Pi_pos_x", "Delta Pi pos x; Dx (GeV); N", 200, -0.05, +0.05);
    
    TCanvas* canv_delta_Pi_pos_y = new TCanvas("canv_delta_Pi_pos_y", "Delta Pi pos y", 1000, 600);
    TH1D* histo_delta_Pi_pos_y = new TH1D("histo_delta_Pi_pos_y", "Delta Pi pos y; Dy (GeV); N", 200, -0.05, +0.05);
    
    TCanvas* canv_delta_Pi_pos_z = new TCanvas("canv_delta_Pi_pos_z", "Delta Pi pos z", 1000, 600);
    TH1D* histo_delta_Pi_pos_z = new TH1D("histo_delta_Pi_pos_z", "Delta Pi pos z; Dz (GeV); N", 200, -0.3,  +0.3);
    
    TCanvas* canv_z_path_reco = new TCanvas("canv_z_path_reco", "K path reco distribution", 1000, 600);
    TH1D* histo_z_path_reco = new TH1D("histo_z_path_reco", "K path reco distribution; Path (m); N", 200, 0, z1);
    
    TCanvas* canv_delta_z_path = new TCanvas("canv_delta_z_path", "Delta K path distribution", 1000, 600);
    TH1D* histo_delta_z_path = new TH1D("histo_delta_z_path", "Delta K path distribution; Path (m); N", 200, -0.8, 0.8);
    
    TCanvas* canv_K_energy_reco = new TCanvas("canv_K_energy_reco", "K energy reco distribution", 1000, 600);
    TH1D* histo_K_energy_reco = new TH1D("histo_K_energy_reco", "K energy reco distribution; E (GeV); N", 200, 0, 3);
    
    TCanvas* canv_delta_K_energy = new TCanvas("canv_delta_K_energy", "Delta K energy distribution", 1000, 600);
    TH1D* histo_delta_K_energy = new TH1D("histo_delta_K_energy", "Delta K energy distribution; E (GeV); N", 200, -0.5, 0.5);
    
    TCanvas* canv_K_mass_reco = new TCanvas("canv_K_mass_reco", "K mass reco distribution", 1000, 600);
    TH1D* histo_K_mass_reco = new TH1D("histo_K_mass_reco", "K mass reco distribution; E (GeV); N", 200, 0.4, 0.6);
    
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
        
        TVector3 Pi_versor_reco_p12 = Pi_versor_reco_from_2hits(z12, hit_p1, hit_p2);
        TVector3 Pi_versor_reco_p34 = Pi_versor_reco_from_2hits(z34, hit_p3, hit_p4);
        double momentum_z_p = p_kick / (Pi_versor_reco_p34.Px() - Pi_versor_reco_p12.Px());
        TVector3 Pi_pos_reco = momentum_z_p * Pi_versor_reco_p12;
        
        TVector3 Pi_versor_reco_n12 = Pi_versor_reco_from_2hits(z12, hit_n1, hit_n2);
        TVector3 Pi_versor_reco_n34 = Pi_versor_reco_from_2hits(z34, hit_n3, hit_n4);
        double momentum_z_n = - p_kick / (Pi_versor_reco_n34.Px() - Pi_versor_reco_n12.Px());
        TVector3 Pi_neg_reco = momentum_z_n * Pi_versor_reco_n12;
        
        // Cut out events with pz > 1, which are poorly reconstructed due to low delta{x,y}
        if (momentum_z_p > momentum_z_cut || momentum_z_n > momentum_z_cut) continue;
        histo_Pi_pos_reco->Fill(Pi_pos_reco.Px(), Pi_pos_reco.Py(), Pi_pos_reco.Pz());
        histo_delta_Pi_pos_x->Fill(Pi_pos_reco.Px() - Pi_pos.Px());
        histo_delta_Pi_pos_y->Fill(Pi_pos_reco.Py() - Pi_pos.Py());
        histo_delta_Pi_pos_z->Fill(Pi_pos_reco.Pz() - Pi_pos.Pz());
        
        // Compute the best approximation of the decay point
        // The point is forced to lay on the z axis and computed minimizing
        // chi2(zP) = d(P, track_p)^2 + d(P, track_n)^2
        
        // FIXME: the exact analytical solution of the problem is used
        // may want to use a numerical method
        double mod3_p = Pi_versor_reco_p12.Mag2();
        double mod3_n = Pi_versor_reco_n12.Mag2();
        double z_path_reco = z1 - (mod3_n*Pi_versor_reco_p12.XYvector()*hit_p1 + mod3_p*Pi_versor_reco_n12.XYvector()*hit_n1)/
                                  (mod3_n*Pi_versor_reco_p12.XYvector().Mod2() + mod3_p*Pi_versor_reco_n12.XYvector().Mod2());
        histo_z_path_reco->Fill(z_path_reco);
        histo_delta_z_path->Fill(z_path_reco - path);
        
        // Construction of K_energy and K_mass
        double Pi_pos_energy_reco = sqrt(Pi_pos_reco.Mag2() + Pi_mass*Pi_mass);
        double Pi_neg_energy_reco = sqrt(Pi_neg_reco.Mag2() + Pi_mass*Pi_mass);
        
        double K_energy_reco = Pi_pos_energy_reco + Pi_neg_energy_reco;
        histo_K_energy_reco->Fill(K_energy_reco);
        histo_delta_K_energy->Fill(K_energy_reco - K_energy);
        
        double K_mass_reco = sqrt(2*Pi_mass*Pi_mass + 2*Pi_pos_energy_reco*Pi_neg_energy_reco - 2*Pi_pos_reco*Pi_neg_reco);
        histo_K_mass_reco->Fill(K_mass_reco);
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
    
    canv_z_path_reco->cd();
    canv_z_path_reco->SetLogy();
    histo_z_path_reco->Draw();
    
    canv_delta_z_path->cd();
    histo_delta_z_path->Draw();
    
    canv_K_energy_reco->cd();
    histo_K_energy_reco->Draw();
    
    canv_delta_K_energy->cd();
    histo_delta_K_energy->Draw();
    
    canv_K_mass_reco->cd();
    histo_K_mass_reco->Fit("gaus");
    histo_K_mass_reco->Draw();
    
    return 0;
}

int main(int, char**)
{
    return analysis();
}
