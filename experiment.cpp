#include <iostream>

#include "TH1.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"

/* Simulated experiment is a K_L->PiPi decay*/

const double K_energy_average /*GeV*/ = 100;
const double K_energy_sigma   /*GeV*/ =  30;
const double K_path_average   /* m */ =  15.34;
const double K_mass           /*GeV*/ =   0.497611;

const double Pi_mass          /*GeV*/ =   0.139570; // charged pi

const double Pi_energy_cm = K_mass/2;
const double Pi_momentum_mod_cm = sqrt(Pi_energy_cm*Pi_energy_cm - Pi_mass*Pi_mass);

const int N_events = 1e5;

TRandom rng;

using namespace std;

inline double generate_K_energy()
{
    return max(K_mass, rng.Gaus(K_energy_average, K_energy_sigma));
}

inline double generate_K_path()
{
    return rng.Exp(K_path_average);
}

inline double gamma_from_K_energy(double K_energy)
{
    return K_energy/K_mass;
}

inline TVector3 generate_Pi_momentum_cm()
{
    double x,y,z;
    rng.Sphere(x, y, z, Pi_momentum_mod_cm);
    return TVector3(x, y, z);
}

int experiment()
{
    rng = TRandom3(12345); /* Fixed init */
    
    TCanvas* canv_K_energy = new TCanvas("canv_K_energy", "K energy distribution", 1000, 600);
    TH1D* histo_K_energy = new TH1D("histo_K_energy", "K energy distribution; GeV; N",
                                    180, K_energy_average - 3* K_energy_sigma, K_energy_average + 3* K_energy_sigma);
    
    TCanvas* canv_K_path = new TCanvas("canv_K_path", "K path distribution", 1000, 600);
    TH1D* histo_K_path = new TH1D("histo_K_path", "K path distribution; m; N",
                                  180, 0, 10*gamma_from_K_energy(K_energy_average)*K_path_average);
    
    TCanvas* canv_Pi_momentum_cm = new TCanvas("canv_Pi_momentum_cm", "Pi momentum cm", 1000, 600);
    TH3D* histo_Pi_momentum_cm = new TH3D("histo_Pi_momentum_cm", "Pi momentum cm; GeV; GeV; GeV",
                                          100, -Pi_momentum_mod_cm, +Pi_momentum_mod_cm,
                                          100, -Pi_momentum_mod_cm, +Pi_momentum_mod_cm,
                                          100, -Pi_momentum_mod_cm, +Pi_momentum_mod_cm);

    TCanvas* canv_Pi1 = new TCanvas("canv_Pi1", "Pi1", 1000, 600);
    TH3D* histo_Pi1 = new TH3D("histo_Pi1", "Pi1; GeV; GeV; GeV",
                                          100,                   0, +K_energy_average,
                                          100, -Pi_momentum_mod_cm, +Pi_momentum_mod_cm,
                                          100, -Pi_momentum_mod_cm, +Pi_momentum_mod_cm);
    
    // Main loop
    
    for (int nev=0; nev<N_events; ++nev) {
        double K_energy = generate_K_energy();
        histo_K_energy->Fill(K_energy);
        
        double path_cm = generate_K_path();
        double path = gamma_from_K_energy(K_energy)*path_cm;
        histo_K_path->Fill(path);
        
        double K_beta_z = sqrt(1-pow(K_mass/K_energy, 2));
        TVector3 K_beta = TVector3(0, 0, K_beta_z);
        TVector3 Pi_momentum_cm = generate_Pi_momentum_cm();
        histo_Pi_momentum_cm->Fill(Pi_momentum_cm.Px(), Pi_momentum_cm.Py(), Pi_momentum_cm.Pz());
        
        TLorentzVector Pi1 = TLorentzVector( Pi_momentum_cm, Pi_energy_cm);
        TLorentzVector Pi2 = TLorentzVector(-Pi_momentum_cm, Pi_energy_cm); // Momentum cons. in cm
        Pi1.Boost(K_beta);
        Pi2.Boost(K_beta);
        histo_Pi1->Fill(Pi1.Px(), Pi1.Py(), Pi1.Pz());
    }
    
    // Histogram drawing
    
    canv_K_energy->cd();
    histo_K_energy->Draw();
    
    canv_K_path->cd();
    canv_K_path->SetLogy();
    histo_K_path->Draw();
    
    canv_Pi_momentum_cm->cd();
    histo_Pi_momentum_cm->Draw();
    
    canv_Pi1->cd();
    histo_Pi1->Draw();
    
    return 0;
}

int main(int, char**)
{
    return experiment();
}
