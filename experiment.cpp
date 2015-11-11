#include <iostream>

#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"

/* Simulated experiment is a K_L->PiPi decay*/

const double K_energy_average /*GeV*/ = 100;
const double K_energy_sigma   /*GeV*/ =  30;
const double K_path_average   /* m */ =  15.34;
const double K_mass           /*GeV*/ =   0.497611;

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

int experiment()
{
    rng = TRandom3(12345); /* Fixed init */
    
    TCanvas* canv_K_energy = new TCanvas("canv_K_energy", "K energy distribution", 1000, 600);
    TH1D* histo_K_energy = new TH1D("histo_K_energy", "K energy distribution; GeV; N",
                                    180, K_energy_average - 3* K_energy_sigma, K_energy_average + 3* K_energy_sigma);
    
    TCanvas* canv_K_path = new TCanvas("canv_K_path", "K path distribution", 1000, 600);
    TH1D* histo_K_path = new TH1D("histo_K_path", "K path distribution; m; N",
                                  180, 0, 10*gamma_from_K_energy(K_energy_average)*K_path_average);
    
    // Main loop
    
    for (int nev=0; nev<N_events; ++nev) {
        double K_energy = generate_K_energy();
        histo_K_energy->Fill(K_energy);
        
        double path_cm = generate_K_path();
        double path = gamma_from_K_energy(K_energy)*path_cm;
        histo_K_path->Fill(path);
    }
    
    // Histogram drawing
    
    canv_K_energy->cd();
    histo_K_energy->Draw();
    
    canv_K_path->cd();
    canv_K_path->SetLogy();
    histo_K_path->Draw();
    
    return 0;
}

int main(int, char**)
{
    return experiment();
}
