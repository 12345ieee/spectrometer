#include <iostream>
#include <fstream>
#include <cmath>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

/* Simulated experiment is a K_L->PiPi decay (CP violating)*/

const double K_energy_average /*GeV*/ =  2; // Around the value of Cronin-Fitch experiment
const double K_energy_sigma   /*GeV*/ =  0.5;
const double K_path_average   /* m */ = 15.34;
const double K_mass           /*GeV*/ =  0.497611;

const double Pi_mass          /*GeV*/ =  0.139570; // charged pi

const double Pi_energy_cm = K_mass/2;
const double Pi_momentum_mod_cm = sqrt(Pi_energy_cm*Pi_energy_cm - Pi_mass*Pi_mass);

/* Detector coordinates */
const double z1 /* m */ = 200;
const double z2 /* m */ = 201;
const double z3 /* m */ = 250;
const double z4 /* m */ = 251;

/* Detector smearing */
const double dx /* m */ =  0.01;
const double dy /* m */ =  0.01;

/* dp kick - applied at the center of the field region */
const double z_drift /* m */ = z3-z2;
const double z_kick /* m */ = (z2+z3)/2;
const double p_kick /* GeV */ = 1e-5*z_drift;

/* Plots constants */
const double hits_bound /* m */ = 100;

const int N_events = 1e5;

TRandom rng;

using namespace std;

ostream& operator<<(ostream& ost, TVector2& vec)
{
    return ost << vec.X() << " " << vec.Y();
}

ostream& operator<<(ostream& ost, TVector3& vec)
{
    return ost << vec.X() << " " << vec.Y() << " " << vec.Z();
}

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

inline TVector2 line_propagation(TVector3 Pi_momentum, double z_end, double z_start, double x_start=0, double y_start=0)
{
    double z_travelled = z_end - z_start;
    double scale = z_travelled/Pi_momentum.Pz();
    return TVector2(Pi_momentum.Px()*scale + x_start, Pi_momentum.Py()*scale + y_start);
}

inline TVector2 apply_chamber_smearing(TVector2 hit, double dx, double dy)
{
    double x = rng.Gaus(hit.X(), dx);
    double y = rng.Gaus(hit.Y(), dy);
    return TVector2(x, y);
}

int experiment()
{
    rng = TRandom3(12345); /* Fixed init */
    
    // Histograms declaration
    TCanvas* canv_K_energy = new TCanvas("canv_K_energy", "K energy distribution", 1000, 600);
    TH1D* histo_K_energy = new TH1D("histo_K_energy", "K energy distribution; E (GeV); N",
                                    180, K_energy_average - 3* K_energy_sigma, K_energy_average + 3* K_energy_sigma);
    
    TCanvas* canv_K_path = new TCanvas("canv_K_path", "K path distribution", 1000, 600);
    TH1D* histo_K_path = new TH1D("histo_K_path", "K path distribution; Path (m); N",
                                  180, 0, 10*gamma_from_K_energy(K_energy_average)*K_path_average);
    
    TCanvas* canv_Pi_momentum_cm = new TCanvas("canv_Pi_momentum_cm", "Pi momentum cm", 1000, 600);
    TH3D* histo_Pi_momentum_cm = new TH3D("histo_Pi_momentum_cm", "Pi momentum cm; Px (GeV); Py (GeV); Pz (GeV)",
                                          100, -Pi_momentum_mod_cm, +Pi_momentum_mod_cm,
                                          100, -Pi_momentum_mod_cm, +Pi_momentum_mod_cm,
                                          100, -Pi_momentum_mod_cm, +Pi_momentum_mod_cm);
    
    TCanvas* canv_Pi_pos = new TCanvas("canv_Pi_pos", "Pi pos", 1000, 600);
    TH3D* histo_Pi_pos = new TH3D("histo_Pi_pos", "Pi pos; Px (GeV); Py (GeV); Pz (GeV)",
                                  100, -Pi_momentum_mod_cm, +Pi_momentum_mod_cm,
                                  100, -Pi_momentum_mod_cm, +Pi_momentum_mod_cm,
                                  100, 0,                   +2*K_energy_average);
    
    TCanvas* canv_hit_p1 = new TCanvas("canv_hit_p1", "Hit pos ch1", 1000, 600);
    TH2D* histo_hit_p1 = new TH2D("histo_hit_p1", "Hit pos ch1; x (m); y (m); N",
                                  200, -hits_bound, +hits_bound, 200, -hits_bound, +hits_bound);
    
    TCanvas* canv_hit_p1_smeared = new TCanvas("canv_hit_p1_smeared", "Hit pos ch1 smeared", 1000, 600);
    TH2D* histo_hit_p1_smeared = new TH2D("histo_hit_p1_smeared", "Hit pos ch1 smeared; x (m); y (m); N",
                                          200, -hits_bound, +hits_bound, 200, -hits_bound, +hits_bound);
    
    TCanvas* canv_hit_p3 = new TCanvas("canv_hit_p13", "Hit pos ch3", 1000, 600);
    TH2D* histo_hit_p3 = new TH2D("histo_hit_p3", "Hit pos ch3; x (m); y (m); N",
                                  200, -hits_bound, +hits_bound, 200, -hits_bound, +hits_bound);
    
    // Out file
    ofstream outfile;
    outfile.open("experiment.txt");
    
    // Main loop
    for (int nev=0; nev<N_events; ++nev) {
        
        // Generate K energy
        double K_energy = generate_K_energy();
        histo_K_energy->Fill(K_energy);
        
        // Generate K path
        double path_cm = generate_K_path();
        double path = gamma_from_K_energy(K_energy)*path_cm;
        histo_K_path->Fill(path);
        
        // Discard K that have not decayed before the first chamber
        if (path > z1) continue;
        
        // Generate cm dynamic
        TVector3 Pi_momentum_cm = generate_Pi_momentum_cm();
        histo_Pi_momentum_cm->Fill(Pi_momentum_cm.Px(), Pi_momentum_cm.Py(), Pi_momentum_cm.Pz());
        
        // Boost to lab frame
        double K_beta_z = sqrt(1 - pow(K_mass/K_energy, 2));
        // double x = K_mass/K_energy;
        // double K_beta_z = 1 + pow(x, 2)/2 - pow(x, 4)/8 + pow(x, 6)/16;
        TVector3 K_beta = TVector3(0, 0, K_beta_z);
        TLorentzVector Pi_pos_4v = TLorentzVector( Pi_momentum_cm, Pi_energy_cm);
        TLorentzVector Pi_neg_4v = TLorentzVector(-Pi_momentum_cm, Pi_energy_cm); // Momentum cons. in cm
        Pi_pos_4v.Boost(K_beta);
        Pi_neg_4v.Boost(K_beta);
        TVector3 Pi_pos = Pi_pos_4v.Vect();
        TVector3 Pi_neg = Pi_neg_4v.Vect();
        histo_Pi_pos->Fill(Pi_pos.Px(), Pi_pos.Py(), Pi_pos.Pz());
        
        // Calculate hit points in first two chambers
        TVector2 hit_p1=line_propagation(Pi_pos, z1, path);
        TVector2 hit_n1=line_propagation(Pi_neg, z1, path);
        TVector2 hit_p2=line_propagation(Pi_pos, z2, path);
        TVector2 hit_n2=line_propagation(Pi_neg, z2, path);
        histo_hit_p1->Fill(hit_p1.X(), hit_p1.Y());
        
        // Apply chamber smearing
        TVector2 hit_p1_smeared = apply_chamber_smearing(hit_p1, dx, dy);
        TVector2 hit_n1_smeared = apply_chamber_smearing(hit_n1, dx, dy);
        TVector2 hit_p2_smeared = apply_chamber_smearing(hit_p2, dx, dy);
        TVector2 hit_n2_smeared = apply_chamber_smearing(hit_n2, dx, dy);
        histo_hit_p1_smeared->Fill(hit_p1_smeared.X(), hit_p1_smeared.Y());
        
        // The magnetic field is schematized as a kick of magnitude p_kick
        // applied at the center of the magnetic field region (z_kick)
        
        // Find position at kick
        TVector2 Pi_pos_position_kick = line_propagation(Pi_pos, z_kick, path);
        TVector2 Pi_neg_position_kick = line_propagation(Pi_neg, z_kick, path);
        TVector3 p_kick_vec = TVector3(p_kick, 0, 0);
        
        // Find hits in last two chambers - p_kick is opposite due to charge
        TVector2 hit_p3=line_propagation(Pi_pos + p_kick_vec, z3, z_kick, Pi_pos_position_kick.X(), Pi_pos_position_kick.Y());
        TVector2 hit_n3=line_propagation(Pi_neg - p_kick_vec, z3, z_kick, Pi_neg_position_kick.X(), Pi_neg_position_kick.Y());
        TVector2 hit_p4=line_propagation(Pi_pos + p_kick_vec, z4, z_kick, Pi_pos_position_kick.X(), Pi_pos_position_kick.Y());
        TVector2 hit_n4=line_propagation(Pi_neg - p_kick_vec, z4, z_kick, Pi_neg_position_kick.X(), Pi_neg_position_kick.Y());
        histo_hit_p3->Fill(hit_p3.X(), hit_p3.Y());
        
        // Apply chamber smearing
        TVector2 hit_p3_smeared = apply_chamber_smearing(hit_p3, dx, dy);
        TVector2 hit_n3_smeared = apply_chamber_smearing(hit_n3, dx, dy);
        TVector2 hit_p4_smeared = apply_chamber_smearing(hit_p4, dx, dy);
        TVector2 hit_n4_smeared = apply_chamber_smearing(hit_n4, dx, dy);
        
        
        // Write to file
        outfile << nev            << " ";
        outfile << hit_p1_smeared << " ";
        outfile << hit_n1_smeared << " ";
        outfile << hit_p2_smeared << " ";
        outfile << hit_n2_smeared << " ";
        outfile << hit_p3_smeared << " ";
        outfile << hit_n3_smeared << " ";
        outfile << hit_p4_smeared << " ";
        outfile << hit_n4_smeared << " ";
        outfile << K_energy       << " ";
        outfile << path           << " ";
        outfile << Pi_pos         << " ";
        outfile << Pi_neg         <<endl;
    }
    
    // Histogram drawing
    
    canv_K_energy->cd();
    histo_K_energy->Draw();
    
    canv_K_path->cd();
    canv_K_path->SetLogy();
    histo_K_path->Draw();
    
    canv_Pi_momentum_cm->cd();
    histo_Pi_momentum_cm->Draw();
    
    canv_Pi_pos->cd();
    histo_Pi_pos->Draw();
        
    canv_hit_p1->cd();
    histo_hit_p1->Draw();
    
    canv_hit_p1_smeared->cd();
    histo_hit_p1_smeared->Draw();
    
    canv_hit_p3->cd();
    histo_hit_p3->Draw();
    
    return 0;
}

int main(int, char**)
{
    return experiment();
}
