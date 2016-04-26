//********************************************************************
// Toy Monte Carlo Simulation of the beta decay detector
// Uses a pre-existing beta decay spectrum
//
// To run, do "root -l bdecay.cpp"
// Make sure you set up your Xterminal so that the graph can be displayed on your screen
//********************************************************************

// C++ libs
#include <iostream>
#include<stdlib.h>	//rand, srand
#include<cmath>

// ROOT libs
#include<TH1D.h>
#include<TFile.h>
#include<TMath.h>
#include<TRandom3>

using namespace std;

////////////////// Parameters ///////////////////////
// Initial nucleus
const int Z_1 = 1;	// Atomic number of initial nucleus (3H)
const double m_1 = 3.0160492;	// Isotope mass (in atomic mass units) 

// Final nucleus
const int Z_2 = 2;	// Atomic number of final nucleus (3He)
const double m_2 = 3.0160293;	// Isotope mass (in atomic mass  units)

// Other parameters
const int charge = -1;
const float Q = 931.5e6*(m_1-m_2);
const int ndivisions = 100;	// Number of divisions in energy histograms
////////////////// End Of Parameters ///////////////

// Physical Constants
const double Pi = 3.14159265;	// Pi
const double alpha = 1./137;	// Structure constant
const float m_p = 938.272046e6;	// Mass (in eV) of proton. From Wikipedia.
const float m_n = 939.5654133e6;	// Mass (in eV) of neutron. From Wikipedia.
const float m_e = 0.511e6;	// Mass (in eV) of electron. From Wikipedia.
const float m_nu= 0.2;	// Mass (in eV) of neutrino. Value was given in our project synopsis

// Functions
float N(float, float, float);		// Distribution of energy, N(T_e)
float F(int, float, int);	// Fermi function, F(Z',T_e)

// Main program
void bdecay2(){
	gStyle->SetOptStat("nemr");	// Makes statistics box appear automatically in histograms

	cout << "Enter the name of the rootfile you will load (e.g. b_decay_histo): ";
	string filename;
	cin >> filename;

	// ROOT rootfile (will contain all histograms)
	TFile *rootfile = new TFile((filename + ".root").c_str(), "read");

	// ROOT Histograms
	TH1D *E_e = (TH1D*)rootfile->Get("E_e");

	// ROOT fit function
	TF1 *func = new TF1("func", "N(x,[0],[1])", 0,Q);
	func->SetParName(0,"m_nu");
	func->SetParameter(0,2);	// We search around a m_nu mass of 2eV (find articles that cite this as upper neutrino mass)
	func->SetParName(1,"C");
	func->SetLineColor(2);	// red

	TCanvas *c1=new TCanvas("E_e","E_e");	// ROOT canvas creation
	E_e->Fit("func","");
	E_e->SetFillColor(4);	//blue
	E_e->Draw();	// Draw histogram
}

// Energy distribution for beta decay
float N(float T_e, float m_nu, float C)
{
	return C*sqrt( pow(T_e,2) + 2*T_e*m_e ) * (T_e + m_e) * (Q-T_e) * sqrt( pow(Q-T_e,2) - pow(m_nu,2) ) * F(Z_2,T_e, charge); // Supposing C=1, taken from http://www2.warwick.ac.uk/fac/sci/physics/research/epp/exp/detrd/amber/betaspectrum/ 
}

// Fermi function
float F(int Z_2, float T_e, int charge)
{
	float eta = (T_e + m_e) * charge * alpha * Z_1 / sqrt(2*T_e*m_e);	// Taken from https://en.wikipedia.org/wiki/Beta_decay#Fermi_function
	return 2. * Pi * eta / (1 - exp(-2*Pi*eta)); // Taken from https://en.wikipedia.org/wiki/Beta_decay#Fermi_function
}
