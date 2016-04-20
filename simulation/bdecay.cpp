//********************************************************************
// Toy Monte Carlo Simulation of the beta decay detector
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
const int Z_1 = 6;	// Atomic number of initial nucleus
const double m_1 = 14.003241;	// Isotope mass (in atomic mass units) 

// Final nucleus
const int Z_2 = 7;	// Atomic number of final nucleus
const double m_2 = 14.0030740052;	// Isotope mass (in atomic mass  units)

// Other parameters
const int nevents = 1e5;// Number of events to generate
const float h = 10;	// Constant used for Von Neuman method. Should range about [1,10]. Too low => cutting distribution, Too high => execution takes too long. Used to estimate the maximum of N(T_e)
int ndivisions = 50;	// Number of divisions in energy histograms
////////////////// End Of Parameters ///////////////

// Physical Constants
const double Pi = 3.14159265;	// Pi
const double alpha = 1./137;	// Structure constant
const float m_p = 938.272046;	// Mass (in MeV) of proton. From Wikipedia.
const float m_n = 939.5654133;	// Mass (in MeV) of neutron. From Wikipedia.
const float m_e = 0.511;	// Mass (in MeV) of electron. From Wikipedia.
const float m_nu= 0.2e-6;	// Mass (in MEV) of neutrino. Value was given in our project synopsis

// Functions
float N(float, int, float);		// Distribution of energy, N(T_e)
float F(int, float, int);	// Fermi function, F(Z',T_e)
float m(float, float);	// Semi-empirical mass formula (Bertulani eq 5.8)

// Main program
void bdecay(){
	// Parameter-dependent values
	const int charge = Z_1 - Z_2;		// Charge of electron/positron (in units of e)
	float Q;	// Q-value of our experiment
	if (charge == -1)	// If negative beta decay...
	{
		Q = 931.5*(m_1-m_2);	// Taken from Bertulani eq 8.9. The 931.5 factor is to convert amu unit to MeV unit.
	}
	else if (charge == +1)	// If positive beta decay...
	{
		Q = 931.5*(m_1-m_2)-2m_e;	// Taken from Bertulani eq 8.9. The 931.5 factor is to convert amu unit to MeV unit.
	}
	else	// This is not a beta decay
	{
		break;
	}

	gStyle->SetOptStat("nemr");	// Makes statistics box appear automatically in histograms

	// ROOT random number generator
	TRandom3 *rand = new TRandom3();	// Generate a random number generator for TRandom3

	// ROOT Histograms
	TH1D *E_e = new TH1D("E_{e}", ";E_{e} [MeV];Intensity", ndivisions, 0, 1.1*Q);	// Energy histogram for electron

	// ROOT rootfile (will contain all histograms)
	TFile *rootfile = new TFile("beta_decay_histograms.root", "recreate");

	cout << "(Generating 1e" << log10(1.*nevents) << " events...)\n";
	cout << "Q = " << Q << " MeV\n";

	int counter=0;	// Counter for the while loop
	while (counter < nevents) {
		// We use Von Neumann acceptance-rejection method, see phys620 course notes (Monte Carlo p. 21)
		double T_e = rand->Uniform(Q);	// Number between 0 and Q, as there is no energy above Q
		double u = rand->Uniform(1);	// Number between 0 and 1
		if (u <= N(T_e, charge, Q) / (h*N(Q/2, charge, Q)))	// "h" is a factor that can be changed so that the sample is more efficient. See phys620 course notes (Monte Carlo p. 21)
		{
			E_e->Fill(T_e);	// Enter electron kinetic energy in histogram to create beta decay spectrum

			// For execution purposes, acts as a "progress bar"
			if (!(++counter % (nevents/10))) {	// Add 1 to counter and take its modulo
				cout << "Current progress: " << 100.*counter/nevents << "%"<< endl;	// Display progress
			}
		}
	}

	E_e->Write();	// Save histogram into the rootfile
	TCanvas *c1=new TCanvas("E_e","E_e");	// ROOT canvas creation
	E_e->SetFillColor(4);	//blue
	E_e->Draw();	// Draw histogram
}

// Energy distribution for beta decay
float N(float T_e, int charge, float Q)
{
	return sqrt( pow(T_e,2) + 2*T_e*m_e ) * (T_e + m_e) * (Q-T_e) * sqrt( pow(Q-T_e,2) - pow(m_nu,2) ) * F(Z_2,T_e, charge); // Supposing C=1, taken from http://www2.warwick.ac.uk/fac/sci/physics/research/epp/exp/detrd/amber/betaspectrum/ 
}

// Fermi function
float F(int Z_2, float T_e, int charge)
{
	float eta = (T_e + m_e) * charge * alpha * Z_1 / sqrt(2*T_e*m_e);	// Taken from https://en.wikipedia.org/wiki/Beta_decay#Fermi_function
	return 2. * Pi * eta / (1 - exp(-2*Pi*eta)); // Taken from https://en.wikipedia.org/wiki/Beta_decay#Fermi_function
}

// Semi-empirical mass formula (Bertulani eq 5.8). Returns mass in MeV
float m(float Z, float A)
{
	// Constants (in MeV) for binding energy B. Taken from Bertulani, eq 5.9
	float a_v = 15.85;
	float a_s = 18.34;
	float a_c = 0.71;
	float a_a = 92.86;
	float a_p = 1.46;

	// Binding energy formula, taken from Bertulani eq 5.8
	float B = a_v*A - a_s*pow(1.*A, 2./3.) - a_c*pow(Z,2)*pow(1.*A,-1./3.) - a_a*pow(Z-1.*A/2.,2)/A - 0.5*(pow(-1,Z)+pow(-1,A-Z))*a_p*pow(1.*A,-1./2);

	return Z*m_p + (A-Z)*m_n - B;	// From Bertulani, eq 4.9
	
}
