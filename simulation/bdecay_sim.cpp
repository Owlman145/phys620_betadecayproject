//********************************************************************
// Toy Monte Carlo Simulation of the beta decay detector
//
// To run, do <root -l 'bdecay_sim.cpp("filename")'>
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
//const double Q = 931.494095e6*(m_1-m_2);
const double Q = 18590; // Katrin Q Value (in eV)
const int nevents = 1e7;// Number of events to generate
const double res = 1;	// Resolution of detector (in eV)
const double h = 0.00002;	// Constant used for Von Neuman method. Should range about [1,10]. Too low => cutting distribution, Too high => execution takes too long. Used to estimate the maximum of N(T_e)
double limit=(Q-25)/Q;	// Number between 0 and 1. Fraction of Q over which we necessitate the energy to be 
const int ndivisions = 100;	// Number of divisions in energy histograms
////////////////// End Of Parameters ///////////////

// Physical Constants
const double Pi = 3.14159265;	// Pi
const double alpha = 1./137;	// Structure constant
const double m_p = 938.272046e6;	// Mass (in eV) of proton. From Wikipedia.
const double m_n = 939.5654133e6;	// Mass (in eV) of neutron. From Wikipedia.
const double m_e = 0.510998910e6;	// Mass (in eV) of electron. From Wikipedia.
const double m_nu= 0.2;	// Mass (in eV) of neutrino. Value was given in our project synopsis

// Functions
double N(double, double, double);		// Distribution of energy, N(T_e)
double F(int, double, int);	// Fermi function, F(Z',T_e)

// Main program
void bdecay_sim(string filename){

	limit= limit*Q;	// Limit above which we want our spectrum

	// ROOT random number generator
	TRandom3 *rand = new TRandom3(time(0));	// Generate a random number generator for TRandom3

	// ROOT Histograms
	TH1D *E_e = new TH1D("E_{e}", ";E_{e} [eV];Intensity", ndivisions, limit, Q);	// True kinetic energy histogram for electron
	E_e->SetName("E_e");
	TH1D *E_e_sm = new TH1D("E_{e}", ";E_{e} [eV];Intensity", ndivisions, limit, Q);	// Smeared kinetic energy histogram for electron
	E_e_sm->SetName("E_e_sm");

	// ROOT rootfile (will contain all histograms)
	TFile *rootfile = new TFile((filename + ".root").c_str(), "recreate");

	cout << "(Generating 1e" << log10(1.*nevents) << " events...)\n";
	cout << "Q = " << Q << " eV\n";

	int counter=0;	// Counter for the while loop
	while (counter < nevents) {
		// We use Von Neumann acceptance-rejection method, see phys620 course notes (Monte Carlo p. 21)
		double T_e = rand->Uniform(limit,Q);	// True kinetic energy of the electron. Number between 0 and Q, as there is no energy above Q

		double u = rand->Uniform(1);	// Number between 0 and 1
		if (u <= N(T_e, m_nu, 1.) / (h*N(Q/2, m_nu, 1)))	// "h" is a factor that can be changed so that the sample is more efficient. See phys620 course notes (Monte Carlo p. 21)
		{
			E_e->Fill(T_e);		// Enter true electron kinetic energy in histogram to create beta decay spectrum

			double T_e_sm = rand->Gaus(T_e,res);	// Smeared kinetic energy
			if (Q<=T_e_sm<=limit*Q){
				E_e_sm->Fill(T_e_sm);	// Enter smeared electron kinetic energy in histogram to create beta decay spectrum
			}
			// For execution purposes, acts as a "progress bar"
			if (!(++counter % (nevents/100))) {	// Add 1 to counter and take its modulo, if we finished a 10% of the job
				cout << "Current progress: ";
				cout << 100.*counter/nevents << "%"<< endl;	// Display progress
				if (!(counter % (nevents/10))) {
					cout << "-----" << endl;	// Output "-----" every 10% events
				}
			}
		}
	}
	E_e->Write();	// Save histogram into the rootfile
	E_e_sm->Write();	// Save histogram into the rootfile
}

// Energy distribution for beta decay
double N(double T_e, double m_nu, double C)
{
	return C*sqrt( pow(T_e,2) + 2*T_e*m_e ) * (T_e + m_e) * (Q-T_e) * sqrt( pow(Q-T_e,2) - pow(m_nu,2) ) * F(Z_2,T_e, charge); // Supposing C=1, taken from http://www2.warwick.ac.uk/fac/sci/physics/research/epp/exp/detrd/amber/betaspectrum/ 
}

// Fermi function
double F(int Z_2, double T_e, int charge)
{
	double eta = (T_e + m_e) * charge * alpha * Z_1 / sqrt(2*T_e*m_e);	// Taken from https://en.wikipedia.org/wiki/Beta_decay#Fermi_function
	return 2. * Pi * eta / (1 - exp(-2*Pi*eta)); // Taken from https://en.wikipedia.org/wiki/Beta_decay#Fermi_function
}
