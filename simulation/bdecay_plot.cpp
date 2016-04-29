//********************************************************************
// Toy Monte Carlo Simulation of the beta decay detector
// Uses a pre-existing beta decay spectrum
//
// To run, do "root -l bdecay2.cpp"
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
const int ndivisions = 100;	// Number of divisions in energy histograms
const double fitmax = Q-25;
const double fitmin = Q-0.2;
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
void gint(TF1*);

// Main program
void bdecay_plot(string filename){
	gStyle->SetOptStat("nemr");	// Makes statistics box appear automatically in histograms

	// ROOT rootfile (will contain all histograms)
	TFile *rootfile = new TFile((filename + ".root").c_str(), "read");

	// ROOT Histograms
	TH1D *E_e = (TH1D*)rootfile->Get("E_e");
	TH1D *E_e_sm = (TH1D*)rootfile->Get("E_e_sm");

	// ROOT fit function
	TF1 *func = new TF1("func", "N(x,[0],[1])",fitmin ,fitmax);
	
	func->SetParName(0,"m_nu");
	func->FixParameter(0,0.2);	// We search around a m_nu mass upper-bounded by 0.3eV (find articles that cite this as upper neutrino mass)
	//func->SetParLimits(0,0.0,0.3);
	//func->SetParLimits(0,0.05,0.3);	// From articles (see bibliography of Jack's article), they estimate the mass to be between 0.05 and 0.3eV)
	func->SetParName(1,"C");
	//func->FixParameter(1,4.46606e-9);
	func->SetLineColor(2);	// red


	TCanvas *c1=new TCanvas("E_e","E_e");	// ROOT canvas creation
	fit = E_e->Fit("func","RMS");
	E_e->SetFillColor(4);	//blue
	E_e->Draw();	// Draw histogram
	cout << "ChiSq = " << fit->Chi2() << endl;

		//double sumundercurve = func->Integral(Q-25,Q);
		//cout << "sumundercurve = " << sumundercurve << endl;

	TCanvas *c2=new TCanvas("E_e_sm","E_e_sm");	// ROOT canvas creation
	func->SetParName(0,"m_nu_sm");
	fit_sm = E_e_sm->Fit("func","RMS");
	E_e_sm->SetFillColor(3);	// green
	E_e_sm->Draw();	// Draw histogram
	cout << "ChiSq = " << fit_sm->Chi2() << endl;
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

void gint(TF1 *g) {
   //default gaus integration method uses 6 points
   //not suitable to integrate on a large domain
   double r1 = g->Integral(0,Q);
   double r2 = g->Integral(Q/2,Q);
   //try with user directives computing more points
   Int_t np = 1000;
   double *x=new double[np];
   double *w=new double[np];
   g->CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
   double r3 = g->IntegralFast(np,x,w,0,5);
   double r4 = g->IntegralFast(np,x,w,0,1000);
   double r5 = g->IntegralFast(np,x,w,0,10000);
   double r6 = g->IntegralFast(np,x,w,0,100000);
   printf("g->Integral(0,5)               = %g\n",r1);
   printf("g->Integral(0,1000)            = %g\n",r2);
   printf("g->IntegralFast(n,x,w,0,5)     = %g\n",r3);
   printf("g->IntegralFast(n,x,w,0,1000)  = %g\n",r4);
   printf("g->IntegralFast(n,x,w,0,10000) = %g\n",r5);
   printf("g->IntegralFast(n,x,w,0,100000)= %g\n",r6);
   delete [] x;
   delete [] w;
}

