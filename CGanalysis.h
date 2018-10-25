//---------------------------------------------------------------------------

#ifndef CGanalysisH
#define CGanalysisH
//---------------------------------------------------------------------------
#include "systemclass.h"
using namespace std;

#define MAXFINDMIN 1000

typedef void(*FMcallback)(int,int,double,int,unsigned int, void *data); //callback function for findemin: parameters: step,#of pair exchanges,typ,time

//class for analysis of the systemclass

class CGcontrol {

public:
	CGcontrol() {es=NULL;}

	void setES(ESystem *e);

	double findemin(int hopptyp,int rhopps,int &im, FMcallback callback=NULL);

	void addeecorr(double *ca,int cl,double avE); //calculate the siteenergy pair correlations
	bool dipol(double &Px,double &Py,double &Pz); //calculate dipol moment of the system

	//additional function: FFT
	void drealft(double *data, unsigned long n, int isign);

	unsigned int getsecs();  //

	//begin{MK}
	double CGfindemin(int hopptyp,int &im);
	double CGrandomizedfindemin(int hopptyp,int &im);
	double findStepInCG();
	double checkAllSPjumps();
	double optimizeCG();
	void makeHistogram( int nrBoxes, int *histogram, int type, double highestE, double lowestE); //type=1: occ, 2: uocc, 3:all
	double checkMicrocrystallinity(string filename);
	void saveHistogram(int *histogram, int nrBoxes, string filename);
	void printCorrelation(int *correlation);
	void saveCorrelation(int *correlation, string filename);
	//end{MK}

	//some compatibilty functions, es is assumed to be valid
	double getMINadden() {int p;return es->getMinE(p,true,false);};//For use in making histograms
	double getMAXadden() {int p;return es->getMaxE(p,true,false);}; //For use in making histograms
	double getsiteadden(int x, int y, int z,bool calc=true) {return es->getsiteen(x,y,z,calc,true);};
	void calcaddenergy() {es->setSPE(true);}; //Addition energies (E(occupied)-E(free))
	double calcbothenergies() {double E=es->calcenergy();es->setSPE(true);return E;}; //sum of siteen+ siteen&siteadden set
	double getMAXoccen(int &p) {return es->getMaxE(p,true,true,1);}; //Returns the energy of the highest occupied site

private:
  int L,N,D;
  double nu,U;
  char *n;
  double *dis,*se;
  ESystem *es;

  int pbc(int a) {PBC(a,L) return a;}
  
 //FFT part
  void dfour1(double *data, unsigned long nn, int isign);

};

int histogram(double *data,int L,int *bins,int N,double dmin, double dmax);


#endif
