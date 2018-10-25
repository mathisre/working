//---------------------------------------------------------------------------

#ifndef MKcurrentH
#define MKcurrentH
//---------------------------------------------------------------------------
#include "systemclass.h"
#include "treeutils.h"
#include "randomc.h"
#include "mersenne.h"
#include "vector"

using namespace std;

#define MAXFINDMIN 1000

//class for current within a systemclass


class MKcurrent {

public:
  MKcurrent(int steps);
  ~MKcurrent();

	void setES(ESystem *e);

    void init(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1, double Bz, bool doTriJumps);

    void setE(double Ex1, double Ey1, double Ez1);
    void setH(double Hx1, double Hy1, double Hz1);
	void setE0(double Ex1, double Ey1, double Ez1);
	void setOmega(double omega1);
	void setWritelines(int writelines1);
	void setRatefun(int ratefun1);

	void runCurrent(int steps, double &E, double &t);
	void runCurrentAC(int steps, double &E, double &t, int &dx);

	//	double getDeltat(){return deltat;}

	//	void toFile(string filename, int steps);
	void jumpsToFileSmall(string filename, int steps);
	void jumpsToFileAC(string filename, int steps);
	void heatMapToFile(string filename, int steps);
	void setMTseed(int seed);
    void updateMap();
    void updateMovement(int i, int j);
    void normalizeMap(std::vector<std::vector<double> > &vector, int steps);
    void writeMapToFile(std::vector<std::vector<double>> vector, string filename);
    void writeGammaToFile(double Hz);
    CRandomMersenne *RanGen;



    std::vector<std::vector<double>> positions = std::vector<std::vector<double>>(N, std::vector<double>(N,0));
    std::vector<std::vector<double>> movement = std::vector<std::vector<double>>(N, std::vector<double>(N,0));

private:
    //	double inline getJL(int dx, int dy, int dz);
    //      savedstate inline *checkstateexistance();
    void inline getjump(int &i, int &j, int &dx, int &dy);
    void inline getjump3sites(int &i, int &j, int &k, int &dx, int &dy, int &dxI, int &dyI);
    bool inline testjump(int i, int j, int dx, int dy);
    bool inline testjump3sites(int i, int j, int k, int dx, int dy, int dxI, int dyI);

	//	double inline calcRate2d(int n1, int n2, int dx, int dy);

	int maxJL, JL2;
	double *jl; //jumpLength
    double *GammaT,*GammaT3Sites; //Tunneling Gamma (see Tsiganskov, Pazy, Laikhmann, Efros
    double GammaTtot, maxProbability, GammaTtot3Sites, totGammaTtot;
    int Nmem,Nmem3Sites, nGamma3Sites;
    int *dxInter,*dyInter,*dxFinal,*dyFinal;
    int testedNumberOf2Site, testedNumberOf3Site, numberOf2Site, numberOf3Site;
	double A, beta; //2/a, for the localization...

	double deltat;

	int *from;
	int *to;
	double *energy,*ts,*des;
	int *MCsteps;
	int *dxs;
    int *dys;

    double Ex, Ey, Ez, E0x, E0y, E0z, Hx, Hy, Hz;
	double tMC,omega;
	int writelines,ratefun;
    double meanJumpLength;
    double meandx, meandy;
    int L,N,D;
	//	char *n;

	ESystem *es;

};

#endif
