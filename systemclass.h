//---------------------------------------------------------------------------

#ifndef systemclassH
#define systemclassH
//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string>

using namespace std;
//---------------------------------------------------------------------------
#define PI    3.1415926535897932384626433832795
#define TWOPI 6.28318530717958647692528676655901
#define SQRT2 1.4142135623730950488016887242097

#define EPS 1E-10
#define FPMIN 1E-30

#define ABS(a) ((a)<0.0?(-(a)):(a))
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

//make a b-periodic
#define PBC(AA,BB) {while(AA>=(BB)) AA-=(BB);while(AA<0) AA+=(BB);}

//fuer ran2 (fortran-version)
#define ma 714025
#define ia 1366
#define ic 150889
#define rm 1.0/ma

static const char b32char[]="0123456789ABCDEFGHIJKLMNOPQRSTUV";

void inline int2b4(unsigned int n,unsigned char *b4)
{//convert n to big endian
 b4[0]=(n>>24)&0xFF;
 b4[1]=(n>>16)&0xFF;
 b4[2]=(n>>8)&0xFF;
 b4[3]=n&0xFF;
}

int inline b42int(unsigned char *b4)
{//and back to the machine dependent format
 return (int)b4[3]+(((int)b4[2])<<8)+(((int)b4[1])<<16)+(((int)b4[0])<<24);
}

string IntToHexStr(unsigned int a,int dig=8);

string hexstr(unsigned char *a,int l);

class ESystem {
public:
	ESystem(int size,int dim, double nu,int rs=3234637);
	~ESystem();

	//configuration
	void resample(int rs = 0);//new disorder config
	void reconfig(int type=0); //new electron config
	void reconfigfromfile(string statefilename); //new electron config
	void readSPEfromfile(string spesfilename); //read spes

	void setgriddis(double maxdist); //enables/disables distortions to the square grid - need to use getindist to calculate inverse distances

	//paramter access
	double getdis(int x,int y,int z);
	char getocc(int x,int y,int z);
    char getocci(int i){return n[i];}
    bool getSPE() {return hasSPE;}
	double inline getindist(int x,int y,int z,int p,int dx,int dy,int dz,int &q,bool calcq=false);
    int length() { return L;}
    int size() { return N;} //N=L^D
    int dim() { return D;}
	double getU() {return U;}
	double nu() {return occnum;}
	void setrmax(double r)            //set maximum interaction distance
        {if((r<0) || (r>0.5*L)) r=0.5*L; //due to PBC there shouldn't be a limit, but
            rM=((int) r); if(rM==L/2) rM--; //sometimes the distance between two given points needs to be calculated -> rmax=L/2 !!!
            rmax=rM+0.5;rmaxi=1.0/r;}        //if rM=L/2 distorted grid distances might not be unique!
    void setU(double du) {U=du;}
    void setT(double temp) {T=temp; printf("T set to %f\n", T);}
    void setE(double ext,bool xpc=true) {E=ext;xpbc=xpc;}
	double *getdisptr() {return dis;}
	double *getseptr() {return siteen;}
	double *getspeptr() {return spe;}
	char *getnptr() {return n;}

	void setseptr(double *siteen1){siteen=siteen1;}
	void setspeptr(double *spe1){spe=spe1;}
	void setnptr(char *nptr){n=nptr;}

	int *getposperm() {return posperm;}
    int getrM() {return rM;}

	//energy function
	double getsiteen(int x, int y, int z,bool calc=true,bool getspe=false);
	double calcenergy(); //recalculates the complete energy
	double calcenergyfromSPE();
	double getenergy();  //just sums the current site energies
	double getMinE(int &p,bool usespe=false,bool occc=false,char occs=1); //calculate the min site energy, if occc, n[i]==occs must be fulfilled fro consideration
	double getMaxE(int &p,bool usespe=false,bool occc=false,char occs=1); //...and max
	double calcSPE(int p);
	void setSPE(bool v) {if(v && !hasSPE) updateSPE();hasSPE=v;} //activate SPE

    int count() {int a=0;for(int i=0;i<N;i++) a=a+n[i];return a;} //# of occupied sites

	//hopping helpers
	double flipsite(int x,int y,int z,int p=-1); //flip n[p], here both coordinates and index can be passed

	double hoppEdiff(int i,int j,int k,int a,int b,int c); //energy difference if (i,j,k) <-> (a,b,c)

	double hoppEdiffij(int i,int j); //energy difference if (i,j,k) <-> (a,b,c)
	double hoppEdiffijkl(int i,int j,int k,int l,  int dx, int dy, int dx2, int dy2);
	double hoppEdiffxyz(int i,int j,int k,int a,int b,int c); //energy difference if (i,j,k) <-> (a,b,c)


	double hopp(int i,int j,int k,int a,int b,int c);

	double foursiteXdiff(int p,int q,int i,int j); //get energy difference for 4-site-exchange; only array indices
	double foursiteX(int p,int q,int i,int j);     //do it

	//dynamics
	double NNBoltzstep(); //NN hopping with Boltzman factors, using T
	double densityLE(double dt, int steps); //continous density approach with Langevin equation

	//state characterization
	void getbinstate(unsigned char *state,int len); //get binary state
	int statecorr(unsigned char *s1,unsigned char *s2,int len,int *pos); //get correlation between two states
	string getstatestr(); //gives a base 32 representation of the electron state
	string statestring(); //give a control # for the current state
	string systemID(bool istate=false); //get system id, indep of n-config.

	//some functions
	double ran2(int iseed);
	double absval(double x) {if(x<0) return (-x);return x;}
	double nsgn(double x) {if(x<0.0) return -1.0;return 1.0;}
	void pshuffle();
	int getppos(int i,int &x,int &y,int &z);
	int inline getdistidx(int x,int y,int z,int a,int b,int c);
	void inline getlatdist(int x,int y,int z,int &a,int &b,int &c);
	int pbc(int a) {PBC(a,L) return a;}


	//particle tracing
	void starttrace(int timesteps,int tracenum,int skipsteps);
	int recordpos(int ts);
    int * getrace(int &np,int &steps) {np=traceN;steps=traceL;return motiontrace;}
	void stoptrace();

	//resistance functions
	void fillit(int *usedres,unsigned char *fill,int pos,bool vert);
	double calcresist(double *resists,bool vert);
	bool clusterconnected(unsigned char *cluster,int *usedres);
	double addnextresist(unsigned char *cluster,int *usedresist,double *resval,int &num);    //must have length N,2N,2N
	bool findminresistpath(double *resists,int *path);
      	//double inline calcSPE(int p);

        void nstofile(string filename);
        void spestofile(string filename);
        void movedtofile(string filename);
	void updateSPEi(int i); 
    void addToSiteen(int i, double de){siteen[i]+=de;}
	void setindistDrag(double d);
	double addenergyFromSecondLayer(ESystem *es2);
	double flipsiteDrag(int x,int y,int z,int p,ESystem *es2);
	double hoppDrag(int i,int j,int k,int a,int b,int c,ESystem *es2);
    bool hasdg,hasSPE;  //has displacementgrid = (disgrid!=NULL); change hasSPE to activate parallel SPE calculation. switched by setSPE


private:
 int L,N,D,rM,NP;
 int InitRS,Rseed,tracepos,traceL,traceN,traceskip;
 bool tracepart,xpbc;
 double occnum,rmax,rmaxi,U,T,E;
 char *n;
 double *density; //cont. density
 double *siteen,*spe;
 double *dis,*indist,*indistDrag,*disgrid; //disorder,inverse square lattice distances, displacement grid (dx,dy,dz) for each lattice position, with |da|<1/4 [a=x,y,z]
 int *partnum; //particle labels on the lattice: N
 int *partpos; //last particle positions: NP
 int *posperm; //position permutation, do not write from external program

 int *motiontrace; //traces all particles, huge array: NP*tracelen
 int *sitesMoved;

 int idum;
 int iff,iyy;
 int *ir;

 //support functions
 void iswap(int &a, int &b) { int t=a; a=b; b=t; }
 void qsort(double *val,int *pa, int beg, int end);
 void sortdis(int * pa);

 void updateSPE();  //calculate SPEs from siteen

 void copyntodensity();
 void transformdensityton();

 void resitrace(double *resi,double *minres,int* pred,int pos,bool vert);
 double addsiteenFromSecondLayer(int x, int y, int z, ESystem *es2);

};

#endif
