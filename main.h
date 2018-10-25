//---------------------------------------------------------------------------
#ifndef mainH
#define mainH
//---------------------------------------------------------------------------
#include <stdio.h>
#include <string>

using namespace std;
//---------------------------------------------------------------------------
#define maxistates 1000
#define maxtpos 100

class params
{

public:	// Anwender-Deklarationen
  params() {printf("lager ny params\n"); fflush(stdout);  tpos=new int[maxtpos];}
 ~params() {delete[] tpos;}
 int rsconfig,rsstate,dim,len,ctype,tracenum,tposn;
 int timesteps,skipsteps,eonsteps,runs,seed2,startrun;
 double disU,fill,fillPassive,temp,inittemp,extf,screen,lattdis,cutoffexp,maxProbability,limitOfNoReturn,starttime;
 int istate,htype,rhopps;
 bool running,wantstop,trace,xpbc;
 //lowT
 int maxJL,maxJL2el,d2el, maxExpMax, maxExpMin;
 double loc, Ex, Ey, Ez,omega,d;
 int initruns, writelines,ratefun;
 string stopfilename,outputprefix, statefilename, spesfilename, statefilenamePassive, spesfilenamePassive;
 int *tpos;

public:		// Anwender-Deklarationen
		void readparams(int n, char* cmdline[]);
};

double findemin(ESystem *es,int size,int hopptyp,int &im,int rhopps=4);

//---------------------------------------------------------------------------
#endif
