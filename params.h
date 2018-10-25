#ifndef PARAMS_H
#define PARAMS_H
//---------------------------------------------------------------------------
#include <stdio.h>
#include <string>
#include <math.h>

#pragma hdrstop


#include "systemclass.h"
#include "CGanalysis.h"
#include "fileutils.h"
#include "stringutils.h"
#include "paramfile.h"
#include "MKcurrent.h"
#include "params.h"
#include "run.h"
using namespace std;
//---------------------------------------------------------------------------
#define maxistates 1000
#define maxtpos 100

class params
{

public:	// Anwender-Deklarationen
  params() {tpos=new int[maxtpos];}
 ~params() {delete[] tpos;}

 int rsconfig,rsstate,dim,len,ctype,tracenum,tposn;
 int timesteps,skipsteps,eonsteps,runs,seed2,startrun;
 double disU,fill,fillPassive,temp,inittemp,extf,screen,lattdis,cutoffexp,maxProbability,limitOfNoReturn,starttime;
 int istate,htype,rhopps;
 bool running,wantstop,trace,xpbc;
 //lowT
 int maxJL,maxJL2el,d2el, maxExpMax, maxExpMin;
 double loc, Ex, Ey, Ez,omega,d, Hx, Hy, Hz;
 int initruns, writelines,ratefun;
 string stopfilename,outputprefix, outputendfix, statefilename, spesfilename, statefilenamePassive, spesfilenamePassive;
 bool doTriJumps;
 int *tpos;

public:		// Anwender-Deklarationen

        void readparams(int n, char* cmdline[]);
};

#endif
