#include "run.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>

#pragma hdrstop


#include "systemclass.h"
#include "CGanalysis.h"
#include "fileutils.h"
#include "stringutils.h"
#include "paramfile.h"
#include "MKcurrent.h"
#include "params.h"

run::run()
{

}


void run::FMcb(int step,int Xn,double E,int type,unsigned int secs,void *data)
{
    if(step==-2) printf("# find local energy minimum started at CPU time ts=%ds\n",secs);
    else if(step==-1) printf("# energies initialized E=%le, time used t=%ds\n",E,secs);
    else if((step==-11) || (step==-12))
     {
       if(step==-11) printf("# *** state: E=%le, after %d pair exchanges, t=%ds\n",E,Xn,secs);
       else printf("# **** last state was a local minimum\n");
       if(data!=NULL)
        { ESystem *es;
          string fn,s;
          bool hspe;
          int *bins,i;
          FILE *f;
          int bn;
          double dmax,dmin;
          es=(ESystem *) data;
          //analyze and save.
          hspe=es->getSPE();
          es->setSPE(true);
          bn=es->size()/50;  //in average we need 50 counts per bin
          if(bn>500) bn=500; //limit to 500
          bn=10*(bn/10);     //truncate to next decade
          bins=new int[bn];
          dmax=es->getMaxE(i,true);dmax=0.1*((int) (dmax*10));
          dmin=es->getMinE(i,true);dmin=0.1*((int) (dmin*10));
          histogram(es->getspeptr(),es->size(),bins,bn,dmin,dmax);
          es->setSPE(hspe);
          fn=es->systemID()+".dat";
          f=FileOpen(fn,fmOpenAppend);
          //else f=FileCreate(fn);
          s=IntToStr(-11-step)+","+IntToStr(Xn)+","+FloatToStr(E)+","+IntToStr(bn)+","+FloatToStr(dmin)+","+FloatToStr(dmax);
          for(i=0;i<bn;i++) s=s+","+IntToStr(bins[i]);
          s=s+"\n";
          FileWrite(f,s.c_str(),s.length());
          FileClose(f);
          delete[] bins;
        }
     }
    else if(step>=0) printf("#  i=%d Xn=%d E=%le (%d) t=%ds\n",step,Xn,E,type,secs);
    fflush(stdout);
}

//---------------------------------------------------------------------------
void run::printHistogram(int *histogram, int nrBoxes)
{
  int x,y;
  for (x=0;x<nrBoxes/10;x++)
    {
      for (y=0;y<10;y++) printf("%d;",histogram[x*10+y]);
      printf("\n");
    }
}

//---------------------------------------------------------------------------
void run::printES(ESystem *es)
{
  int i,j,k, L, D, N, Lb;
  char *buffer,*n;
  D = es->dim();
  L = es->length();
  N = es->size();
  n = es->getnptr();
  Lb = N/L;
  buffer=new char[N+Lb+4];
  i=0;

  for(j=0;j<Lb;j++)
   {
     memmove(&buffer[i],&n[j*L],L);
     for(k=i;k<i+L;k++)
         {
            //buffer[k]+=0x30; //make it ASCII
            if (buffer[k]==0) buffer[k]=' '; else buffer[k] = '*';
         }
     i+=L;
     buffer[i]='\n';i++;
   }
   buffer[i]=0;
   printf("%s",buffer);
   delete[] buffer;
}


//Type 301 - current
void run::runCurrent(class params *p) //2D only so far!!!
{
//  setlocale(LC_ALL,"");

    //Hva er es?, hva gjør cgcontrol?
//  printf("Running current\n");

  ESystem *es;
  CGcontrol *esc;
  MKcurrent *MKcurr;
  string st;

  double  Ex, Ey, Ez, energy,t, Hx, Hy, Hz; // mc,
  int L, D, N, steps, run, runs; //, im, i, s, f

  printf("\n-----Setting up system-----\n\n");

  steps = p->timesteps;
  D = p->dim;
  L = p->len;
  N = L*L;
  Ex = p->Ex; Ey = p->Ey; Ez = p->Ez;
  Hx = p->Hx; Hy = p->Hy; Hz = p->Hz;

  es = new ESystem(L,D,p->fill,p->rsconfig); // Her skjer det mye !
  es->setrmax(p->screen); // -1.0
  es->setU(p->disU); // 1.0
  es->setT(p->temp);

  es->ran2(p->rsstate); //reset random seed, so initial config and state are independent
  // What is "state"?

  esc = new CGcontrol(); //class for analysis of the systemclass
  esc->setES(es);

  es->reconfig(p->istate); //creates a new (random) initial state
  es->setSPE(true); // ??

  energy = es->calcenergy();
//  printf("Energy: %le\n",energy);

//  printES(es);  // Prints map of initial state (electrons and empty sites)

  MKcurr = new MKcurrent(steps);

  MKcurr->setE(Ex, Ey, Ez);
  MKcurr->setH(Hx, Hy, Hz);
  printf("E set to [%le %le %le]\n",Ex,Ey,Ez);
  printf("H set to [%le %le %le]\n",Hx,Hy,Hz);
  fflush(stdout);

  // maxJL = 5 from param. maxJL = max jump length? MKCurr has both maxJL1 and maxJL2 = 2*maxJL1 +1
  // loc = 3 from param, maxprobability is also 3. JL = Jump Length?
  MKcurr->init(D, L, N, p->maxJL, p->loc, 1.0/p->temp,p->maxProbability, p->Hz, p->doTriJumps); // <- calculates rates
  MKcurr->setES(es);
  MKcurr->setWritelines(p->writelines);
  MKcurr->setRatefun(p->ratefun); // ratefun = 1


  MKcurr->setMTseed(p->seed2); // seed for MT RanGen to be used in getJump, this one uses Mersenne
  printf("MKcurrent initialized\n");
//  fflush(stdout);
  MKcurr->writeGammaToFile(p->Hz);

  // set seed for current:
  es->ran2(p->seed2); // again? Why?

  run = 0; runs = p->runs;
  t = 0;
  printf("\n-------Starting timesteps--------\n\n");
  for (run = 0; run < runs; run++)
  {

    MKcurr->runCurrent(steps, energy, t);    
    printf("\n-------Finished with run %d--------\n\n", run);
    printf("For fil\n");
    fflush(stdout);

    MKcurr->jumpsToFileSmall(p->outputprefix + "jumps_" + p->outputendfix + "_run_" + IntToStr(run) + ".dat",steps);
//    MKcurr->heatMapToFile(   p->outputprefix + "map_" + p->outputendfix + "_run_" + IntToStr(run) + ".dat",steps);

//    MKcurr->normalizeMap(MKcurr->positions, steps);
  //  MKcurr->writeMapToFile(MKcurr->positions, p->outputprefix  + "electronMap_" + p->outputendfix + "_run_" + IntToStr(run) + ".dat");
   // MKcurr->writeMapToFile(MKcurr->movement,  p->outputprefix  + "electronMovement_" + p->outputendfix + "_run_" + IntToStr(run) + ".dat");
    fflush(stdout);

//    MKcurr->currentToFile(p->outputprefix+"r"+IntToStr(run)+"curr.dat",steps);
//    es->nstofile(   p->outputprefix + "final_config_" + p->outputendfix + "_run_" + IntToStr(run) +".dat");
//    es->spestofile( p->outputprefix + "SPE_" + p->outputendfix + "_run_" + IntToStr(run) +".dat");
//    es->movedtofile(p->outputprefix + "moved_particles_" + p->outputendfix + "_run_" + IntToStr(run) +".dat");
    printf("Etter fil\n");
  }
 delete esc;
 delete es;
 delete MKcurr;
}



//Type 302 - current: like 301 but output heat release map
void run::runCurrentHeatMap(params *p) //2D only so far!!!
{
//  setlocale(LC_ALL,"");
  printf("Running current\n");

  ESystem *es;
  CGcontrol *esc;
  MKcurrent *MKcurr;
  string st;

  double mc, Ex, Ey, Ez, e,t;
  int L, D, N, im, i, s, f, steps, run, runs;

  steps=p->timesteps;
  D=p->dim;
  L=p->len;
  N=L*L;
  Ex=p->Ex;
  Ey=p->Ey;
  Ez=p->Ez;

  es=new ESystem(L,D,p->fill,p->rsconfig);
  es->setrmax(p->screen);
  es->setU(p->disU);
  es->setT(p->temp);

 es->ran2(p->rsstate); //reset random seed, so initial config and state are independent

 esc=new CGcontrol();
 esc->setES(es);

 es->reconfig(p->istate); //creates a new (random) initial state
 es->setSPE(true);
 e=es->calcenergy();
 printf("Energy: %le\n",e);


 printES(es);

 MKcurr = new MKcurrent(steps);

 MKcurr->setE(Ex, Ey, Ez);
 printf("E set to [%le %le %le]\n",Ex,Ey,Ez);
 MKcurr->init(D,L,N, p->maxJL, p->loc, 1.0/p->temp,p->maxProbability, p->Hz, p->doTriJumps);
 MKcurr->setES(es);
 MKcurr->setMTseed(p->seed2); // seed for MT RanGen to be used in getJump

 printf("New MKcurrent initialized\n");

  run=0; runs=p->runs;
  t=0;
for (run=0;run<runs;run++)
{
  MKcurr->runCurrent(steps,e,t);
  MKcurr->heatMapToFile(p->outputprefix+"r"+IntToStr(run)+"map.dat",steps);
  //MKcurr->currentToFile(p->outputprefix+"r"+IntToStr(run)+"curr.dat",steps);
  es->nstofile(p->outputprefix+"r"+IntToStr(run)+"ns.dat");
  //es->spestofile(p->outputprefix+"r"+IntToStr(run)+"spes.dat");
}
 delete esc;
 delete es;
 delete MKcurr;
}



//ctype 307, start from a state-file
void run::runCurrentStateFile(params *p) //2D only so far!!! Set seed for current independently
{
//  setlocale(LC_ALL,"");
  printf("Running current\n");

  ESystem *es;
  CGcontrol *esc;
  MKcurrent *MKcurr;
  string st;

  double mc, Ex, Ey, Ez, e1,E,t;
  int L, D, N, im, i, s, f, steps, run, runs;

  steps=p->timesteps;
  D=p->dim;
  L=p->len;
  N=L*L;
  Ex=p->Ex;
  Ey=p->Ey;
  Ez=p->Ez;

  es=new ESystem(L,D,p->fill,p->rsconfig);
  es->setrmax(p->screen);
  es->setU(p->disU);
  es->setT(p->temp);

 es->ran2(p->rsstate); //reset random seed, so initial config and state are independent

 esc=new CGcontrol();
 esc->setES(es);

 //es->reconfig(p->istate); //creates a new (random) initial state
 es->reconfigfromfile(p->statefilename);


 es->setSPE(true);
 e1=es->calcenergy();
 E=e1;
 printf("Energy at starting point: %le\n",e1);

 printES(es);

 MKcurr = new MKcurrent(steps);

 MKcurr->setE(Ex, Ey, Ez);
 printf("Exyz set to [%le %le %le]\n",Ex,Ey,Ez);
 MKcurr->init(D,L,N, p->maxJL, p->loc, 1.0/p->temp,p->maxProbability,p->Hz, p->doTriJumps);
 MKcurr->setES(es);
 MKcurr->setMTseed(p->seed2); // seed for MT RanGen to be used in getJump

 printf("New MKcurrent initialized\n");

 // set seed for current:
  es->ran2(p->seed2);

  run=0; runs=p->runs;
  t=0;
  for (run=0;run<runs;run++)
    {
      MKcurr->runCurrent(steps,E,t);
      MKcurr->jumpsToFileSmall(p->outputprefix+"r"+IntToStr(run)+"jumps.dat",steps);
      //MKcurr->currentToFile(p->outputprefix+"r"+IntToStr(run)+"curr.dat",steps);
      es->nstofile(p->outputprefix+"r"+IntToStr(run)+"ns.dat");
      es->spestofile(p->outputprefix+"r"+IntToStr(run)+"spes.dat");
      es->movedtofile(p->outputprefix+"r"+IntToStr(run)+"mov.dat");
    }
  delete esc;
  delete es;
  delete MKcurr;
}




//ctype 309, start from a state-file and spes-file, so that we need not recalculate all spes.
void run::runCurrentStateSpesFile(params *p) //2D only so far!!! Set seed for current independently
{
//  setlocale(LC_ALL,"");
  printf("Running current\n");

  ESystem *es;
  CGcontrol *esc;
  MKcurrent *MKcurr;
  string st;

  double mc, Ex, Ey, Ez, e1,E,t;
  int L, D, N, im, i, s, f, steps, run, runs;

  steps=p->timesteps;
  D=p->dim;
  L=p->len;
  N=L*L;
  Ex=p->Ex;
  Ey=p->Ey;
  Ez=p->Ez;

  es=new ESystem(L,D,p->fill,p->rsconfig);
  es->setrmax(p->screen);
  es->setU(p->disU);
  es->setT(p->temp);

 es->ran2(p->rsstate); //reset random seed, so initial config and state are independent

 esc=new CGcontrol();
 esc->setES(es);

 //es->reconfig(p->istate); //creates a new (random) initial state
 es->reconfigfromfile(p->statefilename);

 es->readSPEfromfile(p->spesfilename);

 //es->setSPE(true);
 // e1=es->calcenergy();
  e1 = es->calcenergyfromSPE();
 E=e1;
 printf("Energy at starting point: %.15le\n",e1);

 // es->setSPE(true);
 //e1=es->calcenergy();
 //E=e1;
 //printf("Energy at starting point: %le\n",e1);

 // printES(es);

 MKcurr = new MKcurrent(steps);

 MKcurr->setE(Ex, Ey, Ez);
 printf("Exyz set to [%le %le %le]\n",Ex,Ey,Ez);
 MKcurr->init(D,L,N, p->maxJL, p->loc, 1.0/p->temp,p->maxProbability,p->Hz, p->doTriJumps);
 MKcurr->setES(es);
 MKcurr->setMTseed(p->seed2); // seed for MT RanGen to be used in getJump

 printf("New MKcurrent initialized\n");
 // set seed for current:
  es->ran2(p->seed2);

  run=0; runs=p->runs;
  t=0;
  for (run=0;run<runs;run++)
    {
      MKcurr->runCurrent(steps,E,t);
      MKcurr->jumpsToFileSmall(p->outputprefix+"r"+IntToStr(run)+"jumps.dat",steps);
      //MKcurr->currentToFile(p->outputprefix+"r"+IntToStr(run)+"curr.dat",steps);
      es->nstofile(p->outputprefix+"r"+IntToStr(run)+"ns.dat");
      es->spestofile(p->outputprefix+"r"+IntToStr(run)+"spes.dat");
      es->movedtofile(p->outputprefix+"r"+IntToStr(run)+"mov.dat");
    }
  delete esc;
  delete es;
  delete MKcurr;
}






