//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string.h>


#pragma hdrstop


#include "systemclass.h"
#include "CGanalysis.h"
#include "fileutils.h"
#include "stringutils.h"
#include "treeutils.h"
#include "paramfile.h"
#include "MKcurrent.h"
#include "params.h"
#include "run.h"
#include "randomc.h"


/* Problem:
 * The integral that depends on the positions of both j and k can't be precalculated
 * Instead, do it numerically via python scipy. Make a script that does this calculation
 * for all combinations. This script should save to a text file. In the c++ script, do an
 * iftest to check if this text file exists. And if it doesn't c++ should call the script
 *
 * These integrals also pop up again but with extra factors in the field dot product!
 * */

int main(int argc, char* argv[])
{
  params *p=new params();

//  printf("Starter lesing av parametre\n");
//  fflush(stdout);
  p->readparams(argc,argv);
//  printf("Slutter lesing av parametre\n");
//  fflush(stdout);
  run runSimul = run();

  //auto start = std::chrono::system_clock::now();
  switch (p->ctype) {
         case 1: runSimul.runCurrent(p);              break;
         case 2: runSimul.runCurrentHeatMap(p);       break;
         case 3: runSimul.runCurrentStateFile(p);     break;
         case 4: runSimul.runCurrentStateSpesFile(p); break;
  }

  delete p;

//  std::chrono::duration<double> diff = end - start;
//  printf("Elapsed time: %.3f min\n",diff.count() / 60.0);

  //int a;scanf("%d",a);
  return 0;
}
