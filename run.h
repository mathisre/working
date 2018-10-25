#ifndef RUN_H
#define RUN_H

#include <stdio.h>
#include <math.h>
#include <string.h>

#pragma hdrstop


#include "systemclass.h"
#include "CGanalysis.h"
#include "fileutils.h"
#include "stringutils.h"
#include "paramfile.h"
#include "MKcurrent.h"
#include "params.h"
#include "paramfile.h"


class run
{
public:
    run();
    void FMcb(int step,int Xn,double E,int type,unsigned int secs,void *data);
    void printHistogram(int *histogram, int nrBoxes);
    void printES(ESystem *es);

    void runCurrent(class params *p);
    void runCurrentHeatMap(params *p);
    void runCurrentStateFile(params *p);
    void runCurrentStateSpesFile(params *p);
};

#endif // RUN_H
