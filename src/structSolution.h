#ifndef STRUTCSOLUTION_H
#define STRUTCSOLUTION_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "structBasicsCuts.h"
#include "gpulib/types.h"
#include "lp.h"
#include <time.h>
#include <sys/time.h>
#include <tgmath.h>


EXTERN_C_BEGIN

typedef int TSNumerator;
typedef int TSDenominator;
typedef int TSConstraints;
typedef int TSPosition;
typedef long double TSTempCreateSolution;
typedef double TSViolation;

typedef struct {
    TSNumerator *sNumerator;
    TSDenominator *sDenominator;
    TSConstraints *sConstraints;
}solutionStructCG1;


typedef struct{
    TSNumerator *sNumerator1;
    TSDenominator *sDenominator1;
    TSNumerator *sNumerator2;
    TSDenominator *sDenominator2;
    TSConstraints *sConstraints;
    TSPosition *sPosition;
}solutionStructCG2;

solutionStructCG1* allocStructSolCG1(int nRuns);

solutionStructCG2* allocStructSolCG2(int numberMaxConst, int nRuns);

TCoefficients returnK(TSTempCreateSolution fa_zero);

int verifyDominanceCG(TCoefficients *v1, TRightSide rhs1, TCoefficients *v2, TRightSide rhs2, TNumberVariables sz);

cutFull *createCutsStrongPhaseOne(cutSmall *constraintsSmall, cutFull *constraintsOriginal, solutionStructCG1 *solution, TNumberConstraints nCuts, int precision, int nRuns, TNumberConstraints nC_initial, TCont cont_ini);

cutFull* createCutsStrongPhaseTwo(cutSmall *constraintsSmall, cutFull *constraintsOriginal,  solutionStructCG2 *solution, int numberMaxConst, int nCuts, int precision, int nThreads,int nBlocks);

cutFull* createCutsCover(cutSmall *constraintsSmall, cutFull *constraintsOriginal, cutCover *cutsCover, int *idc_Cover, int nCuts);

int verifyRepeated(cutFull *originalConstraints, int posCover);


EXTERN_C_END

#endif