#ifndef STRUTCBASICSCUTS_H
#define STRUTCBASICSCUTS_H


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "gpulib/types.h"
#include "lp.h"


EXTERN_C_BEGIN

#define INT 1

#ifdef INT
    typedef int TCoefficients;
    typedef int TElements;
    typedef int TElementsConstraints;
    typedef int TRightSide;
    typedef int TXAsterisc;
    typedef int TNumberVariables;
    typedef int TNumberConstraints;
    typedef int TCont;
    typedef short TInterval;
    typedef int TList;
    typedef int TPosList;
#endif 

typedef double TCoefficientsFull;
typedef double TRightSideFull;
typedef double TXAsteriscFull;
typedef char TNameConstraints;
typedef char TNameVariables;
typedef char TActivedCut;



typedef struct {
    TNumberVariables numberVariables;
    TNumberConstraints numberConstraints;
    TCont cont;
    TCoefficients *Coefficients;
    TElements *Elements;
    TElementsConstraints *ElementsConstraints;
    TRightSide *rightSide;
    TXAsterisc *xAsterisc;
    TNumberConstraints *originalConstraints;
}cutSmall;

typedef struct{
    TNumberVariables numberVariables;
    TNumberConstraints numberConstraints;
    TCont cont;
    TCoefficientsFull *Coefficients;
    TElements *Elements;
    TElementsConstraints *ElementsConstraints;
    TRightSideFull *rightSide;
    TXAsteriscFull *xAsterisc;
}cutFull;

typedef struct{
    TNumberConstraints numberCuts;
    TActivedCut activedCut;
    cutFull *pCuts;
}cutPool;

typedef struct{
    TNumberConstraints numberConstraints;
    TCont cont;
    TCoefficients *Coefficients;
    TElementsConstraints *ElementsConstraints;
    TRightSide *rightSide;
}cutCover;


cutSmall *AllocStrCutSmall(TCont cont, TNumberConstraints nConstraints, TNumberVariables nVariables); //allocation struct for GPU (Small)

cutFull *AllocStrCutFull(TCont cont, TNumberConstraints nConstraints, TNumberVariables nVariables); 

cutCover *AllocStrCover(TCont cont, TNumberConstraints nConstraints); //allocation struct Cuts Cover

cutCover *CopyCutToCover(cutSmall *h_cut); //Copy struct of cutSmall for cutCover

//cutFull *fillStructPerLP(LinearProgram *lp, TNameConstraints **nameConstraints, TNameVariables **nameVariables);
cutFull *fillStructPerLP(LinearProgram *lp, TNameConstraints **nameConstraints, TNameVariables **nameVariables, int *typeVariables, double *lbVariables, double *ubVariables);

int countContraintsValided(LinearProgram *lp); //return of number contraints tranformated for <=

int verifyOfFloatIsInteger(TCoefficientsFull coef);

int *returnVectorTypeContraintsIntOrFloat(cutFull *constraints);

cutSmall *reduceCutFullForCutSmall(cutFull *constraints, int *typeIntOrFloat, int precision);

void showStructSmall(cutSmall *constraintsSmall,TNameConstraints **nameConstraints, TNameVariables **nameVariables);

void showStructFull(cutFull *constraintsFull,TNameConstraints **nameConstraints, TNameVariables **nameVariables);

TCoefficients cutMaxDivisorCommonVector(TCoefficients coefs[], TNumberVariables nElem);

TCoefficients cutMaxDivisorCommonRec(TCoefficients m, TCoefficients n);

void freeStrCutFull(cutFull *cut);

int *calcCover(cutCover *h_cover, int *h_solution, int tolerance);

TNameConstraints **renamedNameConstraints(TNameConstraints **nameConstraints, int typeContraints, TNumberConstraints szNewConstraints, TNumberConstraints szCuts,TNumberConstraints lastCut);

cutFull *removeNegativeCoefficientsAndSort(cutFull *constraintsOriginal, int *convertVector, int precision);

void SortByCoefficients(cutFull *h_cut);

void quicksortCof(TCoefficientsFull *values, int *idc, int began, int end);

cutFull *returnVariablesOriginals(cutFull *constraintsOriginal, int *convertVector, int precision, int nVariablesInitial );

void insertConstraintsLP(LinearProgramPtr lp, cutFull *constraintsOriginal, int nConstrainsInitial, char **nameConstraints);


int verifyCutsValidatedPerSolutionInteger(cutFull *constraintsOriginal, int cut, double *sol, char **nameVariables);

double* readSolFile(const char *name, int nVariables);

void quicksortDouble(double *values, int began, int end);

int *returnBinaryConstraints(cutFull *constraintsOriginal, int *typeVariables);

cutFull *convertBinaryConstraints(cutFull *constraintsOriginal, int *BinaryConstraints, int *typeVariables, double *lbVariables, double *ubVariables);

cutFull *convertBinaryOfOriginalConstraints(cutFull *constraintsOriginal, cutFull *constraintsBinary, int nInitialBinary);

EXTERN_C_END
#endif