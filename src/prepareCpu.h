#ifndef PREPARECPU_H
#define PREPARECPU_H

#include "structSolution.h"
#include <sys/time.h>
#include <omp.h>
#include "lp.h"

EXTERN_C_BEGIN

typedef float TParametersGroup;

cutSmall *removeXAsteriscEqualOne(cutSmall *constraintsSmall, int precision);

cutFull *runCG1_mainCpu(cutFull *constraintsFull, int precision, double timeconst, int limitCoefMAx, TNumberConstraints numberConstraintsInitial, TCont contInitial);

void runCG_Phase_1_cpu(cutSmall *constraintsSmall, solutionStructCG1 *solutionCG1, int precision, double timeconst, int limitCoefMax);

cutFull *runCG2_mainCpu(cutFull *constraintsFull, int precision, double timeLimit, TNumberConstraints numberMaxConst, int nRuns, int maxDenominator, TSPosition posMaior, float pCom, float pCoef, float pSignal, float pFrac);

void runCG_phase_2_cpu(cutSmall *constraintsSmall, solutionStructCG2 *solutionCG2, int precision, double timeLimit, TNumberConstraints *setConstraint, TNumberConstraints numberMaxConst,  int nRuns, int maxDenominator);

void quicksortTParameters(TParametersGroup *values, TNumberConstraints *idc, TNumberConstraints began, TNumberConstraints end);

TNumberConstraints *returnVectorParametersConstraints(cutSmall *constraintsSmall, int precision, TNumberConstraints Constraint, TNumberConstraints numberMaxConst, TParametersGroup pCom, TParametersGroup pCoef, TParametersGroup pSignal, TParametersGroup pFrac);

cutFull *runCC_mainCPu(cutFull *constraintsFull, int precision, int szCoverThread);

void createSolutionsInitial(int *solution, int sz);


void createInitialCoverGRASP(int *solution, int sz, cutSmall *constraintsSmall, int precision, int constraint, float alpha);

//void createCoverGRASP(int *solution, int sz, cutSmall *constraintsSmall, int precision);

int *createCoverGraspIndividual(cutSmall *constraintsSmall, int precision, TNumberConstraints constraint, int numberIteration, int szPoolCutsMax, float alpha, int minimal);

void copyAndVerifyPoolSolution(int *solution, int sz, int *poolSolution, int *numberSolutionAtual);

int *localSearch(int *solution, int sz, cutSmall *constraintsSmall, int precision, TNumberConstraints constraint);

#ifdef DEBUG
//cutFull *runCC_mainCPuDebug(cutFull *constraintsFull, int precision, int szCoverThread, char **nameConstraints, char **nameVariables, double *sol);
cutFull *runCC_mainCPuDebug(cutFull *constraintsFull, int precision, char **nameConstraints, char **nameVariables, double *sol, int szPoolCutsMax, int nIterationGrasp, float alpha, int minimal);
#endif // DEBUG

EXTERN_C_END

#endif