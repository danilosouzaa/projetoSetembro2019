#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

extern "C"
{
#include "lp.h"
#include "prepareCpu.h"
#include "structSolution.h"
#include "structBasicsCuts.h"
}

int main(int argc, const char *argv[])
{
    int i;
    int cg1 = 0, cg2 = 0, cc = 0;
    if (argc < 15)
    {
        printf("Number of parameters invalided\n");
        printf("[1] - Instance Name \n");
        printf("[2] - Precison of xAsterisc\n");
        printf("[3] - Time Maximal of test \n");
        printf("[4] - Max Coef Phase 1\n");
        printf("[5] - Size Group Phase 2\n");
        printf("[6] - number of runs Phase 2\n");
        printf("[7] - max Denominator of Phase 2\n");
        printf("[8-11] - weight's common,  of Phase 2\n");
        printf("[12] - number of cc tests\n");
        printf("[13] - number iteration Local Search Grasp CC\n");
        printf("[14] - alpha GRASP\n");
        printf("[15] - minimal\n");
        return 0;
    }
    for (i = 0; i < argc; i++)
    {
        if (strcmp(argv[i], "-CG1") == 0)
        {
            cg1 = 1;
        }
        if (strcmp(argv[i], "-CG2") == 0)
        {
            cg2 = 1;
        }
        if (strcmp(argv[i], "-CC") == 0)
        {
            cc = 1;
        }
    }
    char nameFileInstance[255] = "../inst/";
    char nameInst[255] = "";
    int precision = atoi(argv[2]);
    double timeMax = atof(argv[3]);
    int limitCoefPhase1 = atoi(argv[4]);
    int sizeGroupConstraintsPhase2 = atoi(argv[5]);
    int nRunsPhase2 = atoi(argv[6]);
    int maxDenominator = atoi(argv[7]);
    float weightVariablesCommon = atof(argv[8]);
    float weightVariablesCoef = atof(argv[9]);
    float weightVariablesSignal = atof(argv[10]);
    float weightVariablesFrac = atof(argv[11]);
    int szPoolCutsMaxCC = atoi(argv[12]);
    int nIterationCCGrasp = atoi(argv[13]);
    float alpha = atof(argv[14]);
    int  minimal = atoi(argv[15]);// 0 non minimal - 1 minimal
    
    strcat(nameInst, argv[1]);
    strcat(nameFileInstance, argv[1]);

    LinearProgram *lp = lp_create();
    lp_read(lp, nameFileInstance);
    //---------------------------------------------------------------------------------------
    //-------------Allocation strings of nameContraints and nameVariables--------------------
    //---------------------------------------------------------------------------------------
    int nVariables = lp_cols(lp);
    int nConstraintsValided = countContraintsValided(lp);
    TNameConstraints **nameConstraints = (TNameConstraints **)malloc(nConstraintsValided * sizeof(TNameConstraints *));
    for (i = 0; i < nConstraintsValided; i++)
    {
        nameConstraints[i] = (TNameConstraints *)malloc(255 * sizeof(TNameConstraints));
    }
    TNameVariables **nameVariables = (TNameVariables **)malloc(nVariables * sizeof(TNameVariables *));
    for (i = 0; i < nVariables; i++)
    {
        nameVariables[i] = (TNameVariables *)malloc(255 * sizeof(TNameVariables));
        lp_col_name(lp, i, nameVariables[i]);
    }
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
//--------------------------------------fill structs-------------------------------------
//---------------------------------------------------------------------------------------

// for(i=0;i<szSol;i++){
//      printf("%lf\n",sol[i]);
//      if(sol[i]!=0.0)
//         getchar();
//  }
//#ifdef DEBUG
    lp_set_max_solutions(lp, 1);
    lp_optimize(lp);
    lp_write_sol(lp, "../sol/solutionTest.sol");
    double *sol = readSolFile("../sol/solutionTest.sol", nVariables);
    for (i=0;i<nVariables;i++){
        printf("x_%d = %f\t", i+1,sol[i]);
    }
    printf("\n");

//#endif // DEBUG
    int *typeVariables = (int *)malloc(sizeof(int) * nVariables);
    double *lbVariables = (double *)malloc(sizeof(double) * nVariables);
    double *ubVariables = (double *)malloc(sizeof(double) * nVariables);
    cutFull *constraintsOriginal = fillStructPerLP(lp, nameConstraints, nameVariables, typeVariables, lbVariables, ubVariables);
#ifdef DEBUG

    // for (i = 0; i < constraintsOriginal->numberVariables; i++)
    // {
    //     if (sol[i] != 0)
    //     {
    //         printf("%d %s = %lf\n", i, nameVariables[i], sol[i]);
    //     }
    // }
    // getchar();
#endif // DEBUG
    TNumberConstraints numberConstraintsInitial = constraintsOriginal->numberConstraints;
    TCont numberContInitial = constraintsOriginal->cont;
    TNumberVariables numberVariablesInitial;
    int numberCutsCG1 = 0, numberCutsCG2 = 0, numberCutsCC = 0, numberAux;
    //showStructFull(constraintsOriginal, nameConstraints, nameVariables);
    lp_write_lp(lp, "teste.lp");
    //getchar();
    double startT = omp_get_wtime();
    double _time = 0;
    _time = ((double)timeMax - (omp_get_wtime() - startT));
    TNumberConstraints numberAuxConstraints;
    TNumberConstraints totalCuts;
    lp_optimize_as_continuous(lp);
    double iniObjSol = lp_obj_value(lp);
    printf("Solution initial: %f \n", iniObjSol);

    while (_time > 1)
    {
//showStructSmall(constraintsSmall,nameConstraints,nameVariables);

//---------------------------------------------------------------------------------------
//-------------------------Call of methods CG and Cover----------------------------------
//---------------------------------------------------------------------------------------
//showStructFull(constraintsOriginal,nameConstraints,nameVariables);
#ifdef DEBUG

        for (i = 0; i < constraintsOriginal->numberConstraints; i++)
        {
            int very = verifyCutsValidatedPerSolutionInteger(constraintsOriginal, i, sol, nameVariables);

            if (very == 0)
            {
                //getchar();
                printf("validado antes: %d %s\n", very, nameConstraints[i]);
            }
        }

#endif // DEBUG
        numberAuxConstraints = constraintsOriginal->numberConstraints;

        if (cg1 == 1)
        {
            numberConstraintsInitial = constraintsOriginal->numberConstraints;
            numberContInitial = constraintsOriginal->cont;

            numberAux = constraintsOriginal->numberConstraints;
            constraintsOriginal = runCG1_mainCpu(constraintsOriginal, precision, timeMax, limitCoefPhase1, numberConstraintsInitial, numberContInitial);
            numberAux = constraintsOriginal->numberConstraints - numberAux;
            nameConstraints = renamedNameConstraints(nameConstraints, 1, constraintsOriginal->numberConstraints, numberAux, numberCutsCG1);
            numberCutsCG1 += numberAux;
#ifdef DEBUG
            for (i = 0; i < constraintsOriginal->numberConstraints; i++)
            {
                int very = verifyCutsValidatedPerSolutionInteger(constraintsOriginal, i, sol, nameVariables);

                if (very == 0)
                {
                    //getchar();
                    printf("validado depois: %d %s\n", very, nameConstraints[i]);
                }
            }
#endif
        }

        if (cg2 == 1)
        {
            numberAux = constraintsOriginal->numberConstraints;
            constraintsOriginal = runCG2_mainCpu(constraintsOriginal, precision, timeMax, sizeGroupConstraintsPhase2, nRunsPhase2, maxDenominator, 0, weightVariablesCommon, weightVariablesCoef, weightVariablesSignal, weightVariablesFrac);
            numberAux = constraintsOriginal->numberConstraints - numberAux;
            nameConstraints = renamedNameConstraints(nameConstraints, 2, constraintsOriginal->numberConstraints, numberAux, numberCutsCG2);
            numberCutsCG2 += numberAux;
        }
        if (cc == 1)
        {
            int *binaryConstraints = returnBinaryConstraints(constraintsOriginal, typeVariables); // verifica quais restrições é possível usar no método
            int constraintsUsed = 0;
            for (i = 0; i < constraintsOriginal->numberConstraints; i++)
            {
                if (binaryConstraints[i] != 0)
                {
                    constraintsUsed++;
                }
                //printf("binaryConstraints: %d\n", binaryConstraints[i]);
            }
         //   printf("Number Constraints Posible: %d\n", constraintsUsed);

            cutFull *constraintsBinary = convertBinaryConstraints(constraintsOriginal, binaryConstraints, typeVariables, lbVariables, ubVariables);
            printf("Constraints Original: %d \n new Constraints: %d\n", constraintsOriginal->numberConstraints, constraintsBinary->numberConstraints);
            int nInitialBinary = constraintsBinary->numberConstraints;
            int *convertVariables = (int *)malloc(sizeof(int) * constraintsBinary->cont);
            numberVariablesInitial = constraintsBinary->numberVariables;
            constraintsBinary = removeNegativeCoefficientsAndSort(constraintsBinary, convertVariables, precision);
            // int j,el_test;
            // double lhs_test;

            /*          for (i = 0; i < constraintsBinary->numberConstraints; i++)
            {
                lhs_test = 0;
                for (j = constraintsBinary->ElementsConstraints[i]; j < constraintsBinary->ElementsConstraints[i + 1]; j++)
                {   
                    el_test = constraintsBinary->Elements[j];
                    lhs_test += constraintsBinary->Coefficients[j] * constraintsBinary->xAsterisc[el_test];
                    printf("%lf x%d = %lf %d + ", constraintsBinary->Coefficients[j], constraintsBinary->Elements[j], constraintsBinary->xAsterisc[el_test], typeVariables[el_test]);
                }
                printf("<= %lf\n", constraintsBinary->rightSide[i]);
                if(lhs_test - 1e-5 >constraintsBinary->rightSide[i]){
                    printf("ERRO na complementação: %lf %lf\n", lhs_test, constraintsBinary->rightSide[i]);
                    getchar();
                }else{
                    printf("OKKK: %lf %lf\n", lhs_test, constraintsBinary->rightSide[i]);
                }

            }*/
            // showStructFull(constraintsBinary,nameConstraints,nameVariables);
            // getchar();
            // for(i=0;i<constraintsBinary->numberVariables;i++){
            //     printf("%s %f %f-- ", nameVariables[i],lbVariables[i], ubVariables[i]);
            // }
            // getchar();

            numberAux = constraintsBinary->numberConstraints;
#ifdef DEBUG
            constraintsBinary = runCC_mainCPuDebug(constraintsBinary, precision, nameConstraints, nameVariables, sol, szPoolCutsMaxCC, nIterationCCGrasp, alpha, minimal);
#else
            constraintsBinary = runCC_mainCPu(constraintsBinary, precision, szPoolCutsMaxCC);
#endif // DEBUG

            // for (i = 0; i < constraintsBinary->numberConstraints; i++)
            // {
            //     printf("%d: ", i);
            //     for (j = constraintsBinary->ElementsConstraints[i]; j < constraintsBinary->ElementsConstraints[i + 1]; j++)
            //     {
            //         printf("%lf x%d + ", constraintsBinary->Coefficients[j], constraintsBinary->Elements[j]);
            //     }
            //     printf("<= %lf\n", constraintsBinary->rightSide[i]);
            // }
            // getchar();

            numberAux = constraintsBinary->numberConstraints - numberAux;

            constraintsBinary = returnVariablesOriginals(constraintsBinary, convertVariables, precision, numberVariablesInitial);

            constraintsOriginal = convertBinaryOfOriginalConstraints(constraintsOriginal, constraintsBinary, nInitialBinary);

            nameConstraints = renamedNameConstraints(nameConstraints, 3, constraintsOriginal->numberConstraints, numberAux, numberCutsCC);
            numberCutsCC += numberAux;

            freeStrCutFull(constraintsBinary);
            //showStructFull(constraintsOriginal,nameConstraints,nameVariables);
            //getchar();
            //showStructFull(constraintsOriginal,nameConstraints,nameVariables);
            free(convertVariables);
            free(binaryConstraints);
        }

        int *verifyTest = (int *)malloc(sizeof(int) * constraintsOriginal->numberConstraints);
        for (i = 0; i < constraintsOriginal->numberConstraints; i++)
        {
            verifyTest[i] = verifyCutsValidatedPerSolutionInteger(constraintsOriginal, i, sol, nameVariables);

            if (verifyTest[i] == 0)
            {
                printf("validado depois: %d %s\n", verifyTest[i], nameConstraints[i]);
                //getchar();
            }
        }

        _time = ((double)timeMax - (omp_get_wtime() - startT));
        totalCuts = constraintsOriginal->numberConstraints - numberAuxConstraints;
        printf("Cuts total: %d\n", totalCuts);
        printf("Depois: %d \n", constraintsOriginal->numberConstraints);
        if (totalCuts >= 0)
        {
            insertConstraintsLP(lp, constraintsOriginal, numberAuxConstraints, nameConstraints, verifyTest);
            lp_write_lp(lp, "danilo.lp");
            lp_optimize_as_continuous(lp);
            double *xTemp = lp_x(lp);
            for (i = 0; i < constraintsOriginal->numberVariables; i++)
            {
                constraintsOriginal->xAsterisc[i] = xTemp[i];
            }

            printf("value: %f\n", lp_obj_value(lp));
        }
        free(verifyTest);
    }
    //  showStructFull(constraintsOriginal,nameConstraints,nameVariables);

    //---------------------------------------------------------------------------------------
    //-------------------------free of Allocation struct-------------------------------------
    //---------------------------------------------------------------------------------------

    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        free(nameConstraints[i]);
    }
    free(nameConstraints);

    for (i = 0; i < constraintsOriginal->numberVariables; i++)
    {
        free(nameVariables[i]);
    }
    free(nameVariables);
#ifdef DEBUG
    free(sol);
#endif // DEBUG
    free(typeVariables);
    free(lbVariables);
    free(ubVariables);
    freeStrCutFull(constraintsOriginal);
    _time = (omp_get_wtime() - startT);
    printf("Time final: %f\n", _time);
    //---------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------
    return 0;
}
