#include "prepareCpu.h"

cutSmall *removeXAsteriscEqualOne(cutSmall *constraintsSmall, int precision)
{
    int i, j;
    TNumberConstraints nNewConstraints = 0, szConst, el;
    TCont newCont = 0;
    TNumberConstraints *constraintsValided = (TNumberConstraints *)malloc(sizeof(TNumberConstraints) * constraintsSmall->numberConstraints);
    for (i = 0; i < constraintsSmall->numberConstraints; i++)
    {
        szConst = 0;
        for (j = constraintsSmall->ElementsConstraints[i]; j < constraintsSmall->ElementsConstraints[i + 1]; j++)
        {
            el = constraintsSmall->Elements[j];
            if (constraintsSmall->xAsterisc[el] != precision)
            {
                szConst++;
            }
        }
        if (szConst >= 2)
        {
            nNewConstraints++;
            newCont += szConst;
            constraintsValided[i] = 1;
        }
        else
        {
            constraintsValided[i] = 0;
        }
    }
    cutSmall *newConstraintsSmall = AllocStrCutSmall(newCont, nNewConstraints, constraintsSmall->numberVariables);
    newConstraintsSmall->ElementsConstraints[0] = 0;

    newCont = 0;
    TRightSide rhsTemp;
    TNumberConstraints aux = 0;
    for (i = 0; i < constraintsSmall->numberConstraints; i++)
    {
        if (constraintsValided[i])
        {
            rhsTemp = constraintsSmall->rightSide[i];
            for (j = constraintsSmall->ElementsConstraints[i]; j < constraintsSmall->ElementsConstraints[i + 1]; j++)
            {
                el = constraintsSmall->Elements[j];
                if (constraintsSmall->xAsterisc[el] != precision)
                {
                    newConstraintsSmall->Coefficients[newCont] = constraintsSmall->Coefficients[j];
                    newConstraintsSmall->Elements[newCont] = constraintsSmall->Elements[j];
                    newCont++;
                }
                else
                {
                    rhsTemp -= constraintsSmall->Coefficients[j];
                }
            }
            newConstraintsSmall->ElementsConstraints[aux + 1] = newCont;
            newConstraintsSmall->rightSide[aux] = rhsTemp;
            newConstraintsSmall->originalConstraints[aux] = constraintsSmall->originalConstraints[i];
            aux++;
        }
    }
    for (i = 0; i < constraintsSmall->numberVariables; i++)
    {
        newConstraintsSmall->xAsterisc[i] = constraintsSmall->xAsterisc[i];
    }
    free(constraintsValided);
    return newConstraintsSmall;
}

cutFull *runCG1_mainCpu(cutFull *constraintsFull, int precision, double timeLimit, int limitCoefMAx, TNumberConstraints numberConstraintsInitial, TCont contInitial)
{
    cutFull *outCutFull;
    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull);
    cutSmall *constraintsSmall = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    free(intOrFloat);
    cutSmall *newConstraintsSmall = removeXAsteriscEqualOne(constraintsSmall, precision);
    free(constraintsSmall);
    int nRuns = newConstraintsSmall->numberConstraints;
    int i;
    solutionStructCG1 *h_solution_cg1 = allocStructSolCG1(nRuns);

    runCG_Phase_1_cpu(newConstraintsSmall, h_solution_cg1, precision, timeLimit, limitCoefMAx);

    int cont = 0;
    for (i = 0; i < newConstraintsSmall->numberConstraints; i++)
    {
        if (h_solution_cg1->sConstraints[i] != -1)
        {
            cont++;
        }
    }
    if (cont > 0)
    {
        printf("Number cuts generated in the phase 1: %d \n", cont);
        outCutFull = createCutsStrongPhaseOne(newConstraintsSmall, constraintsFull, h_solution_cg1, cont, precision, nRuns, numberConstraintsInitial, contInitial);
        free(h_solution_cg1);
        freeStrCutFull(constraintsFull);
        free(newConstraintsSmall);
        return outCutFull;

        //return constraintsFull;
    }
    else
    {
        printf("No cuts generate\n");
        free(newConstraintsSmall);
        free(h_solution_cg1);
        return constraintsFull;
    }
}

void runCG_Phase_1_cpu(cutSmall *constraintsSmall, solutionStructCG1 *solutionCG1, int precision, double timeconst, int limitCoefMax)
{
    TNumberConstraints *constraints;
    int violation = 0;
    int i, j, k;
    constraints = (TNumberConstraints *)malloc(sizeof(TNumberConstraints) * constraintsSmall->numberConstraints);
    for (i = 0; i < constraintsSmall->numberConstraints; i++)
    {
        constraints[i] = i;
    }
    for (k = 0; k < constraintsSmall->numberConstraints; k++)
    {
        int res = constraints[k];
        int sz_u = constraintsSmall->ElementsConstraints[res + 1] - constraintsSmall->ElementsConstraints[res] + 1;
        int *u = (int *)malloc(sizeof(int) * sz_u);
        int n1 = -1, d1 = -1, el, rhs, aux, value_tes;
        aux = 0;

        for (i = 0; i < sz_u; i++)
        {
            u[i] = 0;
        }
        for (j = constraintsSmall->ElementsConstraints[res]; j < constraintsSmall->ElementsConstraints[res + 1]; j++)
        {
            if (constraintsSmall->Coefficients[j] >= 2)
            {
                if (constraintsSmall->Coefficients[j] <= limitCoefMax)
                {
                    u[aux] = constraintsSmall->Coefficients[j];
                }
                else
                {
                    u[aux] = limitCoefMax;
                }
                aux++;
            }
            else
            {
                if (constraintsSmall->Coefficients[j] <= -2)
                {

                    u[aux] = constraintsSmall->Coefficients[j] * (-1);
                    if (u[aux] > limitCoefMax)
                    {
                        u[aux] = limitCoefMax;
                    }
                    aux++;
                }
                else
                {
                    sz_u--;
                }
            }
        }
        if (constraintsSmall->rightSide[res] >= 2)
        {
            u[sz_u - 1] = constraintsSmall->rightSide[res];
        }
        else
        {
            sz_u--;
        }
        value_tes = sz_u;
        for (i = 0; i < sz_u; i++)
        {
            for (j = i + 1; j < sz_u; j++)
            {
                if ((u[i] != 0) && (u[j] != 0) && (i != j))
                {
                    if (u[i] % u[j] == 0)
                    {
                        u[j] = 0;
                        value_tes--;
                    }
                    if (u[j] % u[i] == 0)
                    {
                        u[i] = 0;
                        value_tes--;
                    }
                }
            }
        }
        for (i = 0; i < sz_u; i++)
        {
            for (j = i + 1; j < sz_u; j++)
            {
                if ((u[i] < u[j]))
                {
                    aux = u[i];
                    u[i] = u[j];
                    u[j] = aux;
                }
            }
        }
        sz_u = value_tes;
        int nBest = -1, dBest = -1, violation_best = 0;
        for (j = 0; j < sz_u; j++)
        {
            d1 = u[j];
            n1 = 1;
            while (n1 < d1)
            {
                rhs = 0;
                violation = 0;
                value_tes = 0;
                rhs = constraintsSmall->rightSide[res] * n1;
                if (rhs % d1 != 0)
                {

                    for (i = constraintsSmall->ElementsConstraints[res]; i < constraintsSmall->ElementsConstraints[res + 1]; i++)
                    {
                        el = constraintsSmall->Elements[i];
                        aux = constraintsSmall->Coefficients[i] * n1;
                        if (((aux > 0 && d1 < 0) || (aux < 0 && d1 > 0)) && (aux % d1 != 0))
                        {
                            aux = (aux / d1) - 1;
                        }
                        else
                        {
                            aux = aux / d1;
                        }
                        value_tes += aux * constraintsSmall->xAsterisc[el];
                    }
                    if (((rhs > 0 && d1 < 0) || (rhs < 0 && d1 > 0)) && (rhs % d1 != 0))
                    {
                        rhs = (rhs / d1) - 1;
                    }
                    else
                    {
                        rhs = rhs / d1;
                    }

                    if (value_tes > rhs * precision)
                    {
                        violation = value_tes - (rhs * precision);
                        if (violation > violation_best)
                        {
                            violation_best = violation;
                            nBest = n1;
                            dBest = d1;
                        }
                    }
                }
                n1++;
            }
        }

        if (violation_best != 0)
        {
            solutionCG1->sConstraints[k] = res;
            solutionCG1->sNumerator[k] = nBest;
            solutionCG1->sDenominator[k] = dBest;
        }
        else
        {
            solutionCG1->sConstraints[k] = -1;
            solutionCG1->sNumerator[k] = -1;
            solutionCG1->sDenominator[k] = -1;
        }
        free(u);
    }

    free(constraints);
}

void quicksortTParameters(TParametersGroup *values, TNumberConstraints *idc, TNumberConstraints began, TNumberConstraints end)
{
    TNumberConstraints i, j;
    TParametersGroup pivo, aux;
    TNumberConstraints auxIdc;
    i = began;
    j = end - 1;
    pivo = values[(began + end) / 2];
    while (i <= j)
    {
        while (values[i] > pivo && i < end)
        {
            i++;
        }
        while (values[j] < pivo && j > began)
        {
            j--;
        }
        if (i <= j)
        {
            aux = values[i];
            values[i] = values[j];
            values[j] = aux;

            auxIdc = idc[i];
            idc[i] = idc[j];
            idc[j] = auxIdc;

            i++;
            j--;
        }
    }
    if (j > began)
        quicksortTParameters(values, idc, began, j + 1);
    if (i < end)
        quicksortTParameters(values, idc, i, end);
}

TNumberConstraints *returnVectorParametersConstraints(cutSmall *constraintsSmall, int precision, TNumberConstraints Constraint, TNumberConstraints numberMaxConst, TParametersGroup pCom, TParametersGroup pCoef, TParametersGroup pSignal, TParametersGroup pFrac)
{
    TNumberConstraints *groupConstraints = (TNumberConstraints *)malloc(sizeof(TNumberConstraints) * numberMaxConst);
    TNumberConstraints *countEquals = (TNumberConstraints *)malloc(sizeof(TNumberConstraints) * constraintsSmall->numberConstraints);
    TCoefficients *countCoef = (TCoefficients *)malloc(sizeof(TCoefficients) * constraintsSmall->numberConstraints);
    TNumberConstraints *countSignalDif = (TNumberConstraints *)malloc(sizeof(TNumberConstraints) * constraintsSmall->numberConstraints);
    TXAsterisc *countFrac = (TXAsterisc *)malloc(sizeof(TXAsterisc) * constraintsSmall->numberConstraints);
    TXAsterisc *FractinalSolution = (TXAsterisc *)malloc(sizeof(TXAsterisc) * constraintsSmall->numberVariables);
    TParametersGroup *total = (TParametersGroup *)malloc(sizeof(TParametersGroup) * constraintsSmall->numberConstraints);
    TNumberConstraints *idc = (TNumberConstraints *)malloc(sizeof(TNumberConstraints) * constraintsSmall->numberConstraints);

    int i, j, k, z;

    for (i = 0; i < constraintsSmall->numberVariables; i++)
    {
        if (constraintsSmall->xAsterisc[i] % precision != 0)
        {
            if (constraintsSmall->xAsterisc[i] < 0)
            {
                k = constraintsSmall->xAsterisc[i] / precision;
                FractinalSolution[i] = -1 * (k)*precision - constraintsSmall->xAsterisc[i];
            }
            else
            {
                k = constraintsSmall->xAsterisc[i] / precision;
                FractinalSolution[i] = (k + 1) * precision - constraintsSmall->xAsterisc[i];
            }
        }
        else
        {
            FractinalSolution[i] = 0;
        }
    }

    for (i = 0; i < constraintsSmall->numberConstraints; i++)
    {
        countEquals[i] = 0;
        countCoef[i] = 0;
        countSignalDif[i] = 0;
        countFrac[i] = 0;
        idc[i] = i;
    }
    for (z = constraintsSmall->ElementsConstraints[Constraint]; z < constraintsSmall->ElementsConstraints[Constraint + 1]; z++)
    {
        for (i = 0; i < constraintsSmall->numberConstraints; i++)
        {
            if (i != Constraint)
            {
                for (j = constraintsSmall->ElementsConstraints[i]; j < constraintsSmall->ElementsConstraints[i + 1]; j++)
                {
                    if (constraintsSmall->Elements[z] == constraintsSmall->Elements[j])
                    {
                        countEquals[i]++;
                        countCoef[i] += constraintsSmall->Coefficients[j];
                        if (constraintsSmall->Coefficients[j] * constraintsSmall->Coefficients[z] < 0)
                        {
                            countSignalDif[i]++;
                        }
                        TNumberConstraints el = constraintsSmall->Elements[z];
                        countFrac[i] += FractinalSolution[el];
                    }
                }
            }
        }
    }

    for (i = 0; i < constraintsSmall->numberConstraints; i++)
    {
        total[i] = pCom * (TParametersGroup)countEquals[i] + pCoef * (TParametersGroup)countCoef[i] + pSignal * (TParametersGroup)countSignalDif[i] + pFrac * (TParametersGroup)countFrac[i];
    }
    free(countCoef);
    free(countEquals);
    free(countSignalDif);
    free(countFrac);
    quicksortTParameters(total, idc, 0, constraintsSmall->numberConstraints);
    j = 0;
    TNumberConstraints posMaior = 0;
    if (posMaior + numberMaxConst < constraintsSmall->numberConstraints - 1)
    {
        for (i = posMaior; i < posMaior + numberMaxConst; i++)
        {
            groupConstraints[j] = idc[i];
            j++;
        }
    }
    else
    {
        for (i = constraintsSmall->numberConstraints - (numberMaxConst + 1); i < constraintsSmall->numberConstraints - 1; i++)
        {
            groupConstraints[j] = idc[i];
            j++;
        }
    }
    free(idc);
    free(total);
    free(FractinalSolution);
    return groupConstraints;
}

void shuffle_Set(TNumberConstraints *vec, TNumberConstraints numberMaxConstraints, TNumberConstraints n)
{
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    int i, j, aux;
    float *num_temp = (float *)malloc(sizeof(float) * numberMaxConstraints);                                       //free
    TNumberConstraints *vec_aux = (TNumberConstraints *)malloc(sizeof(TNumberConstraints) * numberMaxConstraints); //free
    aux = n / numberMaxConstraints;
    for (i = 0; i < aux; i++)
    {
        for (j = 0; j < numberMaxConstraints; j++)
        {

            num_temp[j] = (float)(rand() % RAND_MAX);
            vec_aux[j] = vec[i * numberMaxConstraints + j];
        }
        quicksortTParameters(num_temp, vec_aux, 0, numberMaxConstraints);
        //bubble_sort(num_temp,vec_aux,numberMaxConstraints);
        for (j = 0; j < numberMaxConstraints; j++)
        {
            vec[i * numberMaxConstraints + j] = vec_aux[j];
        }
    }
    free(num_temp);
    free(vec_aux);
}

cutFull *runCG2_mainCpu(cutFull *constraintsFull, int precision, double timeLimit, TNumberConstraints numberMaxConst, int nRuns, int maxDenominator, TSPosition posMaior, float pCom, float pCoef, float pSignal, float pFrac)
{
    double startT = omp_get_wtime();
    double _time = 0;
    cutFull *outCutFull;
    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull);
    cutSmall *constraintsSmall = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    free(intOrFloat);
    cutSmall *newConstraintsSmall = removeXAsteriscEqualOne(constraintsSmall, precision);
    free(constraintsSmall);
    int i, j;
    _time = ((double)timeLimit - (omp_get_wtime() - startT));
    if (_time < 1)
    {
        free(newConstraintsSmall);
        return constraintsFull;
    }

    solutionStructCG2 *solutionCG2 = allocStructSolCG2(numberMaxConst, nRuns);

    _time = ((double)timeLimit - (omp_get_wtime() - startT));
    if (_time < 1)
    {
        free(solutionCG2);
        free(newConstraintsSmall);
        return constraintsFull;
    }

    TNumberConstraints *setConstraint = (TNumberConstraints *)malloc(sizeof(TNumberConstraints) * numberMaxConst * nRuns);
    TNumberConstraints *auxSet = (TNumberConstraints *)malloc(sizeof(TNumberConstraints) * numberMaxConst);

    for (i = 0; i < nRuns; i++)
    {
        int p = i % newConstraintsSmall->numberConstraints;
        auxSet = returnVectorParametersConstraints(newConstraintsSmall, precision, p, numberMaxConst, pCom, pCoef, pSignal, pFrac);
        for (j = 0; j < numberMaxConst; j++)
        {
            setConstraint[j + i * numberMaxConst] = auxSet[j];
        }
        free(auxSet);
    }

    shuffle_Set(setConstraint, numberMaxConst, numberMaxConst * nRuns);
    _time = ((double)timeLimit - (omp_get_wtime() - startT));
    if (_time < 1)
    {
        free(newConstraintsSmall);
        free(solutionCG2);
        free(setConstraint);
        return constraintsFull;
    }

    runCG_phase_2_cpu(newConstraintsSmall, solutionCG2, precision, _time, setConstraint, numberMaxConst, nRuns, maxDenominator);

    int cont = 0;
    _time = ((double)timeLimit - (omp_get_wtime() - startT));
    if (_time < 1)
    {
        free(newConstraintsSmall);
        free(solutionCG2);
        free(setConstraint);
        return constraintsFull;
    }
    int k;
    for (k = 0; k < nRuns; k++)
    {
        if (solutionCG2->sConstraints[0 + k * numberMaxConst] != -1)
        {
            cont++;
        }
    }
    if (cont > 0)
    {
        printf("Number of Cuts in the second phase in CPU:%d\n", cont);
        int nT = nRuns;
        int nB = 1;
        outCutFull = createCutsStrongPhaseTwo(newConstraintsSmall, constraintsFull, solutionCG2, numberMaxConst, cont, precision, nT, nB);
    }
    else
    {
        //printf("Nothing  Cuts in phase 2\n");
        free(newConstraintsSmall);
        free(setConstraint);
        free(solutionCG2);

        return constraintsFull;
    }

    free(setConstraint);
    free(solutionCG2);
    free(newConstraintsSmall);
    freeStrCutFull(constraintsFull);
    return outCutFull;
}

void runCG_phase_2_cpu(cutSmall *constraintsSmall, solutionStructCG2 *solutionCG2, int precision, double timeLimit, TNumberConstraints *setConstraint, TNumberConstraints numberMaxConst, int nRuns, int maxDenominator)
{

    int k;
    int qtMp = 20;
    for (k = 0; k < nRuns; k++)
    {
        //printf("k:%d\n",k);
        int mult_1, mult_2, i, j;
        TElements el;
        TRightSide rhs1, rhs2;
        TXAsterisc value_tes, violation = 0;
        TCoefficients aux;
        TSNumerator n1_best = -1, n2_best = -1;
        TSDenominator d1_best = -1, d2_best = -1;
        TSPosition position;
        TNumberConstraints rest_a, rest_b;
        // double timeCurrent = omp_get_wtime();
        TSNumerator Numerator[qtMp];
        TSDenominator Denominator[qtMp];
        TCoefficients *Coef = (TCoefficients *)malloc(sizeof(TCoefficients) * (constraintsSmall->numberVariables));
        TCoefficients *Coef2 = (TCoefficients *)malloc(sizeof(TCoefficients) * (constraintsSmall->numberVariables));
        for (i = 0; i < qtMp; i++)
        {
            Denominator[i] = rand() % maxDenominator + 2;
            Numerator[i] = rand() % (Denominator[i] - 1);
        }

        for (mult_1 = 0; mult_1 < qtMp; mult_1++)
        {
            memset(Coef, 0, sizeof(TCoefficients) * constraintsSmall->numberVariables);
            rhs1 = 0;
            for (rest_a = 0; rest_a < numberMaxConst; rest_a++)
            {
                for (i = constraintsSmall->ElementsConstraints[setConstraint[k * numberMaxConst + rest_a]]; i < constraintsSmall->ElementsConstraints[setConstraint[k * numberMaxConst + rest_a] + 1]; i++)
                {

                    el = constraintsSmall->Elements[i];
                    Coef[el] += constraintsSmall->Coefficients[i] * Numerator[mult_1];
                }
                rhs1 += constraintsSmall->rightSide[setConstraint[k * numberMaxConst + rest_a]] * Numerator[mult_1];
                for (mult_2 = 0; mult_2 < qtMp; mult_2++)
                {
                    memset(Coef2, 0, sizeof(TCoefficients) * constraintsSmall->numberVariables);
                    value_tes = 0;
                    rhs2 = 0;
                    for (rest_b = rest_a + 1; rest_b < numberMaxConst; rest_b++)
                    {
                        for (j = constraintsSmall->ElementsConstraints[setConstraint[k * numberMaxConst + rest_b]]; j < constraintsSmall->ElementsConstraints[setConstraint[k * numberMaxConst + rest_b] + 1]; j++)
                        {
                            el = constraintsSmall->Elements[j];
                            Coef2[el] += constraintsSmall->Coefficients[j] * Numerator[mult_2];
                        }
                        rhs2 += constraintsSmall->rightSide[setConstraint[k * numberMaxConst + rest_b]] * Numerator[mult_2];
                    }
                    for (j = 0; j < constraintsSmall->numberVariables; j++)
                    {
                        if ((Coef[j] * Denominator[mult_2] + Coef2[j] * Denominator[mult_1]) / (Denominator[mult_1] * Denominator[mult_2]) < 0)
                        {
                            aux = (Coef[j] * Denominator[mult_2] + Coef2[j] * Denominator[mult_1]) / (Denominator[mult_1] * Denominator[mult_2]) - 1;
                        }
                        else
                        {
                            aux = (Coef[j] * Denominator[mult_2] + Coef2[j] * Denominator[mult_1]) / (Denominator[mult_1] * Denominator[mult_2]);
                        }
                        value_tes += aux * constraintsSmall->xAsterisc[j];
                    }

                    if ((rhs1 * Denominator[mult_2] + rhs2 * Denominator[mult_1]) / (Denominator[mult_1] * Denominator[mult_2]) < 0)
                    {
                        aux = (rhs1 * Denominator[mult_2] + rhs2 * Denominator[mult_1]) / (Denominator[mult_1] * Denominator[mult_2]) - 1;
                    }
                    else
                    {
                        aux = (rhs1 * Denominator[mult_2] + rhs2 * Denominator[mult_1]) / (Denominator[mult_1] * Denominator[mult_2]);
                    }

                    if ((value_tes > aux * precision) && (value_tes - (aux * precision) > violation))
                    {
                        violation = value_tes - (aux * precision);
                        //printf("violation %d\n",violation);
                        n1_best = Numerator[mult_1];
                        d1_best = Denominator[mult_1];
                        n2_best = Numerator[mult_2];
                        d2_best = Denominator[mult_2];
                        position = rest_a;
                    }
                }
            }
        }

        if (violation > 0)
        {
            for (i = 0; i < numberMaxConst; i++)
            {
                solutionCG2->sConstraints[i + k * numberMaxConst] = setConstraint[i + k * numberMaxConst]; //CPU ja vai ter
            }

            solutionCG2->sPosition[k] = position;
            solutionCG2->sNumerator1[k] = n1_best;
            solutionCG2->sDenominator1[k] = d1_best;
            solutionCG2->sNumerator2[k] = n2_best;
            solutionCG2->sDenominator2[k] = d2_best;
        }
        else
        {

            for (i = 0; i < numberMaxConst; i++)
            {
                solutionCG2->sConstraints[i + k * numberMaxConst] = -1;
            }
            solutionCG2->sPosition[k] = 0;
            solutionCG2->sNumerator1[k] = -1;
            solutionCG2->sDenominator1[k] = -1;
            solutionCG2->sNumerator2[k] = -1;
            solutionCG2->sDenominator2[k] = -1;
        }

        free(Coef);
        free(Coef2);
    }
}

void createSolutionsInitial(int *solution, int sz)
{
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    int i;
    for (i = 0; i < sz; i++)
    {
        solution[i] = rand() % 2;
    }
}

void createInitialCoverGRASP(int *solution, int sz, cutSmall *constraintsSmall, int precision, int constraint, float alpha)
{
    int verifyCompleteSolution = 0;
    int *setId = (int *)malloc(sz * sizeof(int));
    double c_min = 0, c_max = 0, value, value_best;
    int i, aux = 0, el, szAux , lhs;
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));


    for (i = 0; i < sz; i++)
    {
        solution[i] = 0;
    }
    for (i = constraintsSmall->ElementsConstraints[constraint]; i < constraintsSmall->ElementsConstraints[constraint + 1]; i++)
    {
        setId[aux] = i;
        aux++;
    }
    szAux = aux;
    lhs = 0;
    int test =0;
    while (verifyCompleteSolution == 0)
    {   
        aux = 0;
        //printf("coef: ");
        for (i = 0; i < szAux; i++)
        {
            //printf("%d\t", setId[i]);
            el = constraintsSmall->Elements[setId[i]];
            value = (double)constraintsSmall->xAsterisc[el] / (double)constraintsSmall->Coefficients[setId[i]];
            if (aux == 0)
            {
                c_min = value;
                c_max = value;
                aux = 1;
            }
            if (value < c_min)
            {
                c_min = value;
            }
            if (value > c_max)
            {
                c_max = value;
            }
        }

        value_best = c_max - alpha*(c_max - c_min);
        //printf("\nc_max %f c_min %f alpha %f\n", c_max, c_min, alpha);
        //printf("%d-Value Best: %f\n",test, value_best);
        test++;
        int *setTemp  = (int*)malloc(szAux*sizeof(int));
        int *posTemp = (int*)malloc(szAux*sizeof(int));
        int aux_t = 0;
        for(i = 0 ; i<szAux;i++){
            el = constraintsSmall->Elements[setId[i]];
            value = (double)constraintsSmall->xAsterisc[el] / (double)constraintsSmall->Coefficients[setId[i]];
            if(value >= value_best){
                setTemp[aux_t] = setId[i];
                posTemp[aux_t] = i;
                aux_t++;
            }
        }
        int itemAdd = rand() % aux_t;
        el = setTemp[itemAdd] - constraintsSmall->ElementsConstraints[constraint];
        solution[el] = 1;
        szAux--;
        lhs += constraintsSmall->Coefficients[ setTemp[itemAdd] ];
        if(lhs -1e-5 > constraintsSmall->rightSide[constraint]){
            verifyCompleteSolution = 1;
            free(setTemp);
            free(posTemp);
            break;
        }
        free(setTemp);
        int *newSet = (int*)malloc(sizeof(int)*szAux);
        aux_t = 0;
        for(i=0;i<szAux+1;i++){
            if( i != posTemp[itemAdd]){
                newSet[aux_t] = setId[i];
                aux_t++;
            }
        }
        printf("Testes: %d %d\n", aux_t, szAux);
        free(setId);
        setId  = newSet; 
        free(posTemp);
    }
    free(setId);

}

double calcFO(int *solution, cutSmall *constraintsSmall, int precision, TNumberConstraints constraint)
{
    double fo;
    double lhs = 0;
    int aux = 0;
    int i, el;
    for (i = constraintsSmall->ElementsConstraints[constraint]; i < constraintsSmall->ElementsConstraints[constraint + 1]; i++)
    {
        el = constraintsSmall->Elements[i];
        lhs += (double)(constraintsSmall->Coefficients[i] * constraintsSmall->xAsterisc[el] * solution[aux]);
        aux++;
    }
    lhs = lhs / precision;
    fo = lhs - constraintsSmall->rightSide[constraint];
    return fo;
}

double calcViolation(cutSmall *newConstraintsSmall, TNumberConstraints constraints, int precision)
{
    double violation = 0.0;
    int i, el;
    for (i = newConstraintsSmall->ElementsConstraints[constraints]; i < newConstraintsSmall->ElementsConstraints[i + 1]; i++)
    {
        el = newConstraintsSmall->ElementsConstraints[i];
        violation += (double)newConstraintsSmall->Coefficients[i] * (double)newConstraintsSmall->xAsterisc[el];
    }
    violation /= (double)precision;
    violation -= (double)newConstraintsSmall->rightSide[constraints];
    return violation;
}

int verifySolutionCover(int *solution, cutSmall *constraintsSmall, int precision, TNumberConstraints constraint)
{
    double lhs = 0;
    int i, aux = 0;
    for (i = constraintsSmall->ElementsConstraints[constraint]; i < constraintsSmall->ElementsConstraints[constraint + 1]; i++)
    {
        //el = constraintsSmall->Elements[i];
        lhs += constraintsSmall->Coefficients[i] * solution[aux];
        aux++;
    }
    if ((lhs - 1e-5) > constraintsSmall->rightSide[constraint])
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void shuffleVectorInt(int *vec, int sz)
{
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    int j;
    float *num_temp = (float *)malloc(sizeof(float) * sz); //free
    for (j = 0; j < sz; j++)
    {
        num_temp[j] = (float)(rand() % RAND_MAX);
    }
    quicksortTParameters(num_temp, vec, 0, sz);
    free(num_temp);
}

void SortVectorGreedy(int *vec, int sz, cutSmall *constraintsSmall, TNumberConstraints constraint, int precision)
{
    float *num_temp = (float *)malloc(sizeof(float) * sz);
    int i, el, aux;
    for (i = 0; i < sz; i++)
    {
        el = constraintsSmall->ElementsConstraints[constraint] + vec[i];
        aux = constraintsSmall->Elements[el];
        num_temp[i] = (float)precision * (((float)constraintsSmall->xAsterisc[aux]) / ((float)constraintsSmall->Coefficients[el]));
    }
    quicksortTParameters(num_temp, vec, 0, sz);
    free(num_temp);
}

int *localSearch(int *solution, cutSmall *constraintsSmall, int precision, TNumberConstraints constraint)
{
    int i, aux1 = 0, aux2 = 0;
    int sz = constraintsSmall->ElementsConstraints[constraint + 1] - constraintsSmall->ElementsConstraints[constraint];
    int *solutionAux = (int *)malloc(sizeof(int) * sz);
    int *vAtv = (int *)malloc(sizeof(int) * sz);
    int *vNAtv = (int *)malloc(sizeof(int) * sz);
    for (i = 0; i < sz; i++)
    {
        solutionAux[i] = solution[i];
        if (solution[i] == 1)
        {
            vAtv[aux1] = i;
            aux1++;
        }
        else
        {
            vNAtv[aux2] = i;
            aux2++;
        }
    }
    shuffleVectorInt(vAtv, aux1);
    if (aux2 > 0)
    {
        SortVectorGreedy(vNAtv, aux2, constraintsSmall, constraint, precision);
    }
    solutionAux[vAtv[0]] = 0;
    i = 0;
    //int x = 0;
    do
    {
        if (i == aux2)
        {
            free(solutionAux);
            free(vAtv);
            free(vNAtv);
            return solution;
        }
        solutionAux[vNAtv[i]] = 1;
        i++;

    } while (verifySolutionCover(solutionAux, constraintsSmall, precision, constraint) == 0);

    // do
    // {
    //     printf("aqui: %d\n", x);
    //     x++;
    //     if (aux2 > 0)
    //     {
    //         solutionAux[vNAtv[i]] = 1;
    //         i++;
    //         aux2--;
    //     }
    //     else
    //     {
    //         solutionAux[vAtv[j]] = 0;
    //         j++;
    //     }
    //     if ((i == sz) || (j == aux1))
    //     {
    //         free(solutionAux);
    //         free(vAtv);
    //         free(vNAtv);
    //         return solution;
    //     }
    // } while (verifySolutionCover(solutionAux, constraintsSmall, precision, constraint) == 0);

    free(solution);
    free(vAtv);
    free(vNAtv);
    return solutionAux;
}

void copyAndVerifyPoolSolution(int *solution, int sz, int *poolSolution, int szPoolCutsMax, int *numberSolutionAtual)
{
    int i, j;
    int flag = *numberSolutionAtual;
    for (i = 0; i < *numberSolutionAtual; i++)
    {
        for (j = 0; j < sz; j++)
        {
            if (solution[j] != poolSolution[j + i * sz])
            {
                flag--;
                break;
            }
        }
    }
    if (flag == 0)
    {
        for (j = 0; j < sz; j++)
        {
            poolSolution[j + *numberSolutionAtual * sz] = solution[j];
        }
        //printf("Aqui deu certo\n");
        (*numberSolutionAtual)++;
    }
    //else
    //{
    //    printf("Cover repetido!");
    //}
}

int *createCoverGraspIndividual(cutSmall *constraintsSmall, int precision, TNumberConstraints constraint, int numberIteration, int szPoolCutsMax)
{

    int sz = constraintsSmall->ElementsConstraints[constraint + 1] - constraintsSmall->ElementsConstraints[constraint];
    int *solution = (int *)malloc(sizeof(int) * sz);
    int *poolSolution = (int *)malloc(sizeof(int) * sz * szPoolCutsMax);
    int nPoolSolution = 0;
    int ite = 0;
    for (ite = 0; ite < sz * szPoolCutsMax; ite++)
    {
        poolSolution[ite] = 0;
    }

    createInitialCoverGRASP(solution, sz, constraintsSmall, precision, constraint, 1);
    ite = 0;
    while (ite < numberIteration)
    {
        solution = localSearch(solution, constraintsSmall, precision, constraint);
        copyAndVerifyPoolSolution(solution, sz, poolSolution, szPoolCutsMax, &nPoolSolution);
        if (nPoolSolution == szPoolCutsMax)
        {
            break;
        }
        ite++;
    }
    free(solution);
    return poolSolution;
}

cutFull *runCC_mainCPu(cutFull *constraintsFull, int precision, int szCoverThread)
{

    //Cover_gpu *h_cover_new = AllocationStructCover(h_cover->cont*qnt_Cover_per_Thread,h_cover->numberConstraints*qnt_Cover_per_Thread);
    int i, j, k, w;
    double b = 0.0, a_barra = 0.0;
    int qnt;

    int *fillBag;
    // int nConstraintsInitial = h_cut->numberConstrains;
    cutCover *cutsCover;

    //int contInitial;// = cutsCover->cont;
    for (j = 0; j < szCoverThread; j++)
    {
        int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull);
        cutSmall *newConstraintsSmall = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
        //cutSmall *newConstraintsSmall = removeXAsteriscEqualOne(constraintsSmall, precision);
        int *solution = (int *)malloc(sizeof(int) * newConstraintsSmall->cont);
        createSolutionsInitial(solution, newConstraintsSmall->cont);

        // for(i=0;i<newConstraintsSmall->numberConstraints;i++){
        //     for(j=newConstraintsSmall->ElementsConstraints[i];j<newConstraintsSmall->ElementsConstraints[i+1];j++){
        //         int el = newConstraintsSmall->Elements[j];
        //         printf("%d x%d ", newConstraintsSmall->Coefficients[j],el);
        //     }
        //     printf("<= %d\n", newConstraintsSmall->rightSide[i]);
        //     getchar();
        // }

        cutsCover = CopyCutToCover(newConstraintsSmall);

        fillBag = calcCover(cutsCover, solution, 0);

        //contInitial = newConstraintsSmall->cont;
        for (i = 0; i < cutsCover->numberConstraints; i++)
        {
            qnt = 0;
            int el;
            for (k = cutsCover->ElementsConstraints[i]; k < cutsCover->ElementsConstraints[i + 1]; k++)
            {
                qnt += solution[k];
            }
            if ((fillBag[i] <= cutsCover->rightSide[i]) || (qnt <= 1))
            {
                continue;
            }

            int *n_coef = (int *)malloc(sizeof(int) * qnt);
            int *n_el = (int *)malloc(sizeof(int) * qnt);
            qnt = 0;
            for (k = cutsCover->ElementsConstraints[i]; k < cutsCover->ElementsConstraints[i + 1]; k++)
            {
                if (solution[k] == 1)
                {
                    n_coef[qnt] = cutsCover->Coefficients[k];
                    n_el[qnt] = k;
                    qnt++;
                }
            }
            b = (double)cutsCover->rightSide[i];
            double delta = 0;
            double phi = (double)(fillBag[i] - cutsCover->rightSide[i]);
            k = 1;
            a_barra = (double)n_coef[0];
            for (w = 1; w < qnt; w++)
            {

                delta = a_barra - (double)n_coef[w];
                if (((double)k * delta) < phi)
                {
                    a_barra = (double)n_coef[w];
                    phi = phi - (double)k * delta;
                }
                else
                {
                    a_barra = a_barra - (phi / (double)k);
                    phi = 0;
                    break;
                }
                k++;
            }
            if (phi > 0)
            {
                a_barra = (double)b / (double)qnt;
            }
            int *c_menus = (int *)malloc(sizeof(int) * qnt);
            int *c_mais = (int *)malloc(sizeof(int) * qnt);
            double *S_barra = (double *)malloc(sizeof(double) * (qnt + 1));
            int id1 = 0, id2 = 0, id3 = 0;
            for (w = 0; w < qnt; w++)
            {
                if ((double)n_coef[w] <= a_barra)
                {
                    c_menus[id1] = w;
                    id1++;
                }
                else
                {
                    c_mais[id2] = w;
                    id2++;
                }
            }
            S_barra[id3] = 0;
            id3++;
            for (w = 0; w < id2; w++)
            {
                S_barra[id3] = S_barra[id3 - 1] + a_barra;
                id3++;
            }
            for (w = 0; w < id1; w++)
            {
                S_barra[id3] = S_barra[id3 - 1] + (double)n_coef[c_menus[w]];
                id3++;
            }
            int ini = 0;
            for (w = cutsCover->ElementsConstraints[i]; w < cutsCover->ElementsConstraints[i + 1]; w++)
            {
                for (ini = 0; ini < id3 - 1; ini++)
                {
                    if ((cutsCover->Coefficients[w] > S_barra[ini]) && (cutsCover->Coefficients[w] <= S_barra[ini + 1]))
                    {
                        cutsCover->Coefficients[w] = ini;
                        break;
                    }
                    if (ini + 1 == id3 - 1)
                        cutsCover->Coefficients[w] = ini + 1;
                }
            }
            for (w = 0; w < id1; w++)
            {
                el = n_el[c_menus[w]];
                cutsCover->Coefficients[el] = 1;
            }

            cutsCover->rightSide[i] = qnt - 1;
            free(c_menus);
            free(c_mais);
            free(S_barra);
            free(n_coef);
            free(n_el);
        }
        int qnt_cuts_cover = 0;
        int *idc_cover = (int *)malloc(sizeof(int) * cutsCover->numberConstraints);
        for (i = 0; i < cutsCover->numberConstraints; i++)
        {
            idc_cover[i] = 0;
            if (cutsCover->rightSide[i] != newConstraintsSmall->rightSide[i])
            {
                idc_cover[i] = 1;
                qnt_cuts_cover++;
            }
        }
        constraintsFull = createCutsCover(newConstraintsSmall, constraintsFull, cutsCover, idc_cover, qnt_cuts_cover);
        printf("Number cuts cover: %d\n", qnt_cuts_cover);
        free(cutsCover);
        free(newConstraintsSmall);
        //free(constraintsSmall);
        free(intOrFloat);
        free(idc_cover);
        free(solution);
        free(fillBag);
    }

    return constraintsFull;
}

int *creatFillBagGrasp(int *poolSolution, cutSmall *constraintsSmall, int constraint, int szPoolCutsMax, int szConstraint)
{
    int *fillBag = (int *)malloc(sizeof(int) * szPoolCutsMax);
    int i, j, aux;
    for (i = 0; i < szPoolCutsMax; i++)
    {
        fillBag[i] = 0;
        aux = 0;
        for (j = constraintsSmall->ElementsConstraints[constraint]; j < constraintsSmall->ElementsConstraints[constraint + 1]; j++)
        {
            fillBag[i] += constraintsSmall->Coefficients[j] * poolSolution[aux + i * szConstraint];
            aux++;
        }
    }
    return fillBag;
}

#ifdef DEBUG

cutFull *runCC_mainCPuDebug(cutFull *constraintsFull, int precision, char **nameConstraints, char **nameVariables, double *sol, int szPoolCutsMax, int nIterationGrasp)
{
    printf("CC with debug\n");
    //Cover_gpu *h_cover_new = AllocationStructCover(h_cover->cont*qnt_Cover_per_Thread,h_cover->numberConstraints*qnt_Cover_per_Thread);
    int i, j, k, w, itePool = 0;
    double b = 0.0, a_barra = 0.0;
    int qnt;
    int c_XSolution = 0;
    int c_AuxSolution = 0;

    // int nConstraintsInitial = h_cut->numberConstrains;
    cutCover *cutsCover;

    //int contInitial;// = cutsCover->cont;

    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull);
    cutSmall *newConstraintsSmall = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    cutsCover = CopyCutToCover(newConstraintsSmall);

    //contInitial = newConstraintsSmall->cont;
    for (i = 0; i < cutsCover->numberConstraints; i++)
    {
        if (cutsCover->rightSide[i] <= 1)
        {
            continue;
        }
        c_AuxSolution = 0;
        c_XSolution = 0;
        int szConstraint = cutsCover->ElementsConstraints[i + 1] - cutsCover->ElementsConstraints[i];
        cutCover *cutsCoverSolution = AllocStrCover(szConstraint * szPoolCutsMax, szPoolCutsMax);
        int *poolSolution = createCoverGraspIndividual(newConstraintsSmall, precision, i, nIterationGrasp, szPoolCutsMax);
        // for (k = 0; k < szPoolCutsMax; k++)
        // {
        //     for (w = 0; w < szConstraint; w++)
        //     {
        //         printf(" %d %d ", poolSolution[w + k * szConstraint], cutsCover->Coefficients[w + cutsCover->ElementsConstraints[i]]);
        //     }
        //     printf("%d\n", cutsCover->rightSide[i]);
        // }
        int *fillBag = creatFillBagGrasp(poolSolution, newConstraintsSmall, i, szPoolCutsMax, szConstraint);
        cutsCoverSolution->ElementsConstraints[0] = 0;
        for (itePool = 0; itePool < szPoolCutsMax; itePool++)
        {

            qnt = 0;
            int el;
            int caux = 0;
            for (k = cutsCover->ElementsConstraints[i]; k < cutsCover->ElementsConstraints[i + 1]; k++)
            {
                // printf("%d\t ", cutsCover->Coefficients[k]);
                qnt += poolSolution[caux + itePool * szConstraint];
                caux++;
            }
            // printf("rhs:  %d\n",cutsCover->rightSide[i]);
            // getchar();
            //int nC = cutsCover->ElementsConstraints[i + 1] - cutsCover->ElementsConstraints[i];
            // printf("tamanho de C: %d\n Tamanho de N: %d\n",qnt, nC);
            //printf("fill Bag: %d \t rhs %d \t qnt %d \n", fillBag[itePool], cutsCover->rightSide[i], qnt);
            //getchar();
            if (fillBag[itePool] <= cutsCover->rightSide[i]) //|| (qnt <= 0) || (qnt == szConstraint))
            {
                continue;
            }
            //printf("fill Bag: %d \t rhs %d \t qnt %d \n", fillBag[itePool], cutsCover->rightSide[i], qnt);
            int *n_coef = (int *)malloc(sizeof(int) * qnt);
            int *n_el = (int *)malloc(sizeof(int) * qnt);
            qnt = 0;
            caux = 0;
            for (k = cutsCover->ElementsConstraints[i]; k < cutsCover->ElementsConstraints[i + 1]; k++)
            {
                if (poolSolution[caux + itePool * szConstraint] == 1)
                {
                    n_coef[qnt] = cutsCover->Coefficients[k];
                    n_el[qnt] = k;
                    qnt++;
                }
                caux++;
            }
            //printf("qnt posterior: %d\n", qnt);
            b = (double)cutsCover->rightSide[i];
            double delta = 0;
            double phi = ((double)fillBag[itePool] - (double)cutsCover->rightSide[i]);
            k = 1;
            a_barra = (double)n_coef[0];
            for (w = 1; w < qnt; w++)
            {

                delta = a_barra - (double)n_coef[w];
                if ((double)k * delta < phi)
                {
                    a_barra = (double)n_coef[w];
                    phi = phi - (double)k * delta;
                }
                else
                {
                    a_barra = a_barra - (phi / (double)k);
                    phi = 0;
                    break;
                }
                k++;
            }
            if (phi > 0)
            {
                a_barra = (double)b / (double)qnt;
            }

            int *c_menus = (int *)malloc(sizeof(int) * qnt);
            int *c_mais = (int *)malloc(sizeof(int) * qnt);
            double *S_barra = (double *)malloc(sizeof(double) * (qnt + 1));
            double *aux = (double *)malloc(sizeof(double) * qnt);
            int id1 = 0, id2 = 0; //, id3 = 0;

            for (w = 0; w < qnt; w++)
            {
                if ((double)n_coef[w] <= a_barra)
                {
                    c_menus[id1] = w;
                    id1++;
                }
                else
                {
                    c_mais[id2] = w;
                    id2++;
                }
            }
            // S_barra[id3] = 0;
            // id3++;

            for (w = 0; w < qnt; w++)
            {
                if (n_coef[w] < a_barra)
                {
                    aux[w] = n_coef[w];
                }
                else
                {
                    aux[w] = a_barra;
                }
            }
            quicksortDouble(aux, 0, qnt);
            S_barra[0] = 0;
            for (w = 1; w < qnt; w++)
            {
                S_barra[w] = S_barra[w - 1] + aux[w];
            }
            S_barra[qnt] = cutsCover->rightSide[i];

            // for(w = 0; w<id2; w++)
            // {
            //     S_barra[id3] = S_barra[id3-1] + a_barra;
            //     id3++;
            // }
            // for(w = 0; w<id1; w++)
            // {
            //     S_barra[id3] = S_barra[id3-1] + (double)n_coef[ c_menus[w] ];
            //     id3++;
            // }

            // for(w=0;w<=qnt;w++){
            //     printf("S(%d) = %lf\n", w, S_barra[w]);
            // }
            // printf("a_barra: %lf\n", a_barra);
            // getchar();

            int ini = 0;

            for (w = cutsCover->ElementsConstraints[i]; w < cutsCover->ElementsConstraints[i + 1]; w++)
            {
                int flag = 0;
                //cutsCoverSolution->Coefficients[c_XSolution] = cutsCover->Coefficients[w];
                for (ini = 0; ini < qnt; ini++)
                {
                    if ((cutsCover->Coefficients[w] > S_barra[ini]) && (cutsCover->Coefficients[w] <= S_barra[ini + 1]))
                    {
                        cutsCoverSolution->Coefficients[c_XSolution] = ini;
                        flag = 1;
                        break;
                    }
                }
                if (flag == 0)
                {
                    printf("Aqui!");
                    getchar();
                    cutsCoverSolution->Coefficients[c_XSolution] = cutsCover->Coefficients[w];
                }
                c_XSolution++;
            }

            cutsCoverSolution->ElementsConstraints[c_AuxSolution + 1] = c_XSolution;

            for (w = 0; w < id1; w++)
            {
                el = n_el[c_menus[w]];
                caux = el - cutsCover->ElementsConstraints[i];
                cutsCoverSolution->Coefficients[caux + itePool * szConstraint] = 1;
            }
            //printf("Aqui ta quant %d = %d = %d = %d\n",i, qnt, cutsCover->rightSide[i], newConstraintsSmall->rightSide[i]);
            //getchar();
            cutsCoverSolution->rightSide[c_AuxSolution] = qnt - 1;
            c_AuxSolution++;
            free(c_menus);
            free(c_mais);
            free(S_barra);
            free(n_coef);
            free(n_el);
            free(aux);
        }
        int qnt_cuts_cover = 0;
        int *idc_cover = (int *)malloc(sizeof(int) * c_AuxSolution);
        for (j = 0; j < c_AuxSolution; j++)
        {
            //idc_cover[j] = 0;
            //printf("cut: %d  - %d\n",cutsCover->rightSide[i], newConstraintsSmall->rightSide[i]);
            //if (cutsCoverSolution->rightSide[j] != newConstraintsSmall->rightSide[i])
            //{
            idc_cover[j] = 1;
            qnt_cuts_cover++;
            //}
        }
        //printf("Quantidade de Cuts: %d\n", qnt_cuts_cover);
        // for (w = 0; w < c_AuxSolution; w++)
        // {
        //     for (k = cutsCoverSolution->ElementsConstraints[w]; k < cutsCoverSolution->ElementsConstraints[w + 1]; k++)
        //     {
        //         printf("%d ", cutsCoverSolution->Coefficients[k]);
        //     }
        //     printf("< = %d\n ", cutsCoverSolution->rightSide[w]);
        // }

        //constraintsFull = createCutsCoverGrasp(cutsCoverSolution, constraintsFull, newConstraintsSmall, idc_cover, i, qnt_cuts_cover);
        constraintsFull = createCutsCoverGrasp(cutsCoverSolution, constraintsFull, newConstraintsSmall, idc_cover, i, c_AuxSolution);
        free(cutsCoverSolution);
        free(fillBag);
        free(idc_cover);
        free(poolSolution);
    }
    free(cutsCover);
    free(newConstraintsSmall);
    free(intOrFloat);
    // printf("teste testando o teste: %d\n", qnt_cuts_cover);
    // getchar();
    // constraintsFull = createCutsCover(newConstraintsSmall, constraintsFull, cutsCover, idc_cover, qnt_cuts_cover);
    //printf("Number cuts cover: %d\n", qnt_cuts_cover);

    //free(constraintsSmall);

    // showStructFull(contraintsFull, renamedNameConstraints, nameVariables);

    return constraintsFull;
}

#endif // DEBUG