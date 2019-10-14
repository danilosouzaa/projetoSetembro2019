#include "structSolution.h"

solutionStructCG1 *allocStructSolCG1(int nRuns)
{
    size_t size_solution = sizeof(solutionStructCG1) +
                           sizeof(TSNumerator) * (nRuns) +
                           sizeof(TSDenominator) * (nRuns) +
                           sizeof(TSConstraints) * (nRuns);
    solutionStructCG1 *sol;
    sol = (solutionStructCG1 *)malloc(size_solution);
    assert(sol != NULL);
    memset(sol, 0, size_solution);
    sol->sNumerator = (TSNumerator *)(sol + 1);
    sol->sDenominator = (TSDenominator *)(sol->sNumerator + (nRuns));
    sol->sConstraints = (TSConstraints *)(sol->sDenominator + (nRuns));
    return sol;
}

solutionStructCG2 *allocStructSolCG2(int numberMaxConst, int nRuns)
{
    size_t size_solution = sizeof(solutionStructCG2) +
                           sizeof(TSNumerator) * (nRuns) +
                           sizeof(TSDenominator) * (nRuns) +
                           sizeof(TSNumerator) * (nRuns) +
                           sizeof(TSDenominator) * (nRuns) +
                           sizeof(TSConstraints) * (numberMaxConst * nRuns) +
                           sizeof(TSPosition) * (nRuns);
    solutionStructCG2 *sol;
    sol = (solutionStructCG2 *)malloc(size_solution);
    assert(sol != NULL);
    memset(sol, 0, size_solution);
    sol->sNumerator1 = (TSNumerator *)(sol + 1);
    sol->sDenominator1 = (TSDenominator *)(sol->sNumerator1 + nRuns);
    sol->sNumerator2 = (TSNumerator *)(sol->sDenominator1 + nRuns);
    sol->sDenominator2 = (TSDenominator *)(sol->sNumerator2 + nRuns);
    sol->sConstraints = (TSConstraints *)(sol->sDenominator2 + nRuns);
    sol->sPosition = (TSPosition *)(sol->sConstraints + (numberMaxConst * nRuns));
    return sol;
}

TCoefficients returnK(TSTempCreateSolution fa_zero)
{
    TSTempCreateSolution value;
    TCoefficients arround;
    value = 1 / fa_zero;
    arround = value;
    if (fabsl(value - (TSTempCreateSolution)arround) <= 1e-6)
    {
        arround--;
    }
    return arround;
}

int verifyDominanceCG(TCoefficients *v1, TRightSide rhs1, TCoefficients *v2, TRightSide rhs2, TNumberVariables sz)
{
    TCoefficients v1_temp[sz];
    TCoefficients v2_temp[sz];
    int i, cont = 0;
    for (i = 0; i < sz; i++)
    {
        v1_temp[i] = v1[i] * rhs2;
        v2_temp[i] = v2[i] * rhs1;
    }
    rhs2 = rhs2 * rhs1;
    rhs1 = rhs2;
    for (i = 0; i < sz; i++)
    {
        if (v1_temp[i] < v2_temp[i])
        {
            return 0;
        }
        else if (v1_temp[i] > v2_temp[i])
        {
            cont++;
        }
    }
    if (cont > 0)
        return 1;
    else
        return 0;
}

cutFull *createCutsStrongPhaseOne(cutSmall *constraintsSmall, cutFull *constraintsOriginal, solutionStructCG1 *solution, TNumberConstraints nCuts, int precision, int nRuns, TNumberConstraints nC_initial, TCont cont_ini)
{
    TCoefficients *Coef1 = (TCoefficients *)malloc(sizeof(TCoefficients) * (constraintsSmall->numberVariables));
    TCoefficients *CoefXEqualOne = (TCoefficients *)malloc(sizeof(TCoefficients) * (constraintsSmall->numberVariables));
    TCoefficients *Coef_ref = (TCoefficients *)malloc(sizeof(TCoefficients) * (constraintsSmall->numberVariables));
    TSTempCreateSolution *Coef1_mult = (TSTempCreateSolution *)malloc(sizeof(TSTempCreateSolution) * (constraintsSmall->numberVariables));
    TSTempCreateSolution *p_frac = (TSTempCreateSolution *)malloc(sizeof(TSTempCreateSolution) * (constraintsSmall->numberVariables));

    TCoefficients *v_aux = (TCoefficients *)malloc(sizeof(TCoefficients) * (constraintsSmall->numberVariables + 1));

    TSViolation *violation = (TSViolation *)malloc(sizeof(TSViolation) * nCuts);
    TSViolation *violation_mult = (TSViolation *)malloc(sizeof(TSViolation) * nCuts);

    TCoefficients *Coefs_temp = (TCoefficients *)malloc(sizeof(TCoefficients) * (nCuts * constraintsSmall->numberVariables));
    TElements *Elements_temp = (TElements *)malloc(sizeof(TElements) * (nCuts * constraintsSmall->numberVariables));
    TSPosition *Pos_el_temp = (TSPosition *)malloc(sizeof(TSPosition) * (nCuts + 1));
    TRightSide *rhs_temp = (TRightSide *)malloc(sizeof(TRightSide) * (nCuts));
    int i, k, j, ite, aux = 0;
    TSNumerator n1;
    TSDenominator d1;
    TElements el;
    TNumberConstraints constraint, originalConstraint, tam, c_aux = 0;
    TRightSide rhs, rhs_ref, rhsXEqualOne;
    TSTempCreateSolution rhs_mult, lhs;
    int cont_aux = 0, contXEqualOne = 0;

    Pos_el_temp[0] = 0;

    for (i = 0; i < nRuns; i++)
    {
        if (solution->sConstraints[i] != -1) // se gerou um corte
        {
            // getchar();
            memset(Coef1, 0, sizeof(TCoefficients) * constraintsSmall->numberVariables);
            memset(Coef_ref, 0, sizeof(TCoefficients) * constraintsSmall->numberVariables);
            memset(CoefXEqualOne, 0, sizeof(TCoefficients) * constraintsSmall->numberVariables);
            rhsXEqualOne = 0;
            for (j = 0; j < constraintsSmall->numberVariables; j++)
            {
                Coef1_mult[j] = 0;
                p_frac[j] = 0;
            }
            lhs = 0;
            aux = 0;
            tam = 0;
            violation[c_aux] = 0;
            violation_mult[c_aux] = 0;
            constraint = solution->sConstraints[i];
            n1 = solution->sNumerator[i];
            d1 = solution->sDenominator[i];
            rhs = constraintsSmall->rightSide[constraint] * n1;
            rhs_mult = ((TSTempCreateSolution)constraintsSmall->rightSide[constraint] * (TSTempCreateSolution)n1) / (TSTempCreateSolution)d1;
            if ((rhs < 0 && d1 > 0) || (rhs > 0 && d1 < 0))
            {
                rhs = (rhs / d1) - 1;
            }
            else
            {
                rhs = rhs / d1;
            }
            rhs_mult = rhs_mult - rhs;
            aux = 0;
            // getchar();
            originalConstraint = constraintsSmall->originalConstraints[constraint];
            for (j = constraintsOriginal->ElementsConstraints[originalConstraint]; j < constraintsOriginal->ElementsConstraints[originalConstraint + 1]; j++)
            {
                el = constraintsOriginal->Elements[j];
                if (constraintsSmall->xAsterisc[el] != precision)
                {
                    Coef1[el] = (TCoefficients)constraintsOriginal->Coefficients[j] * n1;
                    Coef1_mult[el] = ((TSTempCreateSolution)constraintsOriginal->Coefficients[j] * (TSTempCreateSolution)n1) / (TSTempCreateSolution)d1;
                    if ((Coef1[el] < 0 && d1 > 0) || (Coef1[el] > 0 && d1 < 0))
                    {

                        Coef1[el] = (Coef1[el] / d1) - 1;
                    }
                    else
                    {
                        Coef1[el] = (Coef1[el] / d1);
                    }
                    p_frac[el] = Coef1_mult[el] - Coef1[el];
                    lhs += (TSTempCreateSolution)Coef1[el] * (TSTempCreateSolution)constraintsSmall->xAsterisc[el];
                }
                else
                {
                    if (constraintsOriginal->Coefficients[j] != 0)
                    {
                        CoefXEqualOne[el] = constraintsOriginal->Coefficients[j];
                        rhsXEqualOne += CoefXEqualOne[el];
                        contXEqualOne++;
                    }
                }
            }

            violation[c_aux] = lhs - (rhs * precision);
            if (rhs_mult - 1e-8 > 0)
            {
                k = returnK(rhs_mult);
                assert(rhs_mult <= (1.0 / ((long double)k)) + 1e-8);
                assert(rhs_mult >= (1.0 / ((long double)k + 1)) - 1e-8);
                TNumberVariables *sizeRows = (TNumberVariables *)malloc((k + 1) * sizeof(TNumberVariables));
                TNumberVariables **mat = (TNumberVariables **)malloc((k + 1) * sizeof(TNumberVariables *));
                for (ite = 0; ite < (k + 1); ite++)
                    mat[ite] = (TNumberVariables *)malloc(constraintsOriginal->numberVariables * sizeof(TNumberVariables));

                for (ite = 0; ite < constraintsOriginal->numberVariables; ite++)
                {
                    if (p_frac[ite] - 1e-6 <= rhs_mult)
                    {
                        mat[0][aux] = ite;
                        aux++;
                    }
                }
                sizeRows[0] = aux;
                aux = 0;
                for (j = 1; j <= k; j++)
                {
                    for (ite = 0; ite < constraintsOriginal->numberVariables; ite++)
                    {
                        if (Coef1_mult[ite] == 0.0)
                        {
                            continue;
                        }
                        if ((p_frac[ite] - 1e-6 >= (rhs_mult + (j - 1.0) * (1.0 - rhs_mult) / (TSTempCreateSolution)k)) && (p_frac[ite] /* + 1e-6*/ <= (rhs_mult + j * (1.0 - rhs_mult) / (TSTempCreateSolution)k)))
                        {
                            mat[j][aux] = ite;
                            aux++;
                        }
                    }
                    sizeRows[j] = aux;
                    aux = 0;
                }
                for (ite = 0; ite < sizeRows[0]; ite++)
                {
                    el = mat[0][ite];
                    if (Coef1_mult[el] < 0)
                    {

                        aux = Coef1_mult[el] - 1;
                    }
                    else
                    {
                        aux = Coef1_mult[el];
                    }
                    Coef_ref[el] = aux * (k + 1);
                }
                for (j = 1; j <= k; j++)
                {
                    for (ite = 0; ite < sizeRows[j]; ite++)
                    {
                        el = mat[j][ite];

                        if (Coef1_mult[el] < 0)
                        {
                            aux = Coef1_mult[el] - 1;
                        }
                        else
                        {
                            aux = Coef1_mult[el];
                        }
                        Coef_ref[el] = aux * (k + 1) + j;
                    }
                }
                rhs_ref = (k + 1) * rhs;
                lhs = 0;
                for (ite = 0; ite < constraintsOriginal->numberVariables; ite++)
                {
                    lhs += (TSTempCreateSolution)Coef_ref[ite] * (TSTempCreateSolution)constraintsSmall->xAsterisc[ite];
                }
                violation_mult[c_aux] = lhs - (rhs_ref * precision);

                for (ite = 0; ite < k + 1; ite++)
                {
                    free(mat[ite]);
                }
                free(mat);
                free(sizeRows);
                tam = 0;
                if (verifyDominanceCG(Coef_ref, rhs_ref, Coef1, rhs, constraintsOriginal->numberVariables) == 1)
                {
                    for (j = 0; j < constraintsOriginal->numberVariables; j++)
                    {
                        if (Coef_ref[j] != 0)
                        {
                            v_aux[tam] = Coef_ref[j];
                            tam++;
                        }
                        if (CoefXEqualOne[j] != 0)
                        {
                            v_aux[tam] = CoefXEqualOne[j];
                            tam++;
                        }
                    }
                    v_aux[tam] = rhs_ref + rhsXEqualOne;
                    TCoefficients mdc = cutMaxDivisorCommonVector(v_aux, tam);
                    for (j = 0; j < constraintsOriginal->numberVariables; j++)
                    {

                        if (Coef_ref[j] != 0)
                        {
                            Coefs_temp[cont_aux] = Coef_ref[j] / mdc;
                            Elements_temp[cont_aux] = j;
                            cont_aux++;
                        }
                        if (CoefXEqualOne[j] != 0)
                        {
                            Coefs_temp[cont_aux] = CoefXEqualOne[j] / mdc;
                            Elements_temp[cont_aux] = j;
                            cont_aux++;
                        }
                    }
                    Pos_el_temp[c_aux + 1] = cont_aux;
                    rhs_temp[c_aux] = (rhs_ref + rhsXEqualOne) / mdc;
                }
                else
                {
                    for (j = 0; j < constraintsOriginal->numberVariables; j++)
                    {
                        if (Coef1[j] != 0)
                        {
                            v_aux[tam] = Coef1[j];
                            tam++;
                        }
                        if (CoefXEqualOne[j] != 0)
                        {
                            v_aux[tam] = CoefXEqualOne[j];
                            tam++;
                        }
                    }
                    v_aux[tam] = rhs + rhsXEqualOne;
                    TCoefficients mdc = cutMaxDivisorCommonVector(v_aux, tam);
                    for (j = 0; j < constraintsOriginal->numberVariables; j++)
                    {
                        if (Coef1[j] != 0)
                        {
                            Coefs_temp[cont_aux] = Coef1[j] / mdc;
                            Elements_temp[cont_aux] = j;
                            cont_aux++;
                        }
                        if (CoefXEqualOne[j] != 0)
                        {
                            Coefs_temp[cont_aux] = CoefXEqualOne[j] / mdc;
                            Elements_temp[cont_aux] = j;
                            cont_aux++;
                        }
                    }
                    Pos_el_temp[c_aux + 1] = cont_aux;
                    rhs_temp[c_aux] = (rhs + rhsXEqualOne) / mdc;
                }
            }
            else
            {
                tam = 0;
                for (j = 0; j < constraintsOriginal->numberVariables; j++)
                {
                    if ((Coef1[j] != 0))
                    {
                        v_aux[tam] = Coef1[j];
                        tam++;
                    }
                    if (CoefXEqualOne[j] != 0)
                    {
                        v_aux[tam] = CoefXEqualOne[j];
                        tam++;
                    }
                }
                v_aux[tam] = rhs;
                TCoefficients mdc = cutMaxDivisorCommonVector(v_aux, tam);
                for (j = 0; j < constraintsOriginal->numberVariables; j++)
                {
                    if (Coef1[j] != 0)
                    {
                        Coefs_temp[cont_aux] = Coef1[j] / mdc;
                        Elements_temp[cont_aux] = j;
                        cont_aux++;
                    }
                    if (CoefXEqualOne[j] != 0)
                    {
                        Coefs_temp[cont_aux] = CoefXEqualOne[j] / mdc;
                        Elements_temp[cont_aux] = j;
                        cont_aux++;
                    }
                }
                Pos_el_temp[c_aux + 1] = cont_aux;
                rhs_temp[c_aux] = (rhs + rhsXEqualOne) / mdc;
            }
            c_aux++;
        }
    }

    free(Coef1);
    free(Coef1_mult);
    free(CoefXEqualOne);
    free(Coef_ref);
    free(violation);
    free(violation_mult);
    free(p_frac);
    free(v_aux);

    //DEPENDE DE COMO VAI SER
    //int cont_ini = constraintsOriginal->cont;
    //int nC_initial  = constraintsOriginal->numberConstraints;
    cutFull *cuts_generated;
    cuts_generated = AllocStrCutFull(cont_ini + cont_aux, nC_initial + nCuts, constraintsOriginal->numberVariables);
    //assert(!cuts_generated);

    for (i = 0; i < constraintsOriginal->numberVariables; i++)
    {
        cuts_generated->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }
    cuts_generated->ElementsConstraints[0] = 0;

    for (i = 0; i < nC_initial; i++)
    {
        cuts_generated->rightSide[i] = constraintsOriginal->rightSide[i];
        cuts_generated->ElementsConstraints[i + 1] = constraintsOriginal->ElementsConstraints[i + 1];
    }
    aux = 1;

    for (i = nC_initial; i < cuts_generated->numberConstraints; i++)
    {

        cuts_generated->rightSide[i] = rhs_temp[i - nC_initial];
        cuts_generated->ElementsConstraints[i + 1] = Pos_el_temp[aux] + cont_ini;
        aux++; //+ h_cut->cont;
    }

    for (i = 0; i < cont_ini; i++)
    {
        cuts_generated->Coefficients[i] = constraintsOriginal->Coefficients[i];
        cuts_generated->Elements[i] = constraintsOriginal->Elements[i];
    }
    for (i = cont_ini; i < cuts_generated->cont; i++)
    {
        cuts_generated->Coefficients[i] = Coefs_temp[i - cont_ini];
        cuts_generated->Elements[i] = Elements_temp[i - cont_ini];
    }

    int *validated = (int *)malloc(sizeof(int) * cuts_generated->numberConstraints);
    memset(validated, 0, sizeof(int) * cuts_generated->numberConstraints);
    aux = 1;
    cont_aux = 0;
    int p1, p2, minus_elements = 0;
    k = 0;
    for (i = 0; i < cuts_generated->numberConstraints - 1; i++)
    {
        for (j = i + 1; j < cuts_generated->numberConstraints; j++)
        {
            if (validated[j] == 0)
            {
                if ((cuts_generated->ElementsConstraints[i + 1] - cuts_generated->ElementsConstraints[i]) != (cuts_generated->ElementsConstraints[j + 1] - cuts_generated->ElementsConstraints[j]) || (cuts_generated->rightSide[i] != cuts_generated->rightSide[j]))
                {
                    aux = 0;
                }
                else
                {
                    p1 = cuts_generated->ElementsConstraints[i];
                    p2 = cuts_generated->ElementsConstraints[j];
                    aux = 1;
                    for (k = 0; k < cuts_generated->ElementsConstraints[i + 1] - cuts_generated->ElementsConstraints[i]; k++)
                    {

                        if ((cuts_generated->Coefficients[p1 + k] != cuts_generated->Coefficients[p2 + k]) || (cuts_generated->Elements[p1 + k] != cuts_generated->Elements[p2 + k]))
                        {
                            aux = 0;

                            break;
                        }
                    }
                }
                if ((aux == 1) && (j >= constraintsOriginal->numberConstraints))
                {
                    validated[j] = 1;
                    minus_elements += cuts_generated->ElementsConstraints[j + 1] - cuts_generated->ElementsConstraints[j];
                    cont_aux++;
                }
            }
        }
    }

    //    printf("Number of repeat Strengthening: %d \n", cont_aux);

    cutFull *new_h_cut;

    new_h_cut = AllocStrCutFull(cuts_generated->cont - minus_elements, cuts_generated->numberConstraints - cont_aux, cuts_generated->numberVariables);
    aux = 0;
    cont_aux = 0;
    new_h_cut->ElementsConstraints[0] = 0;
    for (i = 0; i < cuts_generated->numberConstraints; i++)
    {
        if (validated[i] == 0)
        {
            //            if(i>=h_cut->numberConstrains)
            //            {
            //                printf("violation %f\n",violation[i-h_cut->numberConstrains]);
            //            }
            new_h_cut->rightSide[aux] = cuts_generated->rightSide[i];
            for (j = cuts_generated->ElementsConstraints[i]; j < cuts_generated->ElementsConstraints[i + 1]; j++)
            {
                new_h_cut->Coefficients[cont_aux] = cuts_generated->Coefficients[j];
                new_h_cut->Elements[cont_aux] = cuts_generated->Elements[j];
                cont_aux++;
            }
            new_h_cut->ElementsConstraints[aux + 1] = cont_aux;

            aux++;
        }
    }
    for (i = 0; i < new_h_cut->numberVariables; i++)
    {
        new_h_cut->xAsterisc[i] = cuts_generated->xAsterisc[i];
    }
    //
    //    printf("End phase 1!\n");
    free(Coefs_temp);
    free(Elements_temp);
    free(Pos_el_temp);
    free(rhs_temp);
    free(validated);
    freeStrCutFull(cuts_generated);

    return new_h_cut;
}

cutFull *createCutsStrongPhaseTwo(cutSmall *constraintsSmall, cutFull *constraintsOriginal, solutionStructCG2 *solution, int numberMaxConst, int nCuts, int precision, int nThreads, int nBlocks)
{

    TSViolation *value_violation = (TSViolation *)malloc(sizeof(TSViolation) * nCuts);
    TSViolation *value_violation_r1 = (TSViolation *)malloc(sizeof(TSViolation) * nCuts);
    TSViolation *violation_mult = (TSViolation *)malloc(sizeof(TSViolation) * nCuts);
    TCoefficients *Coef1 = (TCoefficients *)malloc(sizeof(TCoefficients) * (constraintsOriginal->numberVariables));
    TCoefficients *Coef2 = (TCoefficients *)malloc(sizeof(TCoefficients) * (constraintsOriginal->numberVariables));
    TCoefficients *Coef_ref = (TCoefficients *)malloc(sizeof(TCoefficients) * (constraintsOriginal->numberVariables));
    TCoefficients *CoefXEqualOne = (TCoefficients *)malloc(sizeof(TCoefficients) * (constraintsOriginal->numberVariables));
    long double *Coef_mult1 = (long double *)malloc(sizeof(long double) * constraintsOriginal->numberVariables);
    long double *p_frac = (long double *)malloc(sizeof(long double) * constraintsOriginal->numberVariables);
    TCoefficients *Coefs_temp = (TCoefficients *)malloc(sizeof(TCoefficients) * (nCuts * constraintsOriginal->numberVariables));
    TElements *Elements_temp = (TElements *)malloc(sizeof(TElements) * (nCuts * constraintsOriginal->numberVariables));
    TElementsConstraints *Pos_el_temp = (TElementsConstraints *)malloc(sizeof(TElementsConstraints) * (nCuts + 1));
    TRightSide *rhs_temp = (TRightSide *)malloc(sizeof(TRightSide) * (nCuts));
    TRightSide rhsXEqualOne;
    int *v_aux = (int *)malloc(sizeof(int) * (constraintsOriginal->numberVariables + 1));
    long double rhs_mult1, rhs_mult2;
    int temp, k_R2;
    int i, k, j, n_1, d_1, n_2, d_2, qnt_1, el, rest_a, rhs1 = 0, rhs2 = 0, aux = 0, cont_aux = 0, ite, tam = 0;
    cutFull *cuts_generated;
    TNumberConstraints originalConstraint;
    int contXEqualOne = 0;
    Pos_el_temp[0] = 0;
    for (i = 0; i < nThreads; i++)
    {
        for (k = 0; k < nBlocks; k++)
        {
            if (solution->sConstraints[0 + (i + k * nThreads) * numberMaxConst] != -1)
            {
                for (ite = 0; ite < constraintsOriginal->numberVariables; ite++)
                {
                    Coef_mult1[ite] = 0;
                    p_frac[ite] = 0;
                }

                rhs1 = 0;
                rhs2 = 0;
                rhs_mult1 = 0;
                rhs_mult2 = 0;
                qnt_1 = solution->sPosition[i + k * nThreads];
                n_1 = solution->sNumerator1[i + k * nThreads];
                d_1 = solution->sDenominator1[i + k * nThreads];
                n_2 = solution->sNumerator2[i + k * nThreads];
                d_2 = solution->sDenominator2[i + k * nThreads];
                if ((d_1 == 0) || (d_2 == 0))
                    continue;

                /*for(rest_a = 0; rest_a <= qnt_1; rest_a++)
                {
                    rhs1 += constraintsSmall->rightSide[ solution->sConstraints[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ] * n_1;
                }

                for(rest_a = qnt_1 + 1; rest_a < numberMaxConst; rest_a++)
                {
                    rhs2 += constraintsSmall->rightSide[ solution->sConstraints[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ] * n_2;
                }
                rhs_mult1 = ((TRightSideFull)rhs1/(TRightSideFull)d_1 + (TRightSideFull)rhs2/(TRightSideFull)d_2);// RHS Double
                if ((rhs1*d_2 + rhs2*d_1)/(d_1*d_2) < 0)
                {
                    temp = (rhs1*d_2 + rhs2*d_1)/(d_1*d_2) - 1;
                }
                else
                {
                    temp = (rhs1*d_2 + rhs2*d_1)/(d_1*d_2);
                }
                rhs_mult2 = rhs_mult1 - temp;
                rhs1 = 0;
                rhs2 = 0;
                rhs_mult1 = 0;
                rhs_mult2 = 0;*/

                memset(Coef1, 0, sizeof(TCoefficients) * constraintsOriginal->numberVariables);
                memset(Coef2, 0, sizeof(TCoefficients) * constraintsOriginal->numberVariables);
                memset(Coef_ref, 0, sizeof(TCoefficients) * constraintsOriginal->numberVariables);
                memset(CoefXEqualOne, 0, sizeof(TCoefficients) * constraintsOriginal->numberVariables);
                rhsXEqualOne = 0;
                contXEqualOne = 0;
                for (rest_a = 0; rest_a <= qnt_1; rest_a++)
                {
                    originalConstraint = constraintsSmall->originalConstraints[solution->sConstraints[rest_a + i * numberMaxConst + k * numberMaxConst * nThreads]];
                    for (j = constraintsOriginal->ElementsConstraints[originalConstraint]; j < constraintsOriginal->ElementsConstraints[originalConstraint + 1]; j++)
                    {
                        el = constraintsOriginal->Elements[j];
                        if (constraintsSmall->xAsterisc[el] != precision)
                        {
                            Coef1[el] += constraintsOriginal->Coefficients[j] * n_1;
                        }
                        else
                        {
                            if (constraintsOriginal->Coefficients[j] != 0)
                            {
                                CoefXEqualOne[el] += constraintsOriginal->Coefficients[j];
                                rhsXEqualOne += CoefXEqualOne[el];
                            }
                        }
                    }
                    rhs1 += constraintsSmall->rightSide[solution->sConstraints[rest_a + i * numberMaxConst + k * numberMaxConst * nThreads]] * n_1;
                }

                for (rest_a = qnt_1 + 1; rest_a < numberMaxConst; rest_a++)
                {
                    originalConstraint = constraintsSmall->originalConstraints[solution->sConstraints[rest_a + i * numberMaxConst + k * numberMaxConst * nThreads]];
                    for (j = constraintsOriginal->ElementsConstraints[originalConstraint]; j < constraintsOriginal->ElementsConstraints[originalConstraint + 1]; j++)
                    {
                        el = constraintsOriginal->Elements[j];
                        if (constraintsSmall->xAsterisc[el] != precision)
                        {
                            Coef2[el] += constraintsOriginal->Coefficients[j] * n_2;
                        }
                        else
                        {
                            if (constraintsOriginal->Coefficients[j] != 0)
                            {
                                CoefXEqualOne[el] += constraintsOriginal->Coefficients[j];
                                rhsXEqualOne += CoefXEqualOne[el];
                            }
                        }
                    }
                    rhs2 += constraintsSmall->rightSide[solution->sConstraints[rest_a + i * numberMaxConst + k * numberMaxConst * nThreads]] * n_2;
                }

                value_violation[cont_aux] = 0;
                value_violation_r1[cont_aux] = 0;
                violation_mult[cont_aux] = 0;

                for (j = 0; j < constraintsOriginal->numberVariables; j++)
                {
                    Coef_mult1[j] = (long double)Coef1[j] / (long double)d_1 + (long double)Coef2[j] / (long double)d_2; //Coef Double
                    if ((Coef1[j] * d_2 + Coef2[j] * d_1) / (d_1 * d_2) < 0)
                    {
                        temp = ((Coef1[j] * d_2 + Coef2[j] * d_1) / (d_1 * d_2)) - 1;
                    }
                    else
                    {
                        temp = (Coef1[j] * d_2 + Coef2[j] * d_1) / (d_1 * d_2);
                    }

                    p_frac[j] = Coef_mult1[j] - temp; //Coef Frac
                    if (temp != 0)
                    {
                        Coef1[j] = temp;
                        value_violation[cont_aux] += (temp)*constraintsSmall->xAsterisc[j];
                    }
                    else
                    {
                        Coef1[j] = 0;
                    }
                }
                rhs_mult1 = ((double)rhs1 / (double)d_1 + (double)rhs2 / (double)d_2); // RHS Double
                if ((rhs1 * d_2 + rhs2 * d_1) / (d_1 * d_2) < 0)
                {
                    temp = (rhs1 * d_2 + rhs2 * d_1) / (d_1 * d_2) - 1;
                }
                else
                {
                    temp = (rhs1 * d_2 + rhs2 * d_1) / (d_1 * d_2);
                }
                rhs_temp[cont_aux] = temp;
                value_violation[cont_aux] -= (rhs_temp[cont_aux]) * precision;
                rhs_mult2 = rhs_mult1 - temp; //Frac RHS
                rhs1 = rhs_temp[cont_aux];
                if (rhs_mult2 > 0)
                {
                    k_R2 = returnK(rhs_mult2);

                    assert(rhs_mult2 <= (1.0 / ((long double)k_R2)) + 1e-8);
                    assert(rhs_mult2 >= (1.0 / ((long double)k_R2 + 1)) - 1e-8);
                    int *sizeRows = (int *)malloc((k_R2 + 1) * sizeof(int));
                    memset(sizeRows, 0, (k_R2 + 1) * sizeof(int));
                    int **mat = (int **)malloc((k_R2 + 1) * sizeof(int *));
                    for (ite = 0; ite < (k_R2 + 1); ite++)
                        mat[ite] = (int *)malloc(constraintsOriginal->numberVariables * sizeof(int));
                    temp = 0;
                    for (ite = 0; ite < constraintsOriginal->numberVariables; ite++)
                    {
                        if ((p_frac[ite] - 1e-6 <= rhs_mult2) && (Coef_mult1[ite] != 0))
                        {
                            mat[0][temp] = ite;
                            temp++;
                        }
                    }
                    sizeRows[0] = temp;
                    temp = 0;
                    for (j = 1; j <= k_R2; j++)
                    {
                        for (ite = 0; ite < constraintsOriginal->numberVariables; ite++)
                        {
                            if (Coef_mult1[ite] == 0)
                            {
                                continue;
                            }

                            if ((p_frac[ite] - 1e-6 > (rhs_mult2 + (j - 1.0) * (1.0 - rhs_mult2) / (long double)k_R2)) && (p_frac[ite] + 1e-6 <= (rhs_mult2 + (j * (1.0 - rhs_mult2) / (long double)k_R2))))
                            {
                                mat[j][temp] = ite;
                                temp++;
                            }
                        }
                        sizeRows[j] = temp;
                        temp = 0;
                    }
                    for (ite = 0; ite < sizeRows[0]; ite++)
                    {
                        el = mat[0][ite];
                        if (Coef_mult1[el] < 0)
                        {

                            temp = Coef_mult1[el] - 1;
                        }
                        else
                        {
                            temp = Coef_mult1[el];
                        }
                        Coef_ref[el] = temp * (k_R2 + 1);
                    }
                    for (j = 1; j <= k_R2; j++)
                    {
                        for (ite = 0; ite < sizeRows[j]; ite++)
                        {
                            el = mat[j][ite];

                            if (Coef_mult1[el] < 0)
                            {
                                temp = Coef_mult1[el] - 1;
                            }
                            else
                            {
                                temp = Coef_mult1[el];
                            }
                            Coef_ref[el] = temp * (k_R2 + 1) + j;
                        }
                    }
                    rhs_temp[cont_aux] = (k_R2 + 1) * rhs_temp[cont_aux];
                    int lhs = 0;

                    for (ite = 0; ite < constraintsOriginal->numberVariables; ite++)
                    {
                        lhs += (long double)Coef_ref[ite] * (long double)constraintsSmall->xAsterisc[ite];
                    }
                    violation_mult[cont_aux] = lhs - (rhs_temp[cont_aux] * precision);
                    for (ite = 0; ite < k_R2 + 1; ite++)
                    {
                        free(mat[ite]);
                    }
                    free(mat);
                    free(sizeRows);
                    tam = 0;
                    if (verifyDominanceCG(Coef_ref, rhs_temp[cont_aux], Coef1, rhs1, constraintsOriginal->numberVariables) == 1)
                    {
                        for (j = 0; j < constraintsOriginal->numberVariables; j++)
                        {
                            if (Coef_ref[j] != 0)
                            {
                                v_aux[tam] = Coef_ref[j];
                                tam++;
                            }
                            if (CoefXEqualOne[j] != 0)
                            {
                                v_aux[tam] = CoefXEqualOne[j];
                                tam++;
                                contXEqualOne++;
                            }
                        }
                        v_aux[tam] = rhs_temp[cont_aux] + rhsXEqualOne;
                        int mdc = cutMaxDivisorCommonVector(v_aux, tam);
                        for (j = 0; j < constraintsOriginal->numberVariables; j++)
                        {
                            if (Coef_ref[j] != 0)
                            {
                                Coefs_temp[aux] = Coef_ref[j] / mdc;
                                Elements_temp[aux] = j;
                                aux++;
                            }
                            if (CoefXEqualOne[j] != 0)
                            {
                                Coefs_temp[aux] = CoefXEqualOne[j] / mdc;
                                Elements_temp[aux] = j;
                                aux++;
                            }
                        }
                        Pos_el_temp[cont_aux + 1] = aux;
                        rhs_temp[cont_aux] = (rhs_temp[cont_aux] + rhsXEqualOne) / mdc;
                    }
                    else
                    {
                        tam = 0;
                        for (j = 0; j < constraintsOriginal->numberVariables; j++)
                        {
                            if (Coef1[j] != 0)
                            {
                                v_aux[tam] = Coef1[j];
                                tam++;
                            }
                            if (CoefXEqualOne[j] != 0)
                            {
                                v_aux[tam] = CoefXEqualOne[j];
                                tam++;
                                contXEqualOne++;
                            }
                        }
                        v_aux[tam] = rhs1 + rhsXEqualOne;
                        int mdc = cutMaxDivisorCommonVector(v_aux, tam);
                        for (j = 0; j < constraintsOriginal->numberVariables; j++)
                        {
                            if (Coef1[j] != 0)
                            {
                                Coefs_temp[aux] = Coef1[j] / mdc;
                                Elements_temp[aux] = j;
                                aux++;
                            }
                            if (CoefXEqualOne[j] != 0)
                            {
                                Coefs_temp[aux] = CoefXEqualOne[j] / mdc;
                                Elements_temp[aux] = j;
                                aux++;
                            }
                        }
                        Pos_el_temp[cont_aux + 1] = aux;
                        rhs_temp[cont_aux] = (rhs1 + rhsXEqualOne) / mdc;
                    }
                }
                else
                {
                    tam = 0;
                    for (j = 0; j < constraintsOriginal->numberVariables; j++)
                    {
                        if (Coef1[j] != 0)
                        {
                            v_aux[tam] = Coef1[j];
                            tam++;
                        }
                        if (CoefXEqualOne[j != 0])
                        {
                            v_aux[tam] = CoefXEqualOne[j];
                            tam++;
                            contXEqualOne++;
                        }
                    }
                    v_aux[tam] = rhs1;
                    int mdc = cutMaxDivisorCommonVector(v_aux, tam);
                    for (j = 0; j < constraintsOriginal->numberVariables; j++)
                    {
                        if (Coef1[j] != 0)
                        {
                            Coefs_temp[aux] = Coef1[j] / mdc;
                            Elements_temp[aux] = j;
                            aux++;
                        }
                        if (CoefXEqualOne[j] != 0)
                        {
                            Coefs_temp[aux] = CoefXEqualOne[j] / mdc;
                            Elements_temp[aux] = j;
                            aux++;
                        }
                    }
                    Pos_el_temp[cont_aux + 1] = aux;
                    rhs_temp[cont_aux] = (rhs1 + rhsXEqualOne) / mdc;
                }
                cont_aux++;
            }
        }
    }

    free(Coef1);
    free(Coef2);
    free(Coef_mult1);
    free(Coef_ref);
    free(CoefXEqualOne);
    free(violation_mult);
    free(value_violation);
    free(value_violation_r1);
    free(v_aux);
    free(p_frac);
    cont_aux = aux;
    cuts_generated = AllocStrCutFull(constraintsOriginal->cont + cont_aux, constraintsOriginal->numberConstraints + nCuts, constraintsOriginal->numberVariables);
    for (i = 0; i < cuts_generated->numberVariables; i++)
    {
        cuts_generated->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }
    cuts_generated->ElementsConstraints[0] = 0;
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        cuts_generated->rightSide[i] = constraintsOriginal->rightSide[i];
        cuts_generated->ElementsConstraints[i + 1] = constraintsOriginal->ElementsConstraints[i + 1];
    }
    aux = 1;

    for (i = constraintsOriginal->numberConstraints; i < cuts_generated->numberConstraints; i++)
    {

        cuts_generated->rightSide[i] = rhs_temp[i - constraintsOriginal->numberConstraints];
        cuts_generated->ElementsConstraints[i + 1] = Pos_el_temp[aux] + constraintsOriginal->cont;
        aux++; //+ h_cut->cont;
    }
    for (i = 0; i < constraintsOriginal->cont; i++)
    {
        cuts_generated->Coefficients[i] = constraintsOriginal->Coefficients[i];
        cuts_generated->Elements[i] = constraintsOriginal->Elements[i];
    }
    for (i = constraintsOriginal->cont; i < cuts_generated->cont; i++)
    {
        cuts_generated->Coefficients[i] = Coefs_temp[i - constraintsOriginal->cont];
        cuts_generated->Elements[i] = Elements_temp[i - constraintsOriginal->cont];
    }
    int *validated = (int *)malloc(sizeof(int) * cuts_generated->numberConstraints);
    memset(validated, 0, sizeof(int) * cuts_generated->numberConstraints);
    aux = 1;
    cont_aux = 0;
    int p1, p2, minus_elements = 0;
    k = 0;
    for (i = 0; i < cuts_generated->numberConstraints - 1; i++)
    {
        for (j = i + 1; j < cuts_generated->numberConstraints; j++)
        {
            if (validated[j] == 0)
            {
                if ((cuts_generated->ElementsConstraints[i + 1] - cuts_generated->ElementsConstraints[i]) != (cuts_generated->ElementsConstraints[j + 1] - cuts_generated->ElementsConstraints[j]) || (cuts_generated->rightSide[i] != cuts_generated->rightSide[j]))
                {
                    aux = 0;
                }
                else
                {
                    p1 = cuts_generated->ElementsConstraints[i];
                    p2 = cuts_generated->ElementsConstraints[j];
                    aux = 1;
                    for (k = 0; k < cuts_generated->ElementsConstraints[i + 1] - cuts_generated->ElementsConstraints[i]; k++)
                    {

                        if ((cuts_generated->Coefficients[p1 + k] != cuts_generated->Coefficients[p2 + k]) || (cuts_generated->Elements[p1 + k] != cuts_generated->Elements[p2 + k]))
                        {
                            aux = 0;

                            break;
                        }
                    }
                }
                if ((aux == 1) && (j >= constraintsOriginal->numberConstraints))
                {
                    validated[j] = 1;
                    minus_elements += cuts_generated->ElementsConstraints[j + 1] - cuts_generated->ElementsConstraints[j];
                    cont_aux++;
                }
            }
        }
    }

    cutFull *new_h_cut;

    new_h_cut = AllocStrCutFull(cuts_generated->cont - minus_elements, cuts_generated->numberConstraints - cont_aux, cuts_generated->numberVariables);
    aux = 0;
    cont_aux = 0;
    new_h_cut->ElementsConstraints[0] = 0;
    for (i = 0; i < cuts_generated->numberConstraints; i++)
    {
        if (validated[i] == 0)
        {
            new_h_cut->rightSide[aux] = cuts_generated->rightSide[i];
            for (j = cuts_generated->ElementsConstraints[i]; j < cuts_generated->ElementsConstraints[i + 1]; j++)
            {
                new_h_cut->Coefficients[cont_aux] = cuts_generated->Coefficients[j];
                new_h_cut->Elements[cont_aux] = cuts_generated->Elements[j];
                cont_aux++;
            }
            new_h_cut->ElementsConstraints[aux + 1] = cont_aux;

            aux++;
        }
    }
    for (i = 0; i < new_h_cut->numberVariables; i++)
    {
        new_h_cut->xAsterisc[i] = cuts_generated->xAsterisc[i];
    }
    freeStrCutFull(cuts_generated);
    free(validated);
    free(Coefs_temp);
    free(Elements_temp);
    free(Pos_el_temp);
    free(rhs_temp);

    return new_h_cut;
}

int verifyRepeated(cutFull *originalConstraints, int posCover)
{
    int i, j, k, cont = 0;
    int sz = originalConstraints->ElementsConstraints[posCover + 1] - originalConstraints->ElementsConstraints[posCover];
    for (i = 0; i < posCover; i++)
    {
        for (j = originalConstraints->ElementsConstraints[posCover]; j < originalConstraints->ElementsConstraints[posCover + 1]; j++)
        {
            if ((sz == (originalConstraints->ElementsConstraints[i + 1] - originalConstraints->ElementsConstraints[i])) && (originalConstraints->rightSide[posCover] == originalConstraints->rightSide[i]))
            {
                for (k = originalConstraints->ElementsConstraints[i]; k < originalConstraints->ElementsConstraints[i + 1]; k++)
                {
                    if ((originalConstraints->Elements[j] == originalConstraints->Elements[k]) && (originalConstraints->Coefficients[j] == originalConstraints->Coefficients[k]))
                    {
                        cont++;
                        break;
                    }
                }
            }
        }
        if (cont == sz)
        {
            return 0;
        }
        else
        {
            cont = 0;
        }
    }
    return 1;
}

cutFull *createCutsCover(cutSmall *constraintsSmall, cutFull *constraintsOriginal, cutCover *cutsCover, int *idc_Cover, int nCuts)
{
    if (nCuts <= 0)
    {
        return constraintsOriginal;
    }
    int i = 0, j = 0, cont = 0, contConstraints = 0;
    for (i = 0; i < cutsCover->numberConstraints; i++)
    {
        if (idc_Cover[i] == 1)
        {
            for (j = cutsCover->ElementsConstraints[i]; j < cutsCover->ElementsConstraints[i + 1]; j++)
            {
                if (cutsCover->Coefficients[j] != 0)
                {
                    cont++;
                }
            }
            contConstraints++;
        }
    }
    cutFull *outCutsNew = AllocStrCutFull(constraintsOriginal->cont + cont, constraintsOriginal->numberConstraints + contConstraints, constraintsOriginal->numberVariables);
    outCutsNew->ElementsConstraints[0] = constraintsOriginal->ElementsConstraints[0];
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        outCutsNew->rightSide[i] = constraintsOriginal->rightSide[i];
        outCutsNew->ElementsConstraints[i + 1] = constraintsOriginal->ElementsConstraints[i + 1];
    }
    for (i = 0; i < constraintsOriginal->cont; i++)
    {
        outCutsNew->Elements[i] = constraintsOriginal->Elements[i];
        outCutsNew->Coefficients[i] = constraintsOriginal->Coefficients[i];
    }
    for (i = 0; i < constraintsOriginal->numberVariables; i++)
    {
        outCutsNew->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }

    int aux = constraintsOriginal->numberConstraints;
    int c_aux = constraintsOriginal->cont;
    outCutsNew->ElementsConstraints[aux] = c_aux;
    for (i = 0; i < cutsCover->numberConstraints; i++)
    {
        if (idc_Cover[i] == 1)
        {
            outCutsNew->rightSide[aux] = cutsCover->rightSide[i];
            for (j = cutsCover->ElementsConstraints[i]; j < cutsCover->ElementsConstraints[i + 1]; j++)
            {
                if (cutsCover->Coefficients[j] != 0)
                {
                    outCutsNew->Elements[c_aux] = constraintsSmall->Elements[j];
                    outCutsNew->Coefficients[c_aux] = cutsCover->Coefficients[j];
                    c_aux++;
                }
            }

            outCutsNew->ElementsConstraints[aux + 1] = c_aux;
            aux++;
        }
    }

    // int *validated = (int*)malloc(sizeof(int)*outCutsNew->numberConstraints);
    // validated[0] = 1;
    // for (i=1;i<outCutsNew->numberConstraints;i++){
    //     validated[i] = verifyRepeated(outCutsNew,i);
    // }
    // cont = 0, contConstraints = 0;
    // for(i = 0; i<outCutsNew->numberConstraints;i++){
    //     if(validated[i]==1){
    //         contConstraints++;
    //         for(j = outCutsNew->ElementsConstraints[i];j<outCutsNew->ElementsConstraints[i+1];j++){
    //             cont++;
    //         }

    //     }
    // }
    // if(contConstraints == outCutsNew->numberConstraints){
    //     freeStrCutFull(constraintsOriginal);
    //     free(validated);
    //     return outCutsNew;
    // }
    // cutFull *cutsNewNoRepetead = AllocStrCutFull(cont,contConstraints,outCutsNew->numberVariables);
    // cont=0, contConstraints = 0 ;
    // cutsNewNoRepetead->ElementsConstraints[0] = 0;
    // for(i=0;i<outCutsNew->numberConstraints;i++){
    //     if(validated[i]==1){
    //         cutsNewNoRepetead->rightSide[contConstraints] = outCutsNew->rightSide[i];
    //         for(j = outCutsNew->ElementsConstraints[i];j<outCutsNew->ElementsConstraints[i+1];j++){
    //             cutsNewNoRepetead->Elements[cont] = outCutsNew->Elements[j];
    //             cutsNewNoRepetead->Coefficients[cont] = outCutsNew->Coefficients[j];
    //             cont++;
    //         }
    //         cutsNewNoRepetead->ElementsConstraints[contConstraints+1] = cont;
    //         contConstraints++;

    //     }
    // }
    // for(i=0;i<cutsNewNoRepetead->numberVariables;i++){
    //     cutsNewNoRepetead->xAsterisc[i] = outCutsNew->xAsterisc[i];
    // }
    // free(validated);
    freeStrCutFull(constraintsOriginal);
    //freeStrCutFull(outCutsNew);

    return outCutsNew;
}

cutFull *createCutsCoverGrasp(cutCover *cutsCover, cutFull *constraintsOriginal, cutSmall *constraintsSmall, int *idc_Cover, int constraint, int nCuts)
{
    if (nCuts <= 0)
    {  
        return constraintsOriginal;
    }
    long int i = 0, j = 0, cont = 0, contConstraints = 0;
    int sz = constraintsOriginal->ElementsConstraints[constraint + 1] - constraintsOriginal->ElementsConstraints[constraint];
   //printf("%d %d %d\n", sz, nCuts, cutsCover->cont);

    for (i = 0; i < nCuts; i++)
    {
        if (idc_Cover[i] == 1)
        {   
            contConstraints++;
            for (j = cutsCover->ElementsConstraints[i]; j < cutsCover->ElementsConstraints[i + 1]; j++)
            {
                if (cutsCover->Coefficients[j] != 0)
                {
                    cont++;
                }
            }
        }
    }
   // printf("cont: %d\n", cont);
    cutFull *outCutsNew = AllocStrCutFull(constraintsOriginal->cont + cont, constraintsOriginal->numberConstraints + contConstraints, constraintsOriginal->numberVariables);
    outCutsNew->ElementsConstraints[0] = constraintsOriginal->ElementsConstraints[0];
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        outCutsNew->rightSide[i] = constraintsOriginal->rightSide[i];
        outCutsNew->ElementsConstraints[i + 1] = constraintsOriginal->ElementsConstraints[i + 1];
    }
    for (i = 0; i < constraintsOriginal->cont; i++)
    {
        outCutsNew->Elements[i] = constraintsOriginal->Elements[i];
        outCutsNew->Coefficients[i] = constraintsOriginal->Coefficients[i];
    }
    for (i = 0; i < constraintsOriginal->numberVariables; i++)
    {
        outCutsNew->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }

    int aux = constraintsOriginal->numberConstraints;

    int c_aux = constraintsOriginal->cont;
    //printf("nCuts %d el con: %d c_aux: %d\n", nCuts, outCutsNew->ElementsConstraints[aux],c_aux);

    //outCutsNew->ElementsConstraints[aux] = c_aux;
    for (i = 0; i < nCuts; i++)
    {
        if (idc_Cover[i] == 1)
        {
            int c_XSolution = constraintsSmall->ElementsConstraints[constraint];
            outCutsNew->rightSide[aux] = cutsCover->rightSide[i];
            for (j = cutsCover->ElementsConstraints[i]; j < cutsCover->ElementsConstraints[i + 1]; j++)
            {
                if (cutsCover->Coefficients[j] != 0)
                {
                    outCutsNew->Elements[c_aux] = constraintsSmall->Elements[c_XSolution];
                    outCutsNew->Coefficients[c_aux] = cutsCover->Coefficients[j];
                    c_aux++;
                }
                c_XSolution++;
            }
            // printf("c_aux: %d\t", c_aux);
            outCutsNew->ElementsConstraints[aux + 1] = c_aux;
            aux++;
        }
        // printf("\n");
    }
    // printf("terminou\n");
    int *validated = (int *)malloc(sizeof(int) * outCutsNew->numberConstraints);
    validated[0] = 1;
    for (i = 1; i < outCutsNew->numberConstraints; i++)
    {
        validated[i] = verifyRepeated(outCutsNew, i);
        // if(validated[i]==0)
        //     printf("Validated = %d\n",validated[i]);
    }
    cont = 0, contConstraints = 0;
    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        if (validated[i] == 1)
        {
            contConstraints++;
            for (j = outCutsNew->ElementsConstraints[i]; j < outCutsNew->ElementsConstraints[i + 1]; j++)
            {
                cont++;
            }
        }
    }
   // printf("xxx: %d %d\n", contConstraints,  constraintsOriginal->numberConstraints);
    if (contConstraints == outCutsNew->numberConstraints)
    {   
       // printf("AQUI!!\n");
        freeStrCutFull(constraintsOriginal);
        free(validated);
        return outCutsNew;
    }
    cutFull *cutsNewNoRepetead = AllocStrCutFull(cont, contConstraints, outCutsNew->numberVariables);
    cont = 0, contConstraints = 0;
    cutsNewNoRepetead->ElementsConstraints[0] = 0;
    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        if (validated[i] == 1)
        {
            cutsNewNoRepetead->rightSide[contConstraints] = outCutsNew->rightSide[i];
            for (j = outCutsNew->ElementsConstraints[i]; j < outCutsNew->ElementsConstraints[i + 1]; j++)
            {
                cutsNewNoRepetead->Elements[cont] = outCutsNew->Elements[j];
                cutsNewNoRepetead->Coefficients[cont] = outCutsNew->Coefficients[j];
                cont++;
            }
            cutsNewNoRepetead->ElementsConstraints[contConstraints + 1] = cont;
            contConstraints++;
        }
    }
    for (i = 0; i < cutsNewNoRepetead->numberVariables; i++)
    {
        cutsNewNoRepetead->xAsterisc[i] = outCutsNew->xAsterisc[i];
    }
    free(validated);
    freeStrCutFull(constraintsOriginal);
    freeStrCutFull(outCutsNew);

    return cutsNewNoRepetead;
}
