#include "structBasicsCuts.h"

cutSmall *AllocStrCutSmall(TCont cont, TNumberConstraints nConstraints, TNumberVariables nVariables)
{
    size_t size_cut = sizeof(cutSmall) +
                      sizeof(TCoefficients) * (cont) +
                      sizeof(TElements) * (cont) +
                      sizeof(TElementsConstraints) * (nConstraints + 1) +
                      sizeof(TRightSide) * (nConstraints) +
                      sizeof(TXAsterisc) * (nVariables) +
                      sizeof(TNumberConstraints) * (nConstraints);

    cutSmall *cut = (cutSmall *)malloc(size_cut);
    assert(cut != NULL);
    memset(cut, 0, size_cut);
    cut->Coefficients = (TCoefficients *)(cut + 1);
    cut->Elements = (TElements *)(cut->Coefficients + cont);
    cut->ElementsConstraints = (TElementsConstraints *)(cut->Elements + cont);
    cut->rightSide = (TRightSide *)(cut->ElementsConstraints + (nConstraints + 1));
    cut->xAsterisc = (TXAsterisc *)(cut->rightSide + (nConstraints));
    cut->originalConstraints = (TNumberConstraints *)(cut->xAsterisc + (nVariables));
    cut->numberVariables = nVariables;
    cut->numberConstraints = nConstraints;
    cut->cont = cont;
    return cut;
}

cutFull *AllocStrCutFull(TCont cont, TNumberConstraints nConstraints, TNumberVariables nVariables)
{

    cutFull *cut = (cutFull *)malloc(sizeof(cutFull));
    cut->Coefficients = (TCoefficientsFull *)malloc(sizeof(TCoefficientsFull) * cont);
    cut->Elements = (TElements *)malloc(sizeof(TElements) * cont);
    cut->ElementsConstraints = (TElementsConstraints *)malloc(sizeof(TElementsConstraints) * (nConstraints + 1));
    cut->rightSide = (TRightSideFull *)malloc(sizeof(TRightSideFull) * nConstraints);
    cut->xAsterisc = (TXAsteriscFull *)malloc(sizeof(TXAsteriscFull) * nVariables);
    cut->numberVariables = nVariables;
    cut->numberConstraints = nConstraints;
    cut->cont = cont;
    /*size_t size_cut = sizeof(cutFull) +
                      sizeof(TCoefficientsFull)*(cont) +
                      sizeof(TElements)*(cont) +
                      sizeof(TElementsConstraints)*(nConstraints+1) +
                      sizeof(TRightSideFull)*(nConstraints) +
                      sizeof(TXAsteriscFull)*(nVariables);


    cutFull *cut = (cutFull*)malloc(size_cut);
    assert(cut!=NULL);
    //memset(cut,0,size_cut);
    cut->Coefficients = (TCoefficientsFull*)(cut+1);
    cut->Elements = (TElements*)(cut->Coefficients + cont);
    cut->ElementsConstraints = (TElementsConstraints*)(cut->Elements + cont);
    cut->rightSide = (TRightSideFull*)(cut->ElementsConstraints + (nConstraints+1));
    cut->xAsterisc = (TXAsteriscFull*)(cut->rightSide + (nConstraints));
    cut->numberVariables = nVariables;
    cut->numberConstraints = nConstraints;
    cut->cont = cont;*/
    return cut;
}

void freeStrCutFull(cutFull *cut)
{
    free(cut->Coefficients);
    free(cut->Elements);
    free(cut->ElementsConstraints);
    free(cut->rightSide);
    free(cut->xAsterisc);
    free(cut);
}

cutCover *AllocStrCover(TCont cont, TNumberConstraints nConstraints)
{
    size_t size_cover = sizeof(cutCover) +
                        sizeof(TCoefficients) * (cont) +
                        sizeof(TElementsConstraints) * (nConstraints + 1) +
                        sizeof(TRightSide) * (nConstraints);
    cutCover *h_cover = (cutCover *)malloc(size_cover);
    assert(h_cover != NULL);
    memset(h_cover, 0, size_cover);
    h_cover->Coefficients = (TCoefficients *)(h_cover + 1);
    h_cover->ElementsConstraints = (TElementsConstraints *)(h_cover->Coefficients + cont);
    h_cover->rightSide = (TRightSide *)(h_cover->ElementsConstraints + (nConstraints + 1));
    h_cover->numberConstraints = nConstraints;
    h_cover->cont = cont;
    return h_cover;
}

cutCover *CopyCutToCover(cutSmall *h_cut)
{
    //Cover_gpu *cover = AllocationStructCover(h_cut->cont,h_cut->numberConstrains);
    cutCover *cover = AllocStrCover(h_cut->cont, h_cut->numberConstraints);
    int i;
    for (i = 0; i < h_cut->cont; i++) //VERIFY WHY H_CUT->CONT AND NOT NCONSTRAINTSINI
    {
        cover->Coefficients[i] = h_cut->Coefficients[i];
    }
    for (i = 0; i < h_cut->numberConstraints; i++)
    {
        cover->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
        cover->rightSide[i] = h_cut->rightSide[i];
    }
    cover->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
    return cover;
}

cutFull *fillStructPerLP(LinearProgram *lp, TNameConstraints **nameConstraints, TNameVariables **nameVariables)
{
    TNumberConstraints nConstraints, nConstraintsValided;
    TNumberVariables nVariables, nNonZeroValided = 0;

    nVariables = lp_cols(lp);
    nConstraints = lp_rows(lp);
    nConstraintsValided = 0;
    int szConst, i, j;
    TElements *idx = (TElements *)malloc(sizeof(TElements) * nVariables);
    double *coef = (double *)malloc(sizeof(double) * nVariables);
    for (i = 0; i < nConstraints; i++)
    {
        if ((lp_sense(lp, i) == 'L') || (lp_sense(lp, i) == 'G'))
        {
            nConstraintsValided++;
            for (j = 0; j < nVariables; j++)
            {
                coef[j] = 0.0;
                idx[j] = 0;
            }
            lp_row(lp, i, idx, coef);
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0)
                {
                    szConst++;
                    nNonZeroValided++;
                }
            }
            if (szConst < 2)
            {
                nConstraintsValided--;
                nNonZeroValided -= szConst;
            }
        }
        else
        {
            nConstraintsValided += 2;
            for (j = 0; j < nVariables; j++)
            {
                coef[j] = 0.0;
                idx[j] = 0;
            }
            lp_row(lp, i, idx, coef);
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0)
                {
                    szConst++;
                    nNonZeroValided += 2;
                }
            }
            if (szConst < 2)
            {
                nConstraintsValided -= 2;
                nNonZeroValided -= 2 * szConst;
            }
        }
    }
    double rhs;
    TXAsteriscFull *xTemp;
    lp_set_print_messages(lp, 0);
    lp_optimize_as_continuous(lp);
    xTemp = lp_x(lp);
    //saveSoluctionFrac(xTemp, numberVariables,nameInst,lp,0);

    cutFull *h_cut = AllocStrCutFull(nNonZeroValided, nConstraintsValided, nVariables);

    //TCoefficientsFull *v_aux = (TCoefficientsFull*)malloc(sizeof(TCoefficientsFull)*(nVariables + 1));

    for (i = 0; i < nVariables; i++)
    {
        coef[i] = 0.0;
        idx[i] = 0;
        h_cut->xAsterisc[i] = xTemp[i];
    }

    int aux = 0;
    h_cut->ElementsConstraints[0] = 0;
    int contador = 0;
    for (i = 0; i < nConstraints; i++)
    {
        if (lp_sense(lp, i) == 'L')
        {
            lp_row(lp, i, idx, coef);
            rhs = lp_rhs(lp, i);
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    szConst++;
                }
            }
            if (szConst < 2)
            {
                continue;
            }
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {

                    h_cut->Coefficients[aux] = coef[j];
                    h_cut->Elements[aux] = idx[j];
                    aux++;
                }

                coef[j] = 0.0;
                idx[j] = 0;
            }

            lp_row_name(lp, i, nameConstraints[contador]);
            h_cut->ElementsConstraints[contador + 1] = aux;
            h_cut->rightSide[contador] = rhs;
            contador++;
        }
        else if (lp_sense(lp, i) == 'G')
        {
            lp_row(lp, i, idx, coef);
            rhs = lp_rhs(lp, i);
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    szConst++;
                }
            }
            if (szConst < 2)
            {
                continue;
            }

            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    h_cut->Coefficients[aux] = -1 * coef[j];
                    h_cut->Elements[aux] = idx[j];
                    aux++;
                }

                coef[j] = 0.0;
                idx[j] = 0;
            }

            lp_row_name(lp, i, nameConstraints[contador]);

            h_cut->ElementsConstraints[contador + 1] = aux;
            h_cut->rightSide[contador] = -1 * rhs;
            contador++;
        }
        else
        {

            lp_row(lp, i, idx, coef);
            rhs = lp_rhs(lp, i);
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    szConst++;
                }
            }
            if (szConst < 2)
            {
                continue;
            }
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    h_cut->Coefficients[aux] = coef[j];
                    h_cut->Elements[aux] = idx[j];
                    aux++;
                }

                //coef[j] = 0.0;
                //idx[j] = 0;
            }
            char nameTemp[255];
            lp_row_name(lp, i, nameTemp);
            strcat(nameTemp, "_1");
            strcpy(nameConstraints[contador], nameTemp);
            h_cut->ElementsConstraints[contador + 1] = aux;
            h_cut->rightSide[contador] = rhs;
            contador++;

            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0.0)
                {
                    h_cut->Coefficients[aux] = -1 * coef[j];
                    h_cut->Elements[aux] = idx[j];
                    aux++;
                }

                coef[j] = 0.0;
                idx[j] = 0;
            }
            strcpy(nameTemp, "");
            lp_row_name(lp, i, nameTemp);
            strcat(nameTemp, "_2");
            strcpy(nameConstraints[contador], nameTemp);
            h_cut->ElementsConstraints[contador + 1] = aux;
            h_cut->rightSide[contador] = -1 * rhs;
            contador++;
        }
    }

    //free(xTemp);
    free(idx);
    free(coef);
    return h_cut;
}

int countContraintsValided(LinearProgram *lp)
{
    TNumberConstraints nConstraints, nConstraintsValided;
    TNumberVariables nVariables;

    nVariables = lp_cols(lp);
    nConstraints = lp_rows(lp);
    nConstraintsValided = 0;
    int szConst, i, j;
    TElements *idx = (TElements *)malloc(sizeof(TElements) * nVariables);
    TCoefficientsFull *coef = (TCoefficientsFull *)malloc(sizeof(TCoefficientsFull) * nVariables);
    for (i = 0; i < nConstraints; i++)
    {
        if ((lp_sense(lp, i) == 'L') || (lp_sense(lp, i) == 'G'))
        {
            nConstraintsValided++;
            for (j = 0; j < nVariables; j++)
            {
                coef[j] = 0.0;
                idx[j] = 0;
            }
            lp_row(lp, i, idx, coef);
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0)
                {
                    szConst++;
                }
            }
            if (szConst < 2)
            {
                nConstraintsValided--;
            }
        }
        else
        {
            nConstraintsValided += 2;
            for (j = 0; j < nVariables; j++)
            {
                coef[j] = 0.0;
                idx[j] = 0;
            }
            lp_row(lp, i, idx, coef);
            szConst = 0;
            for (j = 0; j < nVariables; j++)
            {
                if (coef[j] != 0)
                {
                    szConst++;
                }
            }
            if (szConst < 2)
            {
                nConstraintsValided -= 2;
            }
        }
    }
    free(idx);
    free(coef);
    return nConstraintsValided;
}

int verifyOfFloatIsInteger(TCoefficientsFull coef)
{
    //printf("%f %f", ceil(coef), coef);
    if (ceil(coef) == coef)
    {
        return 1;
    }

    return 0;
}

int *returnVectorTypeContraintsIntOrFloat(cutFull *constraints)
{
    int flag = 1, i, j, aux = 0;
    TCoefficientsFull coeffAux;
    int *vectorTypeInt = (int *)malloc(sizeof(int) * constraints->numberConstraints);
    for (i = 0; i < constraints->numberConstraints; i++)
    {
        flag = 1;
        for (j = constraints->ElementsConstraints[i]; j < constraints->ElementsConstraints[i + 1]; j++)
        {
            coeffAux = constraints->Coefficients[j];
            aux = verifyOfFloatIsInteger(coeffAux);
            if (aux == 0)
            {
                flag = 0;
                break;
            }
        }
        aux = verifyOfFloatIsInteger(constraints->rightSide[i]);
        if (aux == 0)
        {
            flag = 0;
        }
        if (flag)
        {
            vectorTypeInt[i] = 1;
        }
        else
        {
            vectorTypeInt[i] = 0;
        }
    }
    return vectorTypeInt;
}

cutSmall *reduceCutFullForCutSmall(cutFull *constraints, int *typeIntOrFloat, int precision)
{
    TNumberConstraints numberConstraintsSmall = 0;
    TCont cont = 0;
    int szNonZero, i, j, res, cont_aux;
    double aux;
    TElements el = 0;
    fflush(stdin);
    for (i = 0; i < constraints->numberConstraints; i++)
    {
        if (typeIntOrFloat[i])
        {
            numberConstraintsSmall++;
            cont += (constraints->ElementsConstraints[i + 1] - constraints->ElementsConstraints[i]);
        }
    }

    cutSmall *cSmall = AllocStrCutSmall(cont, numberConstraintsSmall, constraints->numberVariables);
    cont = 0;
    cont_aux = 0;
    cSmall->ElementsConstraints[0] = 0;
    for (i = 0; i < constraints->numberConstraints; i++)
    {
        if (typeIntOrFloat[i])
        {
            for (j = constraints->ElementsConstraints[i]; j < constraints->ElementsConstraints[i + 1]; j++)
            {
                cSmall->Coefficients[cont] = (TCoefficients)constraints->Coefficients[j];
                cSmall->Elements[cont] = constraints->Elements[j];
                cont++;
            }
            cSmall->rightSide[cont_aux] = (TRightSide)constraints->rightSide[i];
            cSmall->ElementsConstraints[cont_aux + 1] = cont;
            cSmall->originalConstraints[cont_aux] = i;
            cont_aux++;
        }
    }
    fflush(stdout);
    //printf("teste: %d %d \n", cont, numberConstraintsSmall);
    for (i = 0; i < constraints->numberVariables; i++)
    {
        cSmall->xAsterisc[i] = (int)(constraints->xAsterisc[i] * (double)precision);
    }
    return cSmall;
}

void showStructSmall(cutSmall *constraintsSmall, TNameConstraints **nameConstraints, TNameVariables **nameVariables)
{
    int i, j, ct, el, flag;
    for (i = 0; i < constraintsSmall->numberConstraints; i++)
    {
        flag = 0;
        ct = constraintsSmall->originalConstraints[i];
        printf("%s:\t ", nameConstraints[ct]);
        for (j = constraintsSmall->ElementsConstraints[i]; j < constraintsSmall->ElementsConstraints[i + 1]; j++)
        {
            el = constraintsSmall->Elements[j];
            if (constraintsSmall->Coefficients[j] > 0)
            {
                if (flag == 0)
                {
                    printf("%d %s ", constraintsSmall->Coefficients[j], nameVariables[el]);
                    flag = 1;
                }
                else
                {
                    printf("+ %d %s ", constraintsSmall->Coefficients[j], nameVariables[el]);
                }
            }
            else
            {
                flag = 1;
                printf("%d %s ", constraintsSmall->Coefficients[j], nameVariables[el]);
            }
        }
        printf("<= %d\n", constraintsSmall->rightSide[i]);
    }
}

void showStructFull(cutFull *constraintsFull, TNameConstraints **nameConstraints, TNameVariables **nameVariables)
{
    int i, j, el, flag;
    for (i = 0; i < constraintsFull->numberConstraints; i++)
    {
        flag = 0;
        printf("%s:\t ", nameConstraints[i]);
        for (j = constraintsFull->ElementsConstraints[i]; j < constraintsFull->ElementsConstraints[i + 1]; j++)
        {
            el = constraintsFull->Elements[j];
            if (constraintsFull->Coefficients[j] > 0)
            {
                if (flag == 0)
                {
                    printf("%.2f %s ", constraintsFull->Coefficients[j], nameVariables[el]);
                    flag = 1;
                }
                else
                {
                    printf("+ %.2f %s ", constraintsFull->Coefficients[j], nameVariables[el]);
                }
            }
            else
            {
                flag = 1;
                printf("%.2f %s ", constraintsFull->Coefficients[j], nameVariables[el]);
            }
        }
        printf("<= %.2f\n", constraintsFull->rightSide[i]);
        
    }
}

TCoefficients cutMaxDivisorCommonVector(TCoefficients coefs[], TNumberVariables nElem)
{
    TCoefficients n = coefs[nElem];
    TCoefficients mdc = 1;
    while (nElem > 0)
    {
        TCoefficients m = coefs[nElem - 1];
        mdc = cutMaxDivisorCommonRec(m, n);
        n = mdc;
        nElem--;
    }
    TCoefficients m = coefs[nElem];
    mdc = cutMaxDivisorCommonRec(m, n);
    n = mdc;
    return n;
}

/*calculates the maximum common divisor for two integers.*/
TCoefficients cutMaxDivisorCommonRec(TCoefficients m, TCoefficients n)
{

    TCoefficients t = 0;
    m = m < 0 ? -m : m; /* abs(u) */
    n = n < 0 ? -n : n;
    if (m < n)
    {
        t = m;
        m = n;
        n = t;
    }

    if (n == 0)
        return m;
    else
    {
        TCoefficients resto = m % n;
        return cutMaxDivisorCommonRec(n, resto);
    }
}

int *calcCover(cutCover *h_cover, int *h_solution, int tolerance)
{
    int *fillBag = (int *)malloc(sizeof(int) * h_cover->numberConstraints);
    int i, j, k, counter = 0;

    for (i = 0; i < h_cover->numberConstraints; i++)
    {
        counter = 0;
        for (k = h_cover->ElementsConstraints[i]; k < h_cover->ElementsConstraints[i + 1]; k++)
        {
            counter += h_cover->Coefficients[k] * h_solution[k];
        }
        if (counter + tolerance <= h_cover->rightSide[i])
        {
            for (k = h_cover->ElementsConstraints[i]; k < h_cover->ElementsConstraints[i + 1]; k++)
            {
                if (h_solution[k] == 0)
                {
                    counter += h_cover->Coefficients[k];
                    h_solution[k] = 1;
                    if (counter + tolerance > h_cover->rightSide[i])
                    {
                        break;
                    }
                }
            }
        }
        fillBag[i] = counter;

        //printf("Fill bag: %d Capacity bag: %d\n", fillBag[i], h_cover->rightSide[i]);

        //getchar();
    }
    //printf("number Constraints: %d qnt %d \n", h_cover->numberConstraints, qnt_Cover_per_Thread);
    return fillBag;
}

TNameConstraints **renamedNameConstraints(TNameConstraints **nameConstraints, int typeContraints, TNumberConstraints szNewConstraints, TNumberConstraints szCuts, TNumberConstraints lastCut)
{
    int i;
    TNameConstraints **newNameConstraints = (TNameConstraints **)malloc(szNewConstraints * sizeof(TNameConstraints *));
    for (i = 0; i < szNewConstraints; i++)
    {
        newNameConstraints[i] = (TNameConstraints *)malloc(255 * sizeof(TNameConstraints));
        if (i < szNewConstraints - szCuts)
        {
            strcpy(newNameConstraints[i], nameConstraints[i]);
        }
    }
    for (i = szNewConstraints - szCuts; i < szNewConstraints; i++)
    {
        if (typeContraints == 1)
        {
            char name[255];
            sprintf(name, "CG1(%d)", lastCut);
            lastCut++;
            strcpy(newNameConstraints[i], name);
        }
        if (typeContraints == 2)
        {
            char name[255];
            sprintf(name, "CG2(%d)", lastCut);
            lastCut++;
            strcpy(newNameConstraints[i], name);
        }
        if (typeContraints == 3)
        {
            char name[255];
            sprintf(name, "CC(%d)", lastCut);
            lastCut++;
            strcpy(newNameConstraints[i], name);
        }
    }
    for (i = 0; i < szNewConstraints - szCuts; i++)
    {
        free(nameConstraints[i]);
    }
    free(nameConstraints);
    return newNameConstraints;
}

cutFull *removeNegativeCoefficientsAndSort(cutFull *constraintsOriginal, int *convertVector, int precision)
{
    int i, j;
    //convertVector = (int*)(malloc(sizeof(int)*h_cut->cont));
    //  convertCoef = (int*)(malloc(sizeof(int)*h_cut->cont));
    int qntX = constraintsOriginal->numberVariables;
    int qntNegative = 0;
    TRightSideFull rhs = 0;
    int el = 0;
    for (i = 0; i < constraintsOriginal->cont; i++)
    {
        if (constraintsOriginal->Coefficients[i] < 0)
        {
            qntNegative++;
        }
    }
    qntNegative += constraintsOriginal->numberVariables;
    cutFull *newConstraints = AllocStrCutFull(constraintsOriginal->cont, constraintsOriginal->numberConstraints, qntNegative);

    for (i = 0; i < constraintsOriginal->numberVariables; i++)
    {
        newConstraints->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }

    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        rhs = constraintsOriginal->rightSide[i];
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            if (constraintsOriginal->Coefficients[j] < 0)
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j] * (-1);

                rhs += newConstraints->Coefficients[j];
                newConstraints->Elements[j] = qntX;
                el = constraintsOriginal->Elements[j];
                newConstraints->xAsterisc[qntX] = 1 - constraintsOriginal->xAsterisc[el];
                convertVector[qntX - constraintsOriginal->numberVariables] = constraintsOriginal->Elements[j];
                //convertCoef[qntX-h_cut->numberVariables] = Cut_new->Coefficients[j];
                qntX++;
            }
            else
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j];
                newConstraints->Elements[j] = constraintsOriginal->Elements[j];
                //convertVector[j] = 0 ;
                //convertCoef[j] = 0;
            }
        }
        newConstraints->rightSide[i] = rhs;
        newConstraints->ElementsConstraints[i] = constraintsOriginal->ElementsConstraints[i];
    }
    newConstraints->ElementsConstraints[i] = constraintsOriginal->ElementsConstraints[i];
    SortByCoefficients(newConstraints);
    freeStrCutFull(constraintsOriginal);
    return newConstraints;
}

cutFull *returnVariablesOriginals(cutFull *constraintsOriginal, int *convertVector, int precision, int nVariablesInitial)
{

    cutFull *newConstraints = AllocStrCutFull(constraintsOriginal->cont, constraintsOriginal->numberConstraints, nVariablesInitial);
    int i, j, el;
    TRightSideFull rhs;
    for (i = 0; i < nVariablesInitial; i++)
    {
        newConstraints->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        rhs = constraintsOriginal->rightSide[i];
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            el = constraintsOriginal->Elements[j];
            if (el >= nVariablesInitial)
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j] * (-1);

                rhs -= constraintsOriginal->Coefficients[j];
                newConstraints->Elements[j] = convertVector[el - nVariablesInitial];
                newConstraints->xAsterisc[newConstraints->Elements[j]] = 1 - constraintsOriginal->xAsterisc[el]; //erro aqui
            }
            else
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j];
                newConstraints->Elements[j] = constraintsOriginal->Elements[j];
            }
        }

        newConstraints->rightSide[i] = rhs;
        newConstraints->ElementsConstraints[i] = constraintsOriginal->ElementsConstraints[i];
    }
    newConstraints->ElementsConstraints[i] = constraintsOriginal->ElementsConstraints[i];
    freeStrCutFull(constraintsOriginal);
    return newConstraints;
}

void SortByCoefficients(cutFull *h_cut)
{
    int i = 0;
    for (i = 0; i < h_cut->numberConstraints; i++)
    {
        quicksortCof(h_cut->Coefficients, h_cut->Elements, h_cut->ElementsConstraints[i], h_cut->ElementsConstraints[i + 1]);
    }
}

void quicksortCof(TCoefficientsFull *values, int *idc, int began, int end)
{
    int i, j;
    TCoefficientsFull pivo, aux;
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

            aux = idc[i];
            idc[i] = idc[j];
            idc[j] = aux;

            i++;
            j--;
        }
    }
    if (j > began)
        quicksortCof(values, idc, began, j + 1);
    if (i < end)
        quicksortCof(values, idc, i, end);
}

int verifyRepeatCuts(cutFull *constraintsOriginal, int cutOriginal, int cutCreate)
{
    int i, j, aux = 1;
    int szOri = constraintsOriginal->ElementsConstraints[cutOriginal + 1] - constraintsOriginal->ElementsConstraints[cutOriginal];
    int szCre = constraintsOriginal->ElementsConstraints[cutCreate + 1] - constraintsOriginal->ElementsConstraints[cutCreate];
    if ((szOri != szCre) || (constraintsOriginal->rightSide[cutOriginal] != constraintsOriginal->rightSide[cutCreate]))
    {
        return 0;
    }
    for (i = constraintsOriginal->ElementsConstraints[cutOriginal]; i < constraintsOriginal->ElementsConstraints[cutOriginal + 1]; i++)
    {
        aux = 0;
        for (j = constraintsOriginal->ElementsConstraints[cutCreate]; j < constraintsOriginal->ElementsConstraints[cutCreate + 1]; j++)
        {
            if ((constraintsOriginal->Coefficients[i] == constraintsOriginal->Coefficients[j]) && (constraintsOriginal->Elements[i] == constraintsOriginal->Elements[j]))
            {
                aux = 1;
            }
        }
        if (aux == 0)
        {
            return 0;
        }
    }
    return 1;
}

int *vectorNonRepeteadNonDominated(cutFull *constraintsOriginal, int nConstraintsInitial)
{
    int i, j, k = 0, qntRepeat = 0;
    int *isRepeat = (int *)malloc(sizeof(int) * constraintsOriginal->numberConstraints);
    for (i = 0; i < nConstraintsInitial; i++)
    {
        isRepeat[i] = 0;
    }
    int *v_aux = (int *)malloc(sizeof(int) * constraintsOriginal->numberVariables);
    for (i = nConstraintsInitial; i < constraintsOriginal->numberConstraints; i++)
    {
        k = 0;
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            v_aux[k] = constraintsOriginal->Coefficients[j];
            k++;
        }
        v_aux[k] = constraintsOriginal->rightSide[i];
        int mdc = cutMaxDivisorCommonVector(v_aux, k);
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            constraintsOriginal->Coefficients[j] = constraintsOriginal->Coefficients[j] / mdc;
        }
        constraintsOriginal->rightSide[i] = constraintsOriginal->rightSide[i] / mdc;
        for (j = 0; j < i; j++)
        {
            isRepeat[i] = verifyRepeatCuts(constraintsOriginal, j, i);
            qntRepeat += isRepeat[i];
            if (isRepeat[i] == 1)
            {
                break;
            }
        }
    }
    printf("Number cuts repeat: %d\n", qntRepeat);
    free(v_aux);
    return isRepeat;
}

void insertConstraintsLP(LinearProgramPtr lp, cutFull *constraintsOriginal, int nConstrainsInitial, char **nameConstraints)
{
    int i, j, w;
    int *idx;
    double *Coef;
    int *cutsNonInserted = vectorNonRepeteadNonDominated(constraintsOriginal, nConstrainsInitial);
    for (i = nConstrainsInitial; i < constraintsOriginal->numberConstraints; i++)
    {
        if (cutsNonInserted[i] == 1)
        {
            continue;
        }
        int sz = constraintsOriginal->ElementsConstraints[i + 1] - constraintsOriginal->ElementsConstraints[i];
        idx = (int *)malloc(sizeof(int) * sz);
        Coef = (double *)malloc(sizeof(double) * sz);
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            idx[w] = constraintsOriginal->Elements[j];
            Coef[w] = constraintsOriginal->Coefficients[j];
            w++;
        }
        double rhs = constraintsOriginal->rightSide[i];

        //if(h_cut->typeConstraints[i]==LPC_CCOVER){
        /*         int v = *counterCuts;
        if(h_cut->typeConstraints[i] == LPC_CGGPU){
            sprintf(name, "CGGPU1(%d)",v);
        }else if(h_cut->typeConstraints[i] == LPC_CGGPUR2){
            sprintf(name, "CGGPU2(%d)",v);
        }else{
            sprintf(name, "CCOVER(%d)",v);
        } */

        //        printf("ok: %s\n",name);
        lp_add_row(lp, sz, idx, Coef, nameConstraints[i], 'L', rhs);
        //}
        w = 0;
        free(idx);
        free(Coef);
    }

    free(cutsNonInserted);
}

int verifyCutsValidatedPerSolutionInteger(cutFull *constraintsOriginal, int cut, double *sol, char **nameVariables)
{
    double lhs = 0;
    int j, el;
    for (j = constraintsOriginal->ElementsConstraints[cut]; j < constraintsOriginal->ElementsConstraints[cut + 1]; j++)
    {
        el = constraintsOriginal->Elements[j];
        lhs += constraintsOriginal->Coefficients[j] * (sol[el]);
        //f (sol[el] != 0.0)
        //    printf("%lf * %lf  %s  = %lf %lf\n", constraintsOriginal->Coefficients[j], sol[el], nameVariables[el], constraintsOriginal->Coefficients[j] * (sol[el]), lhs);
    }
    //printf("lhs: %lf rhs %lf\n", lhs, constraintsOriginal->rightSide[cut]);
    //getchar();
    // if(lhs<0){
    //     lhs += 1e-3;
    // }else{
    //     lhs -= 1e-3;
    // }

    if (lhs - 1e-2<= constraintsOriginal->rightSide[cut])
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

double *readSolFile(const char *name, int nVariables)
{
    double *sol = (double *)malloc(sizeof(double) * nVariables);
    int idx, i;
    double value, aux;
    char n[255];
    for (i = 0; i < nVariables; i++)
    {
        sol[i] = 0;
    }
    FILE *f;
    f = fopen(name, "r");
    if (f != NULL)
    {
        int flag = 0;
        while (!feof(f))
        {
            if (flag == 0)
            {
                fscanf(f, "%[^\n]", n);
                flag = 1;
            }
            fscanf(f, "%d %s %lf %lf\n", &idx, n, &value, &aux);
            sol[idx] = value;
        }
    }
    return sol;
}