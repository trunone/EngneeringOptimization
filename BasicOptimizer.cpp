/*
 *	BasicOptimizer.cpp
 *
 *	Description:
 *		Basic optimization class includes
 *			1. Basic line search algorithm
 *			2. Gradient computing
 *			3. Gradient algorithm
 *			4. Conjugate Gradient algorithm
 *			5. QuasiNewton algorithm
 *
 *
 * 	History:
 *	 	Author			Date			Modify Reason
 *		----------------------------------------------------------------
 *		Chi-Yi Tsai		2012/09/05		File Creation
 *
 *
 */
#include "BasicOptimizer.h"

CBasicOptimizer::CBasicOptimizer() {
    this->Self = this;
}

CBasicOptimizer::~CBasicOptimizer() {
}

void CBasicOptimizer::SetLineSearchMethod(LineSearch Line) {
    this->Line = Line;
}

void CBasicOptimizer::SetSolverMethod(SolverMethod OptMethod) {
    this->OptMethod = OptMethod;
}

void CBasicOptimizer::SetDifferentialMethod(DifferentialMethod Differential) {
    this->Differential = Differential;
}

void CBasicOptimizer::SetConjugateFormula(ConjugateFormula Formula) {
    this->CFormula = Formula;
}

void CBasicOptimizer::SetQuasiNewtonFormula(QuasiNewtonFormula Formula) {
    this->QNFormula = Formula;
}

LineSearch CBasicOptimizer::GetLineSearchMethod(void) {
    return this->Line;
}

SolverMethod CBasicOptimizer::GetSolverMethod(void) {
    return this->OptMethod;
}

DifferentialMethod CBasicOptimizer::GetDifferentialMethod(void) {
    return this->Differential;
}

ConjugateFormula CBasicOptimizer::GetConjugateFormula(void) {
    return this->CFormula;
}

QuasiNewtonFormula CBasicOptimizer::GetQuasiNewtonFormula(void) {
    return this->QNFormula;
}

double CBasicOptimizer::OptOneDimensionFunSolver(pfnScaleCostFun pfun, double fX) {
    this->pfun = pfun;
    this->nDim = 1;
    this->f_g = new double[this->nDim]();
    this->f_g[0] = -1;
    this->pfX = &fX;

    CLineSearch *LineSearchMethod = new CLineSearch(this);

    LineSearchMethod->Find_Interval();
// 	LineSearchMethod->LowerBound = 0;
// 	LineSearchMethod->UpperBound = 2;
// 	LineSearchMethod->Interval = 2;
    switch(Line) {
    case GoldenSection:
        this->pfX = LineSearchMethod->GoldenSelection();
        break;
    case Fibonacci:
        this->pfX = LineSearchMethod->Fibonacci();
        break;
    default:
        return 0;
        break;
    }
    delete LineSearchMethod;

    return pfun(this->pfX);
}

double CBasicOptimizer::OptScaleCostFunSolver(pfnScaleCostFun pfun, double *pfX, int nDim) {
    this->pfun = pfun;
    this->nDim = nDim;
    this->pfX = pfX;
    this->f_g = new double(this->nDim);

    CSolverMethod *SolverMethod = new CSolverMethod(this);

    switch(OptMethod) {
    case SteepestDecent:
        SolverMethod->SteepestDecentMethod();
        break;
    case ConjugateGradient:
        SolverMethod->ConjugateGradientMethod();
        break;
    case QuasiNewton:
        return 0;
        break;
    default:
        return 0;
        break;
    }
    delete SolverMethod;
    memcpy(pfX, this->pfX, this->nDim*sizeof(double));

    cout<<this->pfX[0]<<" "<<this->pfX[1]<<endl;

    return pfun(this->pfX);
}
