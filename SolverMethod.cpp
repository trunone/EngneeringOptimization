#include "BasicOptimizer.h"

CBasicOptimizer::CSolverMethod::CSolverMethod(CBasicOptimizer *p) : Super(p),
    FinalRange(0.00001) {
    n = Super->nDim;
}

void CBasicOptimizer::CSolverMethod::SteepestDecentMethod() {
    for(;;) {
        DifferentialMethod();

        if(norm() < FinalRange) break; // if arrive stop condition then end

        LineSearchMethod();
    }
}

void CBasicOptimizer::CSolverMethod::ConjugateGradientMethod() {
    int k = 0;

    for(;;) {

        if (k == 0) {
            DifferentialMethod();

            if(norm() < FinalRange) break; // if arrive stop condition then end
        }

        pfX_past = (double*)memcpy(
                       new double[Super->nDim],
                       Super->pfX,
                       Super->nDim*sizeof(double)
                   );

        LineSearchMethod();

        f_g_past = (double*)memcpy(
                       new double[Super->nDim],
                       Super->f_g,
                       Super->nDim*sizeof(double)
                   );

        DifferentialMethod();

        if(norm() < FinalRange) break; // if arrive stop condition then end

        double beta;

        switch(Super->CFormula) {
        case HestenesStiefel:
            beta = HestenesStiefelFormula();
            break;
        case PolakRibiere:
            beta = PolakRibiereFormula();
            break;
        case FletcherReeves:
            beta = FletcherReevesFormula();
            break;
        case Powell:
            beta = PowellFormula();
            break;
        default:
            break;
        }

        for(int i = 0; i<Super->nDim; i++)
            Super->f_g[i] += beta * f_g_past[i];

        if (k<=n) k++;
        else {
            Super->pfX = (double*)memcpy(
                             new double[Super->nDim],
                             pfX_past,
                             Super->nDim*sizeof(double)
                         );
            k = 0;
        }
        delete [] f_g_past;
        delete [] pfX_past;
    }
}

void CBasicOptimizer::CSolverMethod::DifferentialMethod() {
    CDifferentialMethod *DifferentialMethod = new CDifferentialMethod(Super);
    delete Super->f_g;
    switch(Super->Differential) {
    case Forward:
        Super->f_g = DifferentialMethod->ForwardMethod();
        break;
    case Backward:
        Super->f_g = DifferentialMethod->BackwardMethod();
        break;
    case Central:
        Super->f_g = DifferentialMethod->CentralMethod();
        break;
    default:
        break;
    }
    delete DifferentialMethod;
}

void CBasicOptimizer::CSolverMethod::LineSearchMethod() {
    CLineSearch *LineSearch = new CLineSearch(Super);
    LineSearch->Find_Interval();
    switch(Super->Line) {
    case GoldenSection:
        Super->pfX = LineSearch->GoldenSelection();
        break;
    case Fibonacci:
        Super->pfX = LineSearch->Fibonacci();
        break;
    default:
        break;
    }
    delete LineSearch;
}

// calculate normal of gradiant
double CBasicOptimizer::CSolverMethod::norm() {
    double tmp = 0;
    for(int i = 0; i<Super->nDim; i++)
        tmp += Super->f_g[i]*Super->f_g[i];
    return sqrt(tmp);
}

double CBasicOptimizer::CSolverMethod::HestenesStiefelFormula() {
    double tmp1 = 0;
    double tmp2 = 0;
    for(int i = 0; i<Super->nDim; i++) {
        tmp1 += Super->f_g[i] * (Super->f_g[i] - f_g_past[i]);
        tmp2 += f_g_past[i] * (Super->f_g[i] - f_g_past[i]);
    }
    return tmp1/tmp2;
}

double CBasicOptimizer::CSolverMethod::PolakRibiereFormula() {
    double tmp1 = 0;
    double tmp2 = 0;
    for(int i = 0; i<Super->nDim; i++) {
        tmp1 += Super->f_g[i] * (Super->f_g[i] - f_g_past[i]);
        tmp2 += f_g_past[i] * f_g_past[i];
    }
    return tmp1/tmp2;
}

double CBasicOptimizer::CSolverMethod::FletcherReevesFormula() {
    double tmp1 = 0;
    double tmp2 = 0;
    for(int i = 0; i<Super->nDim; i++) {
        tmp1 += Super->f_g[i] * Super->f_g[i];
        tmp2 += f_g_past[i] * f_g_past[i];
    }
    return tmp1/tmp2;
}

double CBasicOptimizer::CSolverMethod::PowellFormula() {
    double tmp1 = 0;
    double tmp2 = 0;
    for(int i = 0; i<Super->nDim; i++) {
        tmp1 += Super->f_g[i]*(Super->f_g[i]-f_g_past[i]);
        tmp2 += f_g_past[i] * f_g_past[i];
    }
    return (tmp1/tmp2)<0 ? 0 : (tmp1/tmp2);
}
