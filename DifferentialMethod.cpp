#include "BasicOptimizer.h"

CBasicOptimizer::CDifferentialMethod::CDifferentialMethod(CBasicOptimizer* p) : Super(p),
    Epsilon(0.001) {}

double* CBasicOptimizer::CDifferentialMethod::BackwardMethod() {
    double *fY = new double[Super->nDim]();
    double	fY_Tmp = Super->pfun(Super->pfX); // Calculate value of original point
    for(int i=0; i<Super->nDim; i++) {
        // Copy from original point
        double *fX_Tmp = (double*)memcpy(
                             new double[Super->nDim],
                             Super->pfX,
                             Super->nDim*sizeof(double)
                         );
        fX_Tmp[i] += fX_Tmp[i] * this->Epsilon;

        fY[i] = (Super->pfun(fX_Tmp) - fY_Tmp) / ((Super->pfX[i] * this->Epsilon)+FLT_MIN);

        delete [] fX_Tmp;
    }
    return fY; // return direction vector
}

double* CBasicOptimizer::CDifferentialMethod::ForwardMethod() {
    double *fY = new double[Super->nDim]();
    double	fY_Tmp = Super->pfun(Super->pfX);
    for(int i=0; i<Super->nDim; i++) {
        // Copy from original point
        double *fX_Tmp = (double*)memcpy(
                             new double[Super->nDim],
                             Super->pfX,
                             Super->nDim*sizeof(double)
                         );
        fX_Tmp[i] -= fX_Tmp[i] * this->Epsilon;

        fY[i] = (fY_Tmp - Super->pfun(fX_Tmp)) / ((Super->pfX[i] * this->Epsilon)+FLT_MIN);

        delete [] fX_Tmp;
    }
    return fY; // return direction vector
}

double* CBasicOptimizer::CDifferentialMethod::CentralMethod() {
    double *fY = new double[Super->nDim]();
    for(int i=0; i<Super->nDim; i++) {
        // Copy from original point
        double *fX_Tmp1 = (double*)memcpy(
                              new double[Super->nDim],
                              Super->pfX,
                              Super->nDim*sizeof(double)
                          );
        fX_Tmp1[i] += fX_Tmp1[i] * (this->Epsilon/2);

        // Copy from original point
        double *fX_Tmp2 = (double*)memcpy(
                              new double[Super->nDim],
                              Super->pfX,
                              Super->nDim*sizeof(double)
                          );
        fX_Tmp2[i] -= fX_Tmp2[i] * (this->Epsilon/2);

        fY[i] = (Super->pfun(fX_Tmp1) - Super->pfun(fX_Tmp2)) / ((Super->pfX[i] * this->Epsilon)+FLT_MIN);

        delete [] fX_Tmp1;
        delete [] fX_Tmp2;
    }
    return fY; // return direction vector
}
