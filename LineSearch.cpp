#include "BasicOptimizer.h"

CBasicOptimizer::CLineSearch::CLineSearch(CBasicOptimizer* p) : Super(p),
    Delta(0.01),
    IncrementCoeff(1.618),
    Rho(0.382),
    FinalRange(0.0001),
    Epsilon(0.05) {}

void CBasicOptimizer::CLineSearch::Find_Interval() {
    double Coeff = IncrementCoeff; // Coeff of each iteration
    double Alpha_l = 0; // Intial lower bound
    double Alpha_m = Alpha_l + Delta; // Intial middle bound
    double Alpha_u = Alpha_m + Coeff * Delta; // Intial upper bound

    double *real_fX = Get_fX(Alpha_l);
    double f_Alpha_l = Super->pfun(real_fX);
    delete [] real_fX;
    real_fX = Get_fX(Alpha_m);
    double f_Alpha_m = Super->pfun(real_fX);
    delete [] real_fX;
    real_fX = Get_fX(Alpha_u);
    double f_Alpha_u = Super->pfun(real_fX);
    delete [] real_fX;

#ifdef MATLAB_DRAW_2D
    DrawFunction();
    matlab_hold_on();
    DrawInterval(Alpha_l, f_Alpha_l);
    DrawInterval(Alpha_m, f_Alpha_m);
    DrawInterval(Alpha_u, f_Alpha_u);
    matlab_hold_off();
    getchar();
#endif

    if(f_Alpha_l < f_Alpha_m) {
        UpperBound = Alpha_m;
        LowerBound = Alpha_l;
    } else {
        while(!(f_Alpha_l>f_Alpha_m && f_Alpha_m<f_Alpha_u)) {
            Alpha_l = Alpha_m;
            Alpha_m = Alpha_u;
            Coeff *= IncrementCoeff;
            Alpha_u = Alpha_l + Coeff * Delta; // New upper bound

            f_Alpha_l = f_Alpha_m;
            f_Alpha_m = f_Alpha_u;
            real_fX = Get_fX(Alpha_u);
            f_Alpha_u = Super->pfun(Get_fX(Alpha_u));
            delete [] real_fX;

#ifdef MATLAB_DRAW_2D
            DrawFunction();
            matlab_hold_on();
            DrawInterval(Alpha_l, f_Alpha_l);
            DrawInterval(Alpha_m, f_Alpha_m);
            DrawInterval(Alpha_u, f_Alpha_u);
            matlab_hold_off();
            getchar();
#endif
        }
        UpperBound = Alpha_u;
        LowerBound = Alpha_l;
    }
    Interval = UpperBound - LowerBound;
}

double* CBasicOptimizer::CLineSearch::GoldenSelection() {
    double Alpha_a = LowerBound + Rho * Interval; // New lower bound
    double Alpha_b = UpperBound - Rho * Interval; // New upper bound

    double *real_fX = Get_fX(Alpha_a);
    double f_Alpha_a = Super->pfun(real_fX);
    delete [] real_fX;
    real_fX = Get_fX(Alpha_b);
    double f_Alpha_b = Super->pfun(real_fX);
    delete [] real_fX;

#ifdef MATLAB_DRAW_2D
    DrawFunction();
    matlab_hold_on();
    DrawInterval(Alpha_a, f_Alpha_a);
    DrawInterval(Alpha_b, f_Alpha_b);
    matlab_hold_off();
    getchar();
#endif

    while (Interval > FinalRange) {
        if(f_Alpha_a < f_Alpha_b) {
            UpperBound = Alpha_b;
            Interval = UpperBound - LowerBound;

            Alpha_b = Alpha_a;
            Alpha_a = LowerBound + Rho * Interval; // renew lower bound

            f_Alpha_b = f_Alpha_a;
            real_fX = Get_fX(Alpha_a);
            f_Alpha_a = Super->pfun(real_fX);
            delete [] real_fX;
        } else if(f_Alpha_a > f_Alpha_b) {
            LowerBound = Alpha_a;
            Interval = UpperBound - LowerBound;

            Alpha_a = Alpha_b;
            Alpha_b = UpperBound - Rho * Interval; // renew upper bound

            f_Alpha_a = f_Alpha_b;
            real_fX = Get_fX(Alpha_b);
            f_Alpha_b = Super->pfun(real_fX);
            delete [] real_fX;
        } else {
            LowerBound = Alpha_a;
            UpperBound = Alpha_b;
            Interval = UpperBound - LowerBound;

            Alpha_a = LowerBound + Rho * Interval; // renew lower bound
            Alpha_b = UpperBound - Rho * Interval; // renew upper bound

            real_fX = Get_fX(Alpha_a);
            f_Alpha_a = Super->pfun(real_fX);
            delete [] real_fX;
            real_fX = Get_fX(Alpha_b);
            f_Alpha_b = Super->pfun(real_fX);
            delete [] real_fX;
        }
#ifdef MATLAB_DRAW_2D
        DrawFunction();
        matlab_hold_on();
        DrawInterval(Alpha_a, f_Alpha_a);
        DrawInterval(Alpha_b, f_Alpha_b);
        matlab_hold_off();
        getchar();
#endif
    }
    return Get_fX((LowerBound + UpperBound)/2.0);
}

double* CBasicOptimizer::CLineSearch::Fibonacci() {
    vector<unsigned int> vecFib;
    vecFib.push_back(1);
    vecFib.push_back(1);

    double Fn = (1.0 - 2.0 * Epsilon)/(FinalRange/Interval);

    unsigned int N = 1;

    do {
        N++;
        vecFib.push_back(vecFib[N-1]+vecFib[N-2]); // Generate fibonacci squence
    }

    while(Fn >= vecFib[N]);

    double Rho_Fib = 1.0 - ((double)vecFib[N-1]/vecFib[N]);

    double Alpha_a = LowerBound + Rho_Fib * Interval; // New lower bound
    double Alpha_b = UpperBound - Rho_Fib * Interval; // New upper bound

    double f_Alpha_a = Super->pfun(Get_fX(Alpha_a));
    double f_Alpha_b = Super->pfun(Get_fX(Alpha_b));

#ifdef MATLAB_DRAW_2D
    DrawFunction();
    matlab_hold_on();
    DrawInterval(Alpha_a, f_Alpha_a);
    DrawInterval(Alpha_b, f_Alpha_b);
    matlab_hold_off();
    getchar();
#endif

    for(unsigned int i = N; i>1; i--) {
        if(i == 2) {
            Rho_Fib = 0.5 - Epsilon;
        } else {
            Rho_Fib = 1.0 - ((double)vecFib[i-1]/vecFib[i]);
        }

        if(f_Alpha_a < f_Alpha_b) {
            UpperBound = Alpha_b;
            Interval = UpperBound - LowerBound;

            Alpha_b = Alpha_a;
            Alpha_a = LowerBound + Rho_Fib * Interval;

            f_Alpha_b = f_Alpha_a;
            f_Alpha_a = Super->pfun(Get_fX(Alpha_a));
        } else if(f_Alpha_a > f_Alpha_b) {
            LowerBound = Alpha_a;
            Interval = UpperBound - LowerBound;

            Alpha_a = Alpha_b;
            Alpha_b = UpperBound - Rho_Fib * Interval;

            f_Alpha_a = f_Alpha_b;
            f_Alpha_b = Super->pfun(Get_fX(Alpha_b));
        } else {
            LowerBound = Alpha_a;
            UpperBound = Alpha_b;
            Interval = UpperBound - LowerBound;

            Alpha_a = LowerBound + Rho_Fib * Interval;
            Alpha_b = UpperBound - Rho_Fib * Interval;

            f_Alpha_a = Super->pfun(Get_fX(Alpha_a));
            f_Alpha_b = Super->pfun(Get_fX(Alpha_b));
        }
#ifdef MATLAB_DRAW_2D
        DrawFunction();
        matlab_hold_on();
        DrawInterval(Alpha_a, f_Alpha_a);
        DrawInterval(Alpha_b, f_Alpha_b);
        matlab_hold_off();
        getchar();
#endif
    }
    return Get_fX((LowerBound + UpperBound)/2.0);
}

// Map one dimension to multi-dimension
double* CBasicOptimizer::CLineSearch::Get_fX(double delta) {
    // Copy from initial point
    double *fX_Tmp = (double*)memcpy(
                         new double[Super->nDim],
                         Super->pfX,
                         Super->nDim*sizeof(double)
                     );
    // Calculate multi-dimension point
    for(int i = 0; i<Super->nDim; i++) {
        fX_Tmp[i] += delta * (-Super->f_g[i]);
    }
    return fX_Tmp;
}

#ifdef MATLAB_DRAW_2D
void CBasicOptimizer::CLineSearch::DrawInterval(double x, double y) {
    double point[] = {x, y};
    mwArray mwB(2, 1, mxDOUBLE_CLASS);
    mwB.SetData(point, 2);
    matlab_plot_point(mwB);
}

void CBasicOptimizer::CLineSearch::DrawFunction() {
    // Draw fitness function by MATLAB ====================================
    double data[400];
    double step = 0;

    for(unsigned int i = 0; i<400; i+=2, step+=0.01) {
        data[i] = 0 + step;
        data[i+1] = Super->pfun(&data[i]);
    }
    mwArray mwA(2, 200, mxDOUBLE_CLASS);
    mwA.SetData(data, 400);
    matlab_plot(mwA);
}
#endif
