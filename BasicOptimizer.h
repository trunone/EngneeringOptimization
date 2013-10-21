/*
 *	BasicOptimizer.h
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
 *
 * 	History:
 *	 	Author			Date			Modify Reason
 *		----------------------------------------------------------------
 *		Chi-Yi Tsai		2012/09/05		File Creation
 *
 *
 */

//#define MATLAB_DRAW_2D

#include <cmath>
#include <iostream>
#include <vector>
#if defined(MATLAB_DRAW_2D) || defined(MATLAB_DRAW_3D)
#include "draw2mcl.h"
#endif

using namespace std;

enum	LineSearch {GoldenSection = 1, Fibonacci};
enum	SolverMethod {SteepestDecent = 1, ConjugateGradient, QuasiNewton};
enum	DifferentialMethod {Forward = 1, Backward, Central};
enum	ConjugateFormula {PolakRibiere = 1, HestenesStiefel, FletcherReeves, Powell};
enum	QuasiNewtonFormula {DFP = 1, BFGS};

typedef double (*pfnScaleCostFun)(double *pfX);

class CBasicOptimizer {
public:
    CBasicOptimizer();
    ~CBasicOptimizer();

    // return norm value of the gradient vector of the local minimum point
    double	OptScaleCostFunSolver(pfnScaleCostFun pfun, double *pfX, int nDim);

    // return the local minimum point
    double	OptOneDimensionFunSolver(pfnScaleCostFun pfun, double fX);

    void	SetLineSearchMethod(LineSearch Line);
    void	SetSolverMethod(SolverMethod OptMethod);
    void	SetDifferentialMethod(DifferentialMethod Differential);
    void	SetConjugateFormula(ConjugateFormula Formula);
    void	SetQuasiNewtonFormula(QuasiNewtonFormula Formula);

    LineSearch			GetLineSearchMethod(void);
    SolverMethod		GetSolverMethod(void);
    DifferentialMethod	GetDifferentialMethod(void);
    ConjugateFormula	GetConjugateFormula(void);
    QuasiNewtonFormula	GetQuasiNewtonFormula(void);

protected:

private:
    LineSearch			Line;
    SolverMethod		OptMethod;
    DifferentialMethod	Differential;
    ConjugateFormula	CFormula;
    QuasiNewtonFormula	QNFormula;

    CBasicOptimizer *Self; // the pointer of self

    pfnScaleCostFun pfun; //the cost function

    double*	pfX; // Intial point
    double*	f_g; // Direction vector from gradiant
    int		nDim; // Dimension of fucntion

    // Line Search Method =================================================
    class CLineSearch {
    public:
        CLineSearch(CBasicOptimizer*);
        void	Find_Interval(); // Phase 1
        double*	GoldenSelection(); // Phase 2
        double*	Fibonacci(); // Phase 2
        double*	Get_fX(double); // map one dimension to multi-dimesion

#ifdef MATLAB_DRAW_2D
        void	DrawFunction();
        void	DrawInterval(double, double);
#endif

    private:
        // Constants ======================================================
        const double Delta;
        const double IncrementCoeff;
        const double Rho;
        const double FinalRange; // Stop Condition
        const double Epsilon;

        // Data ===========================================================
        CBasicOptimizer *Super; // the pointer of parent class

        double	UpperBound;
        double	LowerBound;
        double	Interval;
    };

    // Differential Method ================================================
    class CDifferentialMethod {
    public:
        CDifferentialMethod(CBasicOptimizer*);
        double* ForwardMethod();
        double* BackwardMethod();
        double* CentralMethod();

    private:
        // Constants ======================================================
        const double Epsilon;

        // Data ===========================================================
        CBasicOptimizer *Super; // the pointer of parent class
    };

    class CSolverMethod {
    public:
        CSolverMethod(CBasicOptimizer*);
        void	SteepestDecentMethod();
        void	ConjugateGradientMethod();
        void	QuasiNewtonMethod();

    private:
        // Constants ======================================================
        const double	FinalRange; // Stop Condition

        // Data ===========================================================
        CBasicOptimizer *Super; // the pointer of parent class

        int		n;

        double *f_g_past;
        double *pfX_past;

        void	DifferentialMethod();
        void	LineSearchMethod();
        double	norm();

        double PolakRibiereFormula();
        double HestenesStiefelFormula();
        double FletcherReevesFormula();
        double PowellFormula();
    };
};
