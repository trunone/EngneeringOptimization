/*
 *	BasicOptimization.cpp
 *
 *	Description:
 *		Test for BasicOptimizer
 *
 * 	History:
 *	 	Author			Date			Modify Reason
 *		----------------------------------------------------------------
 *		Chih-En Wu		2012/10/12		File Creation
 *
 *
 */

#include "OptimalHomography.h"

#include "BasicOptimizer.h"

Mat *array2matrix(double *fX) {
    Mat *CM_Tmp = new Mat(3, 3, CV_64FC1);
    cvSet2D(CM_Tmp, 0, 0, cvScalar(fX[0]));
    cvSet2D(CM_Tmp, 0, 1, cvScalar(fX[1]));
    cvSet2D(CM_Tmp, 0, 2, cvScalar(fX[2]));
    cvSet2D(CM_Tmp, 1, 0, cvScalar(fX[3]));
    cvSet2D(CM_Tmp, 1, 1, cvScalar(fX[4]));
    cvSet2D(CM_Tmp, 1, 2, cvScalar(fX[5]));
    cvSet2D(CM_Tmp, 2, 0, cvScalar(fX[6]));
    cvSet2D(CM_Tmp, 2, 1, cvScalar(fX[7]));
    cvSet2D(CM_Tmp, 2, 2, cvScalar(fX[8]));
    return CM_Tmp;
}

double* matrix2array(Mat* Matrix) {
    double *fX = new double[Matrix->cols*Matrix->rows];
    int k = 0;
    for(int i=0; i<Matrix->rows; i++) {
        for(int j=0; j<Matrix->cols; j++) {
            fX[k] = cvGet2D(Matrix, i, j).val[0];
            k++;
        }
    }
    return fX;
}

double OH_Cost1(double *fX) {
    // Set H into OpenCV matrix
    Mat *H = array2matrix(fX);
    return OH_Cost_Core(g_x1_Normal, g_Ref_X_Normal, H);
}

double OH_Cost2(double *fX) {
    // Set H into OpenCV matrix
    Mat *H = array2matrix(fX);
    return OH_Cost_Core(g_x2_Normal, g_Ref_X_Normal, H);
}

double OH_Cost3(double *fX) {
    // Set H into OpenCV matrix
    Mat *H = array2matrix(fX);
    return OH_Cost_Core(g_x3_Normal, g_Ref_X_Normal, H);
}


double f_1d_1(double *fX) {
    return 0.65 - (0.75/(1+ fX[0]*fX[0])) - 0.65*fX[0]*atan(1/fX[0]);
}

double f_1d_2(double *fX) {
    return fX[0]*fX[0]*fX[0]*fX[0] - 14*fX[0]*fX[0]*fX[0] + 60*fX[0]*fX[0] - 70*fX[0];
}

double f_2d_1(double *fX) {
    return fX[0]-fX[1]+2*fX[0]*fX[0]+2*fX[0]*fX[1]+fX[1]*fX[1];
}

double f_2d_2(double *fX) {
    return (100*(fX[1]-fX[0]*fX[0])*(fX[1]-fX[0]*fX[0])+(1-fX[0])*(1-fX[0]))/2;
}

double f_2d_3(double *fX) {
    return (-3*fX[1])/(fX[0]*fX[0]+fX[1]*fX[1]+1);
}

int main() {
    // Initial MATLAB Library =============================================
#if defined(MATLAB_DRAW_2D) || defined(MATLAB_DRAW_3D)
    if(!mclInitializeApplication(NULL, 0)) {
        cout<<"Could not initialize the application!"<<endl;
        return -1;
    } else {
        cout<<"Initialize the application done!"<<endl;
    }

    if(!draw2mclInitialize()) {
        cout<<"Could not initialize MATLAB!"<<endl;
        return -1;
    } else {
        cout<<"Initialize MATLAB done!"<<endl;
    }
#endif
    // Initial optimization ===============================================
    CBasicOptimizer *BasicOptimizer = new CBasicOptimizer;
    BasicOptimizer->SetLineSearchMethod(GoldenSection);
    BasicOptimizer->SetDifferentialMethod(Central);
    BasicOptimizer->SetSolverMethod(ConjugateGradient);
    BasicOptimizer->SetConjugateFormula(PolakRibiere);

    vector<pfnScaleCostFun> vecCost;
    vecCost.push_back(&f_1d_1);
    vecCost.push_back(&f_1d_2);
    vecCost.push_back(&f_2d_1);
    vecCost.push_back(&f_2d_2);
    vecCost.push_back(&f_2d_3);

    OH_Initial();

    // Start optimization =================================================
// 	double min_point = BasicOptimizer->OptOneDimensionFunSolver(vecCost[0], 0);
// 	double min_value = vecCost[0](&min_point);
//
// 	cout<<min_point<<" "<<min_value<<endl<<endl;
//
// 	min_point = BasicOptimizer->OptOneDimensionFunSolver(vecCost[1], 0);
// 	min_value = vecCost[1](&min_point);
//
// 	cout<<min_point<<" "<<min_value<<endl<<endl;

// 	double fX[] = {0.5, -0.5};
// 	double min_value = BasicOptimizer->OptScaleCostFunSolver(vecCost[2], fX, 2);
//
// 	cout<<min_value<<endl<<endl;
//
// 	min_value = BasicOptimizer->OptScaleCostFunSolver(vecCost[3], fX, 2);
//
// 	cout<<min_value<<endl<<endl;
//
// 	min_value = BasicOptimizer->OptScaleCostFunSolver(vecCost[4], fX, 2);
//
// 	cout<<min_value<<endl<<endl;

    double *fX = matrix2array(CM_g_x1_H);

    cout<<OH_Cost1(fX)<<endl;

    double min_value = BasicOptimizer->OptScaleCostFunSolver(OH_Cost1, fX, 9);

    Mat *CM_Tmp = array2matrix(fX);

    Mat *CM_g_Ref_X_T_Invert = new Mat(CM_g_Ref_X_T->rows, CM_g_Ref_X_T->cols, CV_64FC1);
    cvInvert(CM_g_Ref_X_T, CM_g_Ref_X_T_Invert, CV_LU);
    cvmMul(CM_g_Ref_X_T_Invert, CM_Tmp, CM_Tmp);
    cvmMul(CM_Tmp, CM_g_x1_T, CM_Tmp);

    PrintMatrix(CM_Tmp);

    cout<<endl<<min_value<<endl<<endl;

    delete BasicOptimizer;

    // Terminate MATLAB Library ===========================================
#if defined(MATLAB_DRAW_2D) || defined(MATLAB_DRAW_3D)
    draw2mclTerminate();
    mclTerminateApplication();
#endif

    getchar();

    return 0;
}
