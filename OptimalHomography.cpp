/*
 *	OptimalHomography.c
 *
 *	Description:
 *
 *
 *
 * 	History:
 *	 	Author			Date			Modify Reason
 *		----------------------------------------------------------------
 *		Chih-En Wu		2012/12/28		File Creation
 *
 *
 */

#include "OptimalHomography.h"

void PrintMatrix(Mat *Matrix) {
    for(int i=0; i<Matrix->rows; i++) {
        for(int j=0; j<Matrix->cols; j++)
            cout<<Matrix->at<Vec3b>(i,j)[0]<<" ";
        cout<<endl;
    }
}

double avg_of_col(Mat *Matrix, int col) {
    double tmp = 0;
    for(int i = 0; i<Matrix->rows; i++)
        tmp += Matrix->at<Vec3b>(i, col)[0];
    return tmp/Matrix->rows;
}

double calulate_s(Mat *Matrix) {
    double x_m = avg_of_col(Matrix, 0);
    double y_m = avg_of_col(Matrix, 1);

    double tmp = 0;
    for(int i = 0; i<Matrix->rows; i++) {
        double tmp1 = Matrix->at<Vec3b>(i, 0)[0];
        double tmp2 = Matrix->at<Vec3b>(i, 1)[0];
        tmp += sqrt(((tmp1-x_m)*(tmp1-x_m))+((tmp2-y_m)*(tmp2-y_m)));
    }

    return 1.414213562373095/(tmp/Matrix->rows);
}

Mat *Get_T(Mat *Matrix) {
    // Calculate average of X, Y
    double x_m = avg_of_col(Matrix, 0);
    double y_m = avg_of_col(Matrix, 1);
    // Calculate s of graph
    double tmp = calulate_s(Matrix);

    Mat *CM_Tmp = new Mat(3, 3, CV_64FC1);

    CM_Tmp->at<Vec3b>(1, 0)[0] = 0;
    CM_Tmp->at<Vec3b>(2, 0)[0] = 0;
    CM_Tmp->at<Vec3b>(2, 1)[0] = 0;
    CM_Tmp->at<Vec3b>(0, 1)[0] = 0;
    CM_Tmp->at<Vec3b>(2, 2)[0] = 1;
    CM_Tmp->at<Vec3b>(0, 0)[0] = tmp;
    CM_Tmp->at<Vec3b>(1, 1)[0] = tmp;
    CM_Tmp->at<Vec3b>(0, 2)[0] = -tmp*x_m;
    CM_Tmp->at<Vec3b>(1, 2)[0] = -tmp*y_m;

    return CM_Tmp;
}

double OH_Cost_Core(Mat **Src, Mat **Dst, Mat *H) {
    double tmp = 0;
    for(int i = 0; i<g_nNumPoints; i++) {
        Mat *CM_Tmp = new Mat(Src[i]->rows, Src[i]->cols, CV_64FC1);
        cvmMul(H, Src[i], CM_Tmp);
        cvmSub(Dst[i], CM_Tmp, CM_Tmp);
        tmp += cvNorm(CM_Tmp, NULL, CV_L2)*cvNorm(CM_Tmp, NULL, CV_L2);
    }
    return tmp/g_nNumPoints;
}

// Create graph matrix
Mat *CM_g_Ref_X = new Mat(g_nNumPoints, 2, CV_64FC1);
Mat *CM_g_x1 = new Mat(g_nNumPoints, 2, CV_64FC1);
Mat *CM_g_x2 = new Mat(g_nNumPoints, 2, CV_64FC1);
Mat *CM_g_x3 = new Mat(g_nNumPoints, 2, CV_64FC1);
Mat *CM_g_x4 = new Mat(g_nNumPoints, 2, CV_64FC1);
Mat *CM_g_x5 = new Mat(g_nNumPoints, 2, CV_64FC1);

// Create H matrix
Mat *CM_g_x1_H = new Mat(3, 3, CV_64FC1);
Mat *CM_g_x2_H = new Mat(3, 3, CV_64FC1);
Mat *CM_g_x3_H = new Mat(3, 3, CV_64FC1);
Mat *CM_g_x4_H = new Mat(3, 3, CV_64FC1);
Mat *CM_g_x5_H = new Mat(3, 3, CV_64FC1);

// Create normalize matrix
Mat **g_Ref_X_Normal = new Mat*[g_nNumPoints];
Mat **g_x1_Normal = new Mat*[g_nNumPoints];
Mat **g_x2_Normal = new Mat*[g_nNumPoints];
Mat **g_x3_Normal = new Mat*[g_nNumPoints];
Mat **g_x4_Normal = new Mat*[g_nNumPoints];
Mat **g_x5_Normal = new Mat*[g_nNumPoints];

// Create T matrix
Mat *CM_g_Ref_X_T;
Mat *CM_g_x1_T;
Mat *CM_g_x2_T;
Mat *CM_g_x3_T;
Mat *CM_g_x4_T;
Mat *CM_g_x5_T;

void OH_Initial() {
    // Set graph data into matrix
    cvSetData(CM_g_Ref_X, g_Ref_X, CM_g_Ref_X->step);
    cvSetData(CM_g_x1, g_x1, CM_g_Ref_X->step);
    cvSetData(CM_g_x2, g_x2, CM_g_Ref_X->step);
    cvSetData(CM_g_x3, g_x3, CM_g_Ref_X->step);
    cvSetData(CM_g_x4, g_x4, CM_g_Ref_X->step);
    cvSetData(CM_g_x5, g_x5, CM_g_Ref_X->step);

    // Calculate T matrix
    CM_g_Ref_X_T = Get_T(CM_g_Ref_X);
    CM_g_x1_T = Get_T(CM_g_x1);
    CM_g_x2_T = Get_T(CM_g_x2);
    CM_g_x3_T = Get_T(CM_g_x3);
    CM_g_x4_T = Get_T(CM_g_x4);
    CM_g_x5_T = Get_T(CM_g_x5);

    // Calcualte normalize of graph points
    for(int i = 0; i<g_nNumPoints; i++) {
        g_Ref_X_Normal[i] = new Mat(3, 1, CV_64FC1);
        g_Ref_X_Normal[i]->at<Vec3b>(0, 0)[0] =  CM_g_Ref_X->at<Vec3b>(i, 0)[0];
        cvSet2D(g_Ref_X_Normal[i], 1, 0, cvGet2D(CM_g_Ref_X,i,1));
        cvSet2D(g_Ref_X_Normal[i], 2, 0, cvScalar(1));
        cvmMul(CM_g_Ref_X_T, g_Ref_X_Normal[i], g_Ref_X_Normal[i]);
    }

    for(int i = 0; i<g_nNumPoints; i++) {
        g_x1_Normal[i] = new Mat(3, 1, CV_64FC1);
        cvSet2D(g_x1_Normal[i], 0, 0, cvGet2D(CM_g_x1,i,0));
        cvSet2D(g_x1_Normal[i], 1, 0, cvGet2D(CM_g_x1,i,1));
        cvSet2D(g_x1_Normal[i], 2, 0, cvScalar(1));
        cvmMul(CM_g_x1_T, g_x1_Normal[i], g_x1_Normal[i]);
    }

    for(int i = 0; i<g_nNumPoints; i++) {
        g_x2_Normal[i] = new Mat(3, 1, CV_64FC1);
        cvSet2D(g_x2_Normal[i], 0, 0, cvGet2D(CM_g_x2,i,0));
        cvSet2D(g_x2_Normal[i], 1, 0, cvGet2D(CM_g_x2,i,1));
        cvSet2D(g_x2_Normal[i], 2, 0, cvScalar(1));
        cvmMul(CM_g_x2_T, g_x2_Normal[i], g_x2_Normal[i]);
    }

    for(int i = 0; i<g_nNumPoints; i++) {
        g_x3_Normal[i] = new Mat(3, 1, CV_64FC1);
        cvSet2D(g_x3_Normal[i], 0, 0, cvGet2D(CM_g_x3,i,0));
        cvSet2D(g_x3_Normal[i], 1, 0, cvGet2D(CM_g_x3,i,1));
        cvSet2D(g_x3_Normal[i], 2, 0, cvScalar(1));
        cvmMul(CM_g_x3_T, g_x3_Normal[i], g_x3_Normal[i]);
    }

    for(int i = 0; i<g_nNumPoints; i++) {
        g_x4_Normal[i] = new Mat(3, 1, CV_64FC1);
        cvSet2D(g_x4_Normal[i], 0, 0, cvGet2D(CM_g_x4,i,0));
        cvSet2D(g_x4_Normal[i], 1, 0, cvGet2D(CM_g_x4,i,1));
        cvSet2D(g_x4_Normal[i], 2, 0, cvScalar(1));
        cvmMul(CM_g_x4_T, g_x4_Normal[i], g_x4_Normal[i]);
    }

    for(int i = 0; i<g_nNumPoints; i++) {
        g_x5_Normal[i] = new Mat(3, 1, CV_64FC1);
        cvSet2D(g_x5_Normal[i], 0, 0, cvGet2D(CM_g_x5,i,0));
        cvSet2D(g_x5_Normal[i], 1, 0, cvGet2D(CM_g_x5,i,1));
        cvSet2D(g_x5_Normal[i], 2, 0, cvScalar(1));
        cvmMul(CM_g_x5_T, g_x5_Normal[i], g_x5_Normal[i]);
    }

    // Find homography
// 	cvFindHomography(CM_g_x1, CM_g_Ref_X, CM_g_x1_H);
// 	cvFindHomography(CM_g_x2, CM_g_Ref_X, CM_g_x2_H);
// 	cvFindHomography(CM_g_x3, CM_g_Ref_X, CM_g_x3_H);
// 	cvFindHomography(CM_g_x4, CM_g_Ref_X, CM_g_x4_H);
// 	cvFindHomography(CM_g_x5, CM_g_Ref_X, CM_g_x5_H);

    // Set H as unit matrix
    cvSetIdentity(CM_g_x1_H, cvScalar(1));
    cvSetIdentity(CM_g_x2_H, cvScalar(1));
    cvSetIdentity(CM_g_x3_H, cvScalar(1));
    cvSetIdentity(CM_g_x4_H, cvScalar(1));
    cvSetIdentity(CM_g_x5_H, cvScalar(1));

    // Calculate normalize of H matrixs
    cvmMul(CM_g_Ref_X_T, CM_g_x1_H, CM_g_x1_H);
    Mat *CM_g_x1_T_Invert = new Mat(CM_g_x1_T->rows, CM_g_x1_T->cols, CV_64FC1);
    cvInvert(CM_g_x1_T, CM_g_x1_T_Invert, CV_LU);
    cvmMul(CM_g_x1_H, CM_g_x1_T_Invert, CM_g_x1_H);

    cvmMul(CM_g_Ref_X_T, CM_g_x2_H, CM_g_x2_H);
    Mat *CM_g_x2_T_Invert = new Mat(CM_g_x2_T->rows, CM_g_x2_T->cols, CV_64FC1);
    cvInvert(CM_g_x2_T, CM_g_x2_T_Invert, CV_LU);
    cvmMul(CM_g_x2_H, CM_g_x2_T_Invert, CM_g_x2_H);

    cvmMul(CM_g_Ref_X_T, CM_g_x3_H, CM_g_x3_H);
    Mat *CM_g_x3_T_Invert = new Mat(CM_g_x3_T->rows, CM_g_x3_T->cols, CV_64FC1);
    cvInvert(CM_g_x3_T, CM_g_x3_T_Invert, CV_LU);
    cvmMul(CM_g_x3_H, CM_g_x3_T_Invert, CM_g_x3_H);

    cvmMul(CM_g_Ref_X_T, CM_g_x4_H, CM_g_x4_H);
    Mat *CM_g_x4_T_Invert = new Mat(CM_g_x4_T->rows, CM_g_x4_T->cols, CV_64FC1);
    cvInvert(CM_g_x4_T, CM_g_x4_T_Invert, CV_LU);
    cvmMul(CM_g_x4_H, CM_g_x4_T_Invert, CM_g_x4_H);

    cvmMul(CM_g_Ref_X_T, CM_g_x5_H, CM_g_x5_H);
    Mat *CM_g_x5_T_Invert = new Mat(CM_g_x5_T->rows, CM_g_x5_T->cols, CV_64FC1);
    cvInvert(CM_g_x5_T, CM_g_x5_T_Invert, CV_LU);
    cvmMul(CM_g_x5_H, CM_g_x5_T_Invert, CM_g_x5_H);
}
