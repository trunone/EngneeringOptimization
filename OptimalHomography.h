/*
 *	OptimalHomography.h
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
#include <iostream>

#include "opencv/cv.h"
#include "opencv/highgui.h"

#include "MatchingPoints.h"

using namespace std;
using namespace cv;

extern Mat *CM_g_Ref_X;
extern Mat *CM_g_x1;
extern Mat *CM_g_x2;
extern Mat *CM_g_x3;
extern Mat *CM_g_x4;
extern Mat *CM_g_x5;

extern Mat *CM_g_Ref_X_T;
extern Mat *CM_g_x1_T;
extern Mat *CM_g_x2_T;
extern Mat *CM_g_x3_T;
extern Mat *CM_g_x4_T;
extern Mat *CM_g_x5_T;

extern Mat *CM_g_x1_H;
extern Mat *CM_g_x2_H;
extern Mat *CM_g_x3_H;
extern Mat *CM_g_x4_H;
extern Mat *CM_g_x5_H;

extern Mat **g_Ref_X_Normal;
extern Mat **g_x1_Normal;
extern Mat **g_x2_Normal;
extern Mat **g_x3_Normal;
extern Mat **g_x4_Normal;
extern Mat **g_x5_Normal;

extern Mat *CM_g_Ref_X_T;
extern Mat *CM_g_x1_T;
extern Mat *CM_g_x2_T;
extern Mat *CM_g_x3_T;
extern Mat *CM_g_x4_T;
extern Mat *CM_g_x5_T;

extern "C" {

//void PrintMatrix(CvMat *Matrix);
    void PrintMatrix(Mat *Matrix);

//double avg_of_col(CvMat *, int);
    double avg_of_col(Mat *, int);

//double calulate_s(CvMat *);
    double calulate_s(Mat *);

//CvMat *Get_T(CvMat *);
    Mat *Get_T(Mat *);

    void OH_Initial();

//double OH_Cost_Core(CvMat **, CvMat **, CvMat *);
    double OH_Cost_Core(Mat **, Mat **, Mat *);

};