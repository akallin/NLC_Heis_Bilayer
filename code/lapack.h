/// 
/// @file lapack.h
/// @Synopsis Definition of the wrapper function for CLAPACK diag
/// @author Ivan Gonzalez
/// @date 2007-02-06
/// 
///
/// Blitz removed, and file stripped down for NLCE: Roger Melko, May 2013
#ifndef LAPACK_H
#define LAPACK_H

//#include<complex>
//#include<vector>
#include"Lanczos_07.h"

/*****************************************************************************/
///
/// Extern declaration of CLAPACK functions
///
extern "C" {
        
//    int zheev_(char *jobz, char *uplo, int *n, complex<double> 
//	    *a, int *lda, double *w, complex<double> *work, int *lwork, 
//	    double *rwork, int *info);
// 
//    int zheevd_(char *jobz, char *uplo, int *n, 
//	    complex<double> *a, int *lda, double *w, complex<double> *work, 
//	    int *lwork, double *rwork, int *lrwork, int *iwork, 
//    	    int *liwork, int *info);

    int dsyev_(char *jobz, char *uplo, int *n, double *a,
	    int *lda, double *w, double *work, int *lwork, int *info);

    int dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, 
        int *lda, double *s, double *u, int *ldu,
        double *vt, int *ldvt, double *work, int *lwork, 
	    int *info);

    int dgesdd_(char *jobz, int *m, int *n, double *a,
        int *lda, double *s, double *u, int *ldu, 
	    double *vt, int *ldvt, double *work, int *lwork, 
	    int *iwork, int *info);

  }
/*****************************************************************************/
///
/// Function to take a real my_Matrix  and diag it with lapack
///
//void diagWithLapack_R(Array<double, 2>& DMpart, vector<double>& EigenVals);
void diagWithLapack_R(double *a,  vector<double>& EigenVals, int &rows_, int &cols_);
void svdWithLapack_simple(double *a,  vector<double>& EigenVals, int &rows_, int &cols_);
void svdWithLapack_divide();

#endif
