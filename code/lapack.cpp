/// 
/// @file lapack.h
/// @Synopsis Definition of the wrapper function for CLAPACK diag
/// @author Ivan Gonzalez
/// @date 2007-02-06
/// 
/// Definitions for LAPACK related function declared in lapack.h
/// This is a quote of a post in comp.lang.c++:
/// (search for "preventing multiple definition")  
///
///-----------------------------------------------------------------
/// In fact, my last answer is wrong. What really happened here
/// is the following:
/// You can't DEFINE (not declare) a function in a file that is
/// included by other files unless this function is inlined
/// (explicitly or inside a class). What you should do in this case
/// is the following:
/// 1) Put the definition of the function f() in a separate file
///   and compile it to an object file as you did with a.cpp an
///   b.cpp.
/// 2) Create a header file with the declaration of the function
///   f() and include it by a.cpp, b.cpp and foo.c
/// - OR -
/// 1) Inline the function f() in your file 'c'. 
///-----------------------------------------------------------------
///
#include"lapack.h"
/*****************************************************************************/
///
/// Function to take a real my_Matrix  and diag it with lapack
///
void diagWithLapack_R(double *a, vector<double>& EigenVals, int &rows_, int &cols_)
{
      //int rows_=DMpart.size(); //assuming it is square here
      //int cols_=DMpart.size();
      ///
      /// CLAPACK function to diagonalize an Hermitian matrix
      ///      
      char jobz='V';
      char uplo='U';
      int n=cols_;
      int lda=rows_;
      int info;

      int elems=rows_*cols_;
      //
      // DMPart is hermitian so we use it to get the fortran array
      //
      //double* a;
      //a=( DMpart.transpose(secondDim,firstDim) ).data();
      // 
      // Output 
      //
      for (int j=0; j<elems; j++) cout<<j<<" "<<" "<<a[j]<<endl;
      //
      // Prepare to do a workspace query 
      //
      int lwork=-1;
      int lwork_=1;
      double *work_query= new double [lwork_];
      // 
      // Workspace query
      //
      double w[n];
      int info_=dsyev_(&jobz, &uplo, &n, a, &lda, w, work_query, 
          &lwork, &info);
      // 
      // Get sizes of the workspace and reallocate
      //
      lwork_=(int)fabs((work_query[0]));
      //cout<<"after query\n"<<" work_query[0] "<<lwork<<endl;
      
      delete[] work_query;
      double *work= new double [lwork_];
      // 
      // Call to dsyev_
      //
      info_=dsyev_(&jobz, &uplo, &n, a, &lda, w, work, 
          &lwork_, &info);
      //
      // Free all
      //
      delete[] work;
      //
      // Transpose the DM part (a is row ordered)
      //
      //DMpart.transposeSelf(secondDim,firstDim);
      // 
      // Output 
      //
      for(int i=0; i<n; i++) EigenVals.push_back(w[i]);
      // 
      // Output 
      //
      for (int j=0; j<elems; j++) cout<<j<<" "<<" "<<a[j]<<endl;
      for(int i=0; i<EigenVals.size(); i++) cout<<EigenVals[i]<<endl;
      
      cout<<"Info "<<info<<endl;
};
