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
///
/// This was stripped away in May 2013 by Roger Melko for use in extracting the 
/// entanglement entropy with NLCE
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
      //for (int j=0; j<elems; j++) cout<<j<<" "<<" "<<a[j]<<endl;
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
      //for (int j=0; j<elems; j++) cout<<j<<" "<<" "<<a[j]<<endl;
      //for(int i=0; i<EigenVals.size(); i++) cout<<EigenVals[i]<<endl;
      
      //cout<<"Info "<<info<<endl;
};
/*****************************************************************************/
///
/// Function to take a real my_Matrix  and diag it with lapack
///

void svdWithLapack_simple(double *a, vector<double>& EigenVals, int &rows_, int &cols_)
{
      char jobu='N';
      char jobvt='N';
      int m=rows_;
      int n=cols_;
      int lda=rows_;
      int info;
      int SVs = m;
      if(m>n){SVs = n;}

      int elems=rows_*cols_;
      
      // Prepare to do a workspace query 
      int lwork=-1;
      double work_query;

       // the singular values ( of length min(m, n) )
      double *s = new double [SVs];
      double *u; int ldu=m;
      double *vt; int ldvt=n;
   
      // Workspace query
      int info_ = dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu,
          vt, &ldvt, &work_query, &lwork, &info);

      
      // Get sizes of the workspace and reallocate
      lwork = (int)work_query;
      double *work= new double [lwork];
 
      // Call to dgesvd_
      info_ = dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu,
              vt, &ldvt, work, &lwork, &info);

      // Output 
      // is the length n or m???????
      for(int i=0; i<n; i++){ EigenVals.push_back(s[i]*s[i]);} // cout << s[i] << endl;}
    
      // Free all
      delete[] work;
      delete[] s;

};
/*****************************************************************************/
///
/// Function to take a real my_Matrix  and diag it with lapack
///

void svdWithLapack_divide(double *a,  vector<double>& EigenVals, int &rows_, int &cols_)
{
      char jobz='N';
      int m=rows_;
      int n=cols_;
      int lda=rows_;
      int info;
      int SVs = m;
      if(m>n){SVs = n;}
      int iwork[8*n];
      
      int elems=rows_*cols_;
      
      // Prepare to do a workspace query 
      int lwork=-1;
      double work_query;

       // the singular values ( of length min(m, n) )
      double *s = new double [SVs];
      double *u; int ldu=m;
      double *vt; int ldvt=n;
   
      // Workspace query
     dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &work_query, &lwork, iwork, &info);
      
      // Get sizes of the workspace and reallocate
      lwork = (int)work_query;
      double *work= new double [lwork];

      
      // Call to dgesvd_
      dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info );

      /* Check for convergence */
      if( info > 0 ) {
          printf( "The algorithm computing SVD failed to converge.\n" );
      }

      // Output 
      // is the length n or m???????
      for(int i=0; i<n; i++){ EigenVals.push_back(s[i]*s[i]);} // cout << s[i] << endl;}
    
      // Free all
      delete[] work ;
      delete[] s;

};

