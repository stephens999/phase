<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0037)http://math.nist.gov/tnt/src/fchol.cc -->
<HTML><HEAD>
<META http-equiv=Content-Type content="text/html; charset=windows-1252">
<META content="MSHTML 5.50.4134.600" name=GENERATOR></HEAD>
<BODY><XMP>

// Test Cholesky module  

#include <iostream>

#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/fmat.h"
#include "tnt/cholesky.h"
#include "tnt/trisolve.h"
#include "tnt/transv.h"         /* transpose views */

using namespace std;
using namespace TNT;

int main()
{
    Fortran_Matrix<double> A;

    cin >> A;                   /* A should be symmetric positive definite */

    Subscript N = A.num_rows();
    assert(N == A.num_cols());

    Vector<double> b(N, 1.0);   // b= [1,1,1,...]
    Fortran_Matrix<double> L(N, N);


    cout << "A: " << A << endl;
    
    if (Cholesky_upper_factorization(A, L) !=0)
    {
        cout << "Cholesky did not work." << endl;
        exit(1);
    }
    

    cout << L << endl;

    // solve Ax =b, as L*L'x =b
    //
    //  let y=L'x, then
    //
    //
    //   solve L y = b;
    //   solve L'x = y;

    Vector<double> y = Lower_triangular_solve(L, b);
    Vector<double> x= Upper_triangular_solve(Transpose_view(L), y);

    cout << "x: " << x << endl;
    cout << "Residual A*x-b: " << A*x-b << endl;

	return 0;
}
</XMP></BODY></HTML>
