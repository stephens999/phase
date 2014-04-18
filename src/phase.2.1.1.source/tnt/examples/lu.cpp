<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0034)http://math.nist.gov/tnt/src/lu.cc -->
<HTML><HEAD>
<META http-equiv=Content-Type content="text/html; charset=windows-1252">
<META content="MSHTML 5.50.4134.600" name=GENERATOR></HEAD>
<BODY><XMP>
// Solve a linear system using LU factorization.
//
// Usage: a.out < matrix.dat
//
// where matrix.dat is an ASCII file consisting of the
// matrix size (M,N) followed by its values.  For example,
//
//  3  2
//  8.1  1.2  4.3
//  1.3  4.3  2.9
//  0.4  1.3  6.1
//


#include <iostream>

#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/fmat.h"
#include "tnt/lu.h"

using namespace std;
using namespace TNT;

int main()
{
    Fortran_Matrix<double> A;

    cin >> A;


    Subscript N = A.dim(1);
    assert(N == A.dim(2));

    Vector<double> b(N, 1.0);   // b= [1,1,1,...]
    Vector<Subscript> index(N);



    cout << "Original Matrix A: " << A << endl;
    
    Fortran_Matrix<double> T(A);
    if (LU_factor(T, index) !=0)
    {
        cout << "LU_factor() failed." << endl;
        exit(1);
    }

    Vector<double> x(b);
    if (LU_solve(T, index, x) != 0)
    {
        cout << "LU_Solve() failed." << endl;
        exit(1);
    }
    cout << "Solution x for Ax=b, where b=[1,1,...] " <<endl;
    cout << " x: " << x << endl;

    cout << "A*x should be the vector [1,1,...] "  <<endl;
    cout     << "residual [A*x - b]: " << matmult(A, x)  - b << endl;
    
	return 0;
}
</XMP></BODY></HTML>
