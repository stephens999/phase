<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0037)http://math.nist.gov/tnt/src/tmult.cc -->
<HTML><HEAD>
<META http-equiv=Content-Type content="text/html; charset=windows-1252">
<META content="MSHTML 5.50.4134.600" name=GENERATOR></HEAD>
<BODY><XMP>
#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/fmat.h"
#include "tnt/cmat.h"

using namespace std;
using namespace TNT;

int main()
{
    Matrix<double> A;
    Matrix<double> B;   

    cin >> A;
    cin >> B;

    Matrix<double> C1;
    Matrix<double> C2;
    Matrix<double> C3;


    C1= A*B;
    C2 = matmult(A,B);
    matmult(C3, A, B);

    cout << "A*B: " << A*B << endl;;
    cout << "matmult(A,B): " << matmult(A,B) << endl;
    cout << "matmult(C, A, B): " << C3 << endl;

    Fortran_Matrix<double> A2( A.num_cols(), A.num_rows(),  &A(1,1));
    A2 = transpose(A2);

    Fortran_Matrix<double> B2( B.num_cols(), B.num_rows(),  &B(1,1));
    B2 = transpose(B2);


    cout << endl << endl;
    cout << "C A * B: " << A2 * B2 << endl;
    cout << "C matmult(A, B) : " << matmult(A2, B2) << endl;
    
    Fortran_Matrix<double> C4;
    matmult(C4, A2, B2);
    cout << "Fortran matmult(C, A, B) : " << C4 << endl;


	return 0;
}
    
</XMP></BODY></HTML>
