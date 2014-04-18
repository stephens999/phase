<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0039)http://math.nist.gov/tnt/src/ftransp.cc -->
<HTML><HEAD>
<META http-equiv=Content-Type content="text/html; charset=windows-1252">
<META content="MSHTML 5.50.4134.600" name=GENERATOR></HEAD>
<BODY><XMP>// Examples of using matrix transposes
//
 
#include "tnt/tnt.h"
#include "tnt/fmat.h"
#include "tnt/transv.h"

using namespace std;
using namespace TNT;

int main()
{
    Fortran_Matrix<double> A(2, 4, 
                          " 1  2  0  4 "
                          " 2  0  9  7 ");

    Fortran_Matrix<double> B(4, 4, 
                          " 9  4  1  2"
                          " 1  3  4  9"
                          " 2  8  3  1"
                          " 0  1  0  0");

    
    Vector<double> x(4, " 1  2  3  4 ");

    cout << "A: " << A  << endl;
    cout << "Transpose_view(A): " << Transpose_view(A) << endl;

    cout << "B: " << B << endl;

    cout << "B' * [1 2 3 4]'  : " << Transpose_view(B) * x << endl;

	return 0;
}
    
</XMP></BODY></HTML>
