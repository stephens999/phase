<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0040)http://math.nist.gov/tnt/src/fcomplex.cc -->
<HTML><HEAD>
<META http-equiv=Content-Type content="text/html; charset=windows-1252">
<META content="MSHTML 5.50.4134.600" name=GENERATOR></HEAD>
<BODY><XMP>

#include <complex>
#include "tnt/tnt.h"
#include "tnt/fmat.h"

using namespace std;
using namespace TNT;

int main()
{
    complex<double> u(0.0, 0.0);
    complex<double> v(1.0, 2.0);

    Fortran_Matrix<complex<double> > A(5, 7, u);
    Fortran_Matrix<complex<double> > B(5, 7, v);
    
    A = u;
    B = v;

    Index1D  I(2,4);
    Index1D  J(3,5);

    A(I,J) = B(I,J);

    cout << A << "\n";

	return 0;
}
</XMP></BODY></HTML>
