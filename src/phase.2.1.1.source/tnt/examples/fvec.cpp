<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0036)http://math.nist.gov/tnt/src/fvec.cc -->
<HTML><HEAD>
<META http-equiv=Content-Type content="text/html; charset=windows-1252">
<META content="MSHTML 5.50.4134.600" name=GENERATOR></HEAD>
<BODY><XMP>
#include <iostream>
#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/fmat.h"

using namespace std;
using namespace TNT;

int main()
{
    const int N=6;
    Vector<double> x(N), y(N);
    Fortran_Matrix<double> A(N,N);

    A(1,1) = 0.1;
    x(2) = 1.1;
    y(3) = 3.2;

    A(1,2) =  x(2) * A(1,1) + y(3);
    
    std::cout << A(1,2) << "\n";

	return 0;
}

</XMP></BODY></HTML>
