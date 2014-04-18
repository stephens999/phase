<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0039)http://math.nist.gov/tnt/src/ctriang.cc -->
<HTML><HEAD>
<META http-equiv=Content-Type content="text/html; charset=windows-1252">
<META content="MSHTML 5.50.4134.600" name=GENERATOR></HEAD>
<BODY><XMP>//  triangular matrices

#include <iostream>
#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/cmat.h"
#include "tnt/triang.h"
#include "tnt/trisolve.h"

using namespace std;
using namespace TNT;

int main()
{
    Matrix<double> A(4, 4, 
                                   " 2  1  1  1 "
                                   " 2  3  1  1 "
                                   " 6  7  4  1 "
                                   " 9  2  3  5 " );

    cout << "A : " << A << endl;

    cout <<"Lower triangular part of A: " <<
        Lower_triangular_view(A) << endl;

    cout << "Bottom left 2x2 corner: " <<
        Lower_triangular_view(A)(Index1D(3,4), Index1D(1,2)) << endl;

    cout << "Unit lower triangular part of A: " <<
        Unit_lower_triangular_view(A) << endl;


    Vector<double> b(4, "1 2 3 4");

    Vector<double> x;
    

    x = linear_solve(Unit_lower_triangular_view(A), b);

    cout << "solution to Unit_lower_triangular_view(A)*x = [1 2 3 4]':  " 
            << x << endl;

    cout << "A*x - b : " << Unit_lower_triangular_view(A)* x - b << endl;


/*
    x = linear_solve(Lower_triangular_view(A), b);

    cout << "solution to Lower_triangular_view(A)*x = [1 2 3 4]':  " 
            << x << endl;

    cout << "A*x - b : " << Lower_triangular_view(A)* x - b << endl;
*/

   	return 0; 
}
    
</XMP></BODY></HTML>
