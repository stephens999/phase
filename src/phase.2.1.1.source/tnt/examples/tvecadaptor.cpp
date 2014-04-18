<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0043)http://math.nist.gov/tnt/src/tvecadaptor.cc -->
<HTML><HEAD>
<META http-equiv=Content-Type content="text/html; charset=windows-1252">
<META content="MSHTML 5.50.4134.600" name=GENERATOR></HEAD>
<BODY><XMP>

//demonstate matrix/vector declarations, and element access.

#include <iostream>
#include "tnt/tnt.h"
#include <vector>      /* STL vector<> */
#include "tnt/vecadaptor.h"

using namespace std;
using namespace TNT;

typedef TNT::Vector_Adaptor<vector<double> >   Vec;

int main()
{
    const int N=6;
    Vector_Adaptor<vector<double> > x(N, "1.0 2.0 3.0 4.0 5.0 6.0");
    Vec  y(N, 1.0);

    cout << x << "\n";
    cout << y << "\n";

    
    Index1D I(2,4);

    // illustrate 1-based indexing

    x(2) = 1.1;
    y(3) = 3.2;

    x(I) = y(I+1);


    cout << "x: " << x << "\n";
    cout << "y: " << y << "\n";


	return 0;
}

</XMP></BODY></HTML>
