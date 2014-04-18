<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0037)http://math.nist.gov/tnt/src/tcmat.cc -->
<HTML><HEAD>
<META http-equiv=Content-Type content="text/html; charset=windows-1252">
<META content="MSHTML 5.50.4134.600" name=GENERATOR></HEAD>
<BODY><XMP>



// Illustrate matrix initialization, basic matrix routines, and simple
// matrix I/O.
 
#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/cmat.h"

using namespace TNT;

int main()
{
    Matrix<double> A(2, 4, 
                                   " 1  2  0  4 "
                                   " 2  0  9  7 ");

    Matrix<double> B(4, 2, 
                                   " 9  4 "
                                   " 1  3 "
                                   " 2  8 "
                                   " 0  1 ");

    Matrix<double> C(2, 4, 
                                   " 3  4  8  1 "
                                   " 7  8  0  1 ");

    std::cout << A*B << "\n";

    std::cout << A+C << "\n";

	return 0;

}
    
</XMP></BODY></HTML>
