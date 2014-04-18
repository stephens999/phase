<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0039)http://math.nist.gov/tnt/src/vecinit.cc -->
<HTML><HEAD>
<META http-equiv=Content-Type content="text/html; charset=windows-1252">
<META content="MSHTML 5.50.4134.600" name=GENERATOR></HEAD>
<BODY><XMP>
#include "tnt/tnt.h"
#include "tnt/vec.h"

using namespace std;
using namespace TNT;

int main()
{
    Vector<double> A;                                   /* 1 */

    Vector<double> B(4);                                /* 2 */

    Vector<double> C(4, 1.1);                           /* 3 */

    double t[4] = {1.1, 2.2, 3.3, 4.4};                 /* 4 */
    Vector<double> D(4, &t[0]);

    Vector<double> E(D);                                /* 5 */

    Vector<double> F(4, "1.3  2.7  6.3  1.4");          /* 6 */

    
    cout << A << endl;
    cout << B << endl;
    cout << C << endl;
    cout << D << endl;
    cout << E << endl;
    cout << F << endl;

	return 0;
}
</XMP></BODY></HTML>
