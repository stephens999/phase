#include "tnt/tnt.h"
#include "tnt/cmat.h"

using namespace std;
using namespace TNT;

int main()
{
    complex<double> u(0.0, 0.0);
    complex<double> v(1.0, 2.0);

    Matrix<complex<double> > A(5, 7, u);
    Matrix<complex<double> > B(5, 7, v);
    
    A = u;
    B = v;

    Index1D  I(2,4);
    Index1D  J(3,5);

    A(I,J) = B(I,J);

    cout << A << "\n";

	return 0;
}
