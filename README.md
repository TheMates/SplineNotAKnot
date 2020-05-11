# SplineNotAKnot
C++ cubic spline library with not-a-knot end conditions and extrapolation


Piecewise cubic polynomial interpolation with not-a-knot boundray condition. Imitates MATLAB **spline()** function.

Utilizes Eigen library.
https://gitlab.com/libeigen/eigen

### Example
```C++
#include "SplineNaK.h"
...
void main()
{
    std::vector<double> x{0., 1., 3., 4.};
    std::vector<dobule> y{0., 0., 2., 2.};
    
    SplineNak::Spline s;
    s.setPoints(x,y);
    
    double ys = s(1.5);     //get interpolated value
    ys = s(5.2);            //extrapolated value
    
    ...
}
```
