#include "ellipsoid.hpp"
#include <iostream>

/** 
 * Problem 1: 
 */
struct micp1
{
  using Vec = std::valarray<double>;
  
  Info4EM<Vec> operator()(const Vec& x0)
  {
    Vec x = x0;
    x[0] = round(x0[0]); // nearest x0[0]
    Vec a(x0.size());
    double f;

    /* C++98 style:
    double ga1[] = {1.1, 1.0};
    g = Vec(ga1, 2);
    f = (g*x).sum() + 3.3;
    if (f > 0.0) {
      Info4Ellipsoid<Vec> res = {false, x, g, f};
      return res;
    }
    */
	
    a = {1.1, 1.0};
    f = (a*x).sum() + 3.3;
    if (f > 0.0) return {false, a, f, x};

    a = {-3.0, 2.0};
    f = (a*x).sum() + 5.5;
    if (f > 0.0) return {false, a, f, x};

    a = {1.1, 2,2};
    Vec z = x - a; 
    return { true, 2.0*z, (z*z).sum(), x };
  }
};

int main()
{
  /* C++98 styly:
  typedef std::valarray<double> Vec;
  */
  using Vec = std::valarray<double>;
  
  micp1 P;
  Vec x = {0.0, 0.0};
  Vec bx(x.size()); // for output
  double bf;        // for output
  ellipsoid E(x, 40.0);
  
  STATUS status = ellipsoid_dc(E, P, bx, bf, 100, 1e-4);
  if (status == FOUND) {
    std::cout << "(" << bx[0] << ", " << bx[1] << ")" << std::endl;
    std::cout << bf << std::endl;
  }
}
