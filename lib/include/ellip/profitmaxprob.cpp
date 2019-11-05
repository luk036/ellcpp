#include "profitmaxprob.hpp"
#include <iostream>

/** 
 * Profit Maximization Problem, taken from:
 *
 * [1] Aliabadi, Hossein, and Maziar Salahi. "Robust Geometric Programming 
 * Approach to Profit Maximization with Interval Uncertainty." Computer 
 * Science Journal of Moldova 21.1 (2013): 61.
 */
template<>
profit_max<double>::profit_max() 
: gp_base<double>(), p{20.0}, alpha{0.1}, beta{0.4}, 
  v1{10.0}, v2{35.0}, k{40}, A{40.0}, cd(3) { gp_setup();}

template<>
profit_max<double>::~profit_max() {}

int main()
{  
  profit_max<double> P;
  std::valarray<double> z = {0.0,1.0,1.0};
  z[0] = log(P.obj(z));
  double bf;        // for output
  ellipsoid E{z, 150};
  
  STATUS status = ellipsoid_dc(E, P, z, bf, 1000, 1e-8);
  if (status == FOUND) {
    std::cout << exp(z[0]) << std::endl;
    std::cout << P.obj(z) << std::endl;
  }
}
