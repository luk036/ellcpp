#ifndef _POSYNOMIAL_HPP
#define _POSYNOMIAL_HPP

#include "monomial.hpp"
#include <valarray>
#include <vector>
#include <cassert>

/**
 * Reference:
 *  S.P. Boyd, S.J. Kim and L.Vandenberghe and A. Hassibi. 
 *  A Tutorial on Geometric Programming. Available at
 *  http://www.standford.edu/~boyd/gp_tutorial.html 
 */

/** 
 * Posynomial function. A sum of one or more monomials, i.e., a
 * function of the form
 *
 *  f(x) = \sum_{k=1}^{K} c_k x_1^{a_{1k}} x_2^{a_{2k}} ... x_n^{a_{nk}},
 *
 * where $c_k$ > 0, is called a posynomial function or, more simply, a
 * posynomial (with $K$ terms, in the variables $x_1$, ..., $x_n$. The
 * term `posynomial' is meant to suggest a combination of `positive'
 * and `polynomial'.
 *
 * Any monomial is also a posynomial. Posynomials are close under
 * addition, multiplication, and positive scaling. Posynomials can be
 * divided by monomials (with the result also a posynomial): If $f$ is
 * a posynomial and $g$ is a monomial, then $f/g$ is a posynomial.
 * Note that <_Tp> could be <double> or <AAF> 
 */
template <typename _Tp>
class posynomial
{
  using Vec = std::valarray<double>;
  using _Self = posynomial<_Tp>;

public:
  /** Constructor */
  explicit posynomial(size_t n, size_t N) : 
  _M(N, monomial<_Tp>(n)) {}
  
  /** Constructor */
  posynomial(const monomial<_Tp>& m) : _M(1, m) {}
  
  /** Constructor (for AAF -> double) */
  template <typename Up, class Map>
  posynomial(const posynomial<Up>& posyn, const Map& polarity)
    : _M(posyn._M.size(), monomial<_Tp>(posyn._M[0]._a.size()))
    {
      for (size_t i=0; i<_M.size(); ++i) {
	_M[i] = monomial<_Tp>(posyn._M[i], polarity);
      }
    }

  /** Destructor */
  ~posynomial() {}

  /** Add and assign */
  _Self& operator+=(const monomial<_Tp>& m) 
  { 
    _M.push_back(m); 
    return *this; 
  }

  /** Add and assign */
  _Self& operator+=(const _Self& P) 
  { 
    for (size_t i=0; i<P._M.size(); ++i) 
      _M.push_back(P._M[i]); 
    return *this; 
  }

  /** Multiply and assign */
  _Self& operator*=(const monomial<_Tp>& m) 
  { for (size_t i=0; i<_M.size(); ++i) _M[i] *= m; 
    return *this; }

  /** Divide and assign */
  _Self& operator/=(const monomial<_Tp>& m)
  { for (size_t i=0; i<_M.size(); ++i) _M[i] /= m; 
    return *this; }

  /** Multiply and assign */
  _Self& operator*=(const _Tp& c) 
  { 
    for (size_t i=0; i<_M.size(); ++i) _M[i] *= c; 
    return *this; 
  }

  /** Divide and assign */
  _Self& operator/=(const _Tp& c)
  { 
    for (size_t i=0; i<_M.size(); ++i) _M[i] /= c; 
    return *this; 
  }

  /** Multiply. @todo simplify the result. */
  _Self operator*(const _Self& P) const
  { 
     _Self res(_M[0]._a.size(), _M.size() * P._M.size());
     size_t k=0;
     for (size_t i=0; i<_M.size(); ++i) {
       for (size_t j=0; j<P._M.size(); ++j) {
	 res._M[k++] = _M[i] * P._M[j];
       }
     }
     return res;
  }

  /** @todo should combine the following two functions into one,
      and eniminate _p */
  /** Function evaluation of log(f(exp(y))). */
  template <typename _Up>
  _Tp lse(const _Up& y) const
  {
    assert(_M[0]._a.size() == y.size());
    const size_t N = _M.size();
    std::valarray<_Tp> p(_Tp(0.0), N);

    for (size_t i=0; i<N; ++i) 
      p[i] = _M[i].lse(y);

    if (N == 1) // monomial
      return p[0];

    // f <- log(sum_i(exp(y_i)))
    _Tp sum(0.0);
    for (size_t i=0; i<N; ++i) {
      sum += exp(p[i]);
    }

    //_Tp sum = exp(p).sum();
    return log(sum);
  }


  /** Gradient of log(f(exp(y))). Precondition: call f(y) previously */
  template <typename _Up>
  std::valarray<_Tp> log_exp_gradient(const _Up& y) const
  {
    assert(_M[0]._a.size() == y.size());
    const size_t n = y.size();
    const size_t N = _M.size();
    std::valarray<_Tp> g(n);

    if (N == 1) {  // monomial
      const Vec& gt = _M[0].gradient(y);
      for (size_t i=0; i<n; ++i) {
	g[i] = gt[i];
      }
      return g;
    }

    
    // g = Aj' * (exp(yj)./sum(exp(yj)));
    // Note that exp(yj) has been previous calculated in _p during the
    // function evaluation.
    std::valarray<_Tp> z(N);

    _Tp sum(0.0);
    for (size_t i=0; i<N; ++i) {
      z[i] = exp(_M[i].lse(y));
      sum += z[i];
    }

    z /= sum;
    
    for (size_t i=0; i<n; ++i) {
      g[i] = 0.;
      for (size_t l=0; l<N; ++l) {
	g[i] += _M[l]._a[i] * z[l];
      }
    }
 
    return g;
  }

  /** function value and gradient of log(f(exp(y))).  
      Note that for AAF, the two quantities have to be evaluated *at the same time*
	  in order to maintain the noise symbols consistency */
  template <typename _Up>
  _Tp log_exp_fvalue_with_gradient(const _Up& y, std::valarray<_Tp>& g) const
  {
    assert(_M[0]._a.size() == y.size());
    const size_t n = y.size();
    const size_t N = _M.size();
    std::valarray<_Tp> p(N);
    
    for (size_t i=0; i<N; ++i) 
      p[i] = _M[i].lse(y);
    
    if (N == 1) {  // i.e. monomial
      const Vec& gt = _M[0].gradient(y);
      for (size_t i=0; i<n; ++i) {
	g[i] = gt[i];
      }
      return p[0];
    }
    
    // f <- log(sum_i(exp(y_i)))
    _Tp sum(0.0);
    for (size_t i=0; i<N; ++i) {
      p[i] = exp(p[i]);
      sum += p[i];
    }

    // g = Aj' * (exp(yj)./sum(exp(yj)));
    p /= sum;
    
    for (size_t i=0; i<n; ++i) {
      g[i] = 0.0;
      for (size_t l=0; l<N; ++l) {
	g[i] += _M[l]._a[i] * p[l];
      }
    }
    
    return log(sum);
  }


  /** function value and gradient of log(f(exp(y))). */
  template <typename _Up>
  std::valarray<_Tp> lse_gradient(const _Up& y, _Tp& f) const
  {
    assert(_M[0]._a.size() == y.size());
    const size_t n = y.size();
    const size_t N = _M.size();
    std::valarray<_Tp> p(N);
    
    for (size_t i=0; i<N; ++i) 
      p[i] = _M[i].lse(y);
    
    if (N == 1) {  // i.e. monomial
      f = p[0];
      return _M[0].lse_gradient(y);
    }
    
    // f <- log(sum_i(exp(y_i)))
    _Tp sum(0.0);
    for (size_t i=0; i<N; ++i) {
      p[i] = exp(p[i]);
      sum += p[i];
    }

    // g = Aj' * (exp(yj)./sum(exp(yj)));
    p /= sum;
    
    std::valarray<_Tp> g(n);
    for (size_t i=0; i<n; ++i) {
      g[i] = 0.0;
      for (size_t l=0; l<N; ++l) {
	g[i] += _M[l]._a[i] * p[l];
      }
    }
    
    f = log(sum);    
    return g;
  }

public:
  std::vector<monomial<_Tp> > _M;    /**< vector of monomials */

public:
  posynomial(const _Self& Q) : _M(Q._M) {}

public:
  _Self& operator=(const _Self& Q) 
  { 
    _M = Q._M; 
    return *this; 
  }
};

/** Add */
template <typename _Tp>
posynomial<_Tp> operator+(const monomial<_Tp>& m1, const monomial<_Tp>& m2)
{ return posynomial<_Tp>(m1) += m2; }

/** Add */
template <typename _Tp>
posynomial<_Tp> operator+(posynomial<_Tp> p, const monomial<_Tp>& m)
{ 
  p += m;
  return p; 
}

/** Multiply */
template <typename _Tp>
posynomial<_Tp> operator*(posynomial<_Tp> p, const monomial<_Tp>& m)
{ 
  p *= m;
  return p; 
}

/** Divide */
template <typename _Tp>
posynomial<_Tp> operator/(posynomial<_Tp> p, const monomial<_Tp>& m)
{ 
  p /= m;
  return p; 
}

/** Multiply */
template <typename _Tp>
posynomial<_Tp> operator*(posynomial<_Tp> p, const _Tp& c)
{ 
  p *= c;
  return p; 
}

/** Divide */
template <typename _Tp>
posynomial<_Tp> operator/(posynomial<_Tp> p, const _Tp& c)
{ 
  p /= c;
  return p; 
}

#endif
