#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP 1

//#include "mat.hpp"
#include "chol_ext.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>
namespace bnu = boost::numeric::ublas;

class lmi_oracle {
  using Mat = bnu::symmetric_matrix<double, bnu::upper>;
  using Vec = bnu::vector<double>;
  using Arr = bnu::vector<Mat>;

private:
  const Arr &_F;

public:
  explicit explicit lmi_oracle(const Arr &F) : _F{F} {}

  auto chk_mtx(Mat &A, const Vec &x) const {
    auto sub_mat = [](const Mat &M, size_t p) {
      return bnu::project(M, bnu::range(0, p), bnu::range(0, p));
    };

    auto n = x.size();
    auto g = Vec(n);
    auto fj = -1.0;
    for (auto i = 0u; i < n; ++i) {
      A += this->_F[i] * x[i];
    }
    chol_ext Q(A);
    if (Q.is_sd()) {
      return std::make_tuple(g, fj);
    }
    auto v = Q.witness();
    auto p = v.size();
    fj = -bnu::inner_prod(v, bnu::prod(sub_mat(A, p), v));
    for (auto i = 0u; i < n; ++i) {
      g[i] = -bnu::inner_prod(v, bnu::prod(sub_mat(this->_F[i], p), v));
    }
    return std::make_tuple(g, fj);
  }

  auto chk_spd_t(const Vec &x, double t) const {
    auto A = this->_F[x.size()];
    auto m = A.size1();
    for (auto i = 0u; i < m; ++i) {
      A(i, i) += t;
    }
    return this->chk_mtx(A, x);
  }

  auto chk_spd(const Vec &x) const {
    auto A = this->_F[x.size()];
    return this->chk_mtx(A, x);
  }

  auto operator()(const Vec &x, double t) const {
    auto [g, fj] = this->chk_spd_t(x, t);
    if (fj < 0.0) {
      t -= 1.0;
    }
    return std::make_tuple(g, fj, t);
  }
};

#endif