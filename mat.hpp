#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_MAT_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_MAT_HPP 1

#include <valarray>

// Symmetric Matrix
template <typename T> class Mat : public std::valarray<std::valarray<T>> {

  using Vec = std::valarray<T>;

public:
  Mat(const Vec &v, size_t n) : std::valarray<Vec>(v, n) {}
  explicit Mat(size_t n) : std::valarray<Vec>(Vec(n), n) {}

  Vec operator*(const Vec &v) const {
    auto n = v.size(); // depend on the size of v
    Vec res(n);
    for (auto i = 0u; i < n; ++i) {
      res[i] = ((*this)[i] * v).sum();
    }
    return res;
  }

  Mat operator*(const T &x) const {
    auto n = this->size(); // depend on the size of v
    Mat res(n);
    for (auto i = 0u; i < n; ++i) {
      res[i] = (*this)[i] * x;
    }
    return res;
  }

  Mat &operator+=(const Mat &B) {
    const auto n = this->size(); // depend on the size of v
    for (auto i = 0u; i < n; ++i) {
      (*this)[i] += B[i];
    }
    return *this;
  }
};

#endif