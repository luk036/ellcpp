#include <ellcpp/oracles/lowpass_oracle.hpp>
#include <ellcpp/utility.hpp>
#include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using ParallelCut = std::tuple<Arr, Arr>;

/*!
 * @brief
 *
 * @param[in] x
 * @param[in] Spsq
 * @return auto
 */
auto lowpass_oracle::operator()(
    const Arr& x, double& Spsq) const -> std::tuple<ParallelCut, bool>
{
    using xt::linalg::dot;

    // 1. nonnegative-real constraint
    // case 1,
    if (x[0] < 0)
    {
        auto g = zeros(x);
        g[0] = -1.;
        auto f = Arr {-x[0]};
        return {{std::move(g), std::move(f)}, false};
    }

    // case 2,
    // 2. passband constraints
    auto N = this->_Ap.shape()[0];
    // for (k in chain(range(i_As, N), range(i_As))) {
    for (size_t i = 0U, k = this->_i_Ap; i != N; ++i, ++k)
    {
        if (k == N)
        {
            k = 0; // round robin
        }
        auto v = dot(xt::view(this->_Ap, k, xt::all()), x)();
        if (v > this->_Upsq)
        {
            // f = v - Upsq;
            Arr g = xt::view(this->_Ap, k, xt::all());
            Arr f {v - this->_Upsq, v - this->_Lpsq};
            this->_i_Ap = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
        if (v < this->_Lpsq)
        {
            // f = Lpsq - v;
            Arr g = -xt::view(this->_Ap, k, xt::all());
            Arr f {-v + this->_Lpsq, -v + this->_Upsq};
            this->_i_Ap = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
    }

    // case 3,
    // 3. stopband constraint
    N = this->_As.shape()[0];
    // Arr w = xt::zeros<double>({N});
    auto fmax = -1.e100; // std::numeric_limits<double>::min()
    size_t imax {0};
    // for (k in chain(range(i_As, N), range(i_As))) {
    for (size_t i = 0U, k = this->_i_As; i != N; ++i, ++k)
    {
        if (k == N)
        {
            k = 0; // round robin
        }
        auto v = dot(xt::view(this->_As, k, xt::all()), x)();
        if (v > Spsq)
        {
            // f = v - Spsq;
            Arr g = xt::view(this->_As, k, xt::all());
            // f = (v - Spsq, v);
            Arr f {v - Spsq, v};
            this->_i_As = k + 1; // k or k+1
            return {{std::move(g), std::move(f)}, false};
        }
        if (v < 0)
        {
            // f = v - Spsq;
            Arr g = -xt::view(this->_As, k, xt::all());
            Arr f {-v, -v + Spsq};
            this->_i_As = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
        if (v > fmax)
        {
            fmax = v;
            imax = k;
        }
    }

    // case 4,
    // 1. nonnegative-real constraint
    N = this->_Anr.shape()[0];
    for (size_t i = 0U, k = this->_i_Anr; i != N; ++i, ++k)
    {
        if (k == N)
        {
            k = 0; // round robin
        }
        auto v = dot(xt::view(this->_Anr, k, xt::all()), x)();
        if (v < 0.)
        {
            Arr f {-v};
            Arr g = -xt::view(this->_Anr, k, xt::all());
            this->_i_Anr = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
    }

    // Begin objective function
    // Spsq, imax = w.max(), w.argmax(); // update best so far Spsq
    Spsq = fmax;
    Arr f {0., fmax}; // ???
    // f = 0
    Arr g = xt::view(this->_As, imax, xt::all());
    return {{std::move(g), std::move(f)}, true};
}
