#include <cassert>
#include <complex>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor-fftw/basic.hpp>  // rfft, irfft
#include <xtensor-fftw/helper.hpp> // rfftscale
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp> // xt::arange
#include <xtensor/xio.hpp>
#include <xtensor/xmath.hpp> // xt::sin, cos

using Arr = xt::xarray<double, xt::layout_type::row_major>;

/**
 * @brief Spectral factorization
 *
 *    Spectral factorization using Kolmogorov 1939 approach.
 *      (code follows pp. 232-233, Signal Analysis, by A. Papoulis)
 *
 *      Computes the minimum-phase impulse response which satisfies
 *      given auto-correlation.
 *
 *      Input:
 *        r: top-half of the auto-correlation coefficients
 *           starts from 0th element to end of the auto-corelation
 *           should be passed in as a column vector
 *      Output
 *        h: impulse response that gives the desired auto-correlation
 *
 * @return auto
 */
// auto spectral_fact(const Arr& r)
// {
//     // length of the impulse response sequence
//     const auto n = r.size();

//     // over-sampling factor
//     const auto mult_factor = 100;  // should have mult_factor*(n) >> n
//     const auto m = mult_factor * n;
//     const auto PI = std::acos(-1);

//     // computation method:
//     // H(exp(jTw)) = alpha(w) + j*phi(w)
//     // where alpha(w) = 1/2*ln(R(w)) and phi(w) = Hilbert_trans(alpha(w))

//     // compute 1/2*ln(R(w))
//     // w = 2*pi*[0:m-1]/m
//     Arr w = xt::linspace<double>(0, 2 * PI, m);
//     // R = [ones(m, 1) 2*cos(kron(w', [1:n-1]))]*r
//     Arr Bn = xt::linalg::outer(w, xt::arange(1, n));
//     Arr An = 2 * xt::cos(Bn);
//     // Arr R = np.hstack((np.ones((m, 1)), An)).dot(r)
//     Arr A = xt::concatenate(xt::xtuple(xt::ones<double>({m, 1}), An), 1);
//     Arr R = A.dot(r);

//     auto alpha = 0.5 * xt::log(xt::abs(R));

//     // find the Hilbert transform
//     auto alphatmp = xt::fftw::rfft(alpha);
//     // alphatmp(floor(m/2)+1: m) = -alphatmp(floor(m/2)+1: m)
//     auto ind = m / 2;
//     auto alphatmp[ind:m] = -alphatmp[ind:m];
//     auto alphatmp[0] = 0;
//     auto alphatmp[ind] = 0;

//     // multiply by i*k
//     std::complex<double> i {0, 1};
//     // auto k = xt::fftw::rfftscale<double>(sin.shape()[0], dx);
//     // xt::xarray<std::complex<double>> temp= xt::eval(i * alphatmp);
//     auto phi = xt::fftw::irfft(xt::eval(i * alphatmp)));

//     // now retrieve the original sampling
//     // index = find(np.reminder([0:m-1], mult_factor) == 0)
//     index = np.arange(m, step=mult_factor)
//     alpha1 = alpha[index]
//     phi1 = phi[index]

//     // compute the impulse response (inverse Fourier transform)
//     h = np.real(np.fft.ifft(np.exp(alpha1 + 1j * phi1), n))

//     return h
// }

// def inverse_spectral_fact(h):
//     """[summary]

//     Arguments:
//         h ([type]): [description]

//     Returns:
//         [type]: [description]
//     """
//     n = len(h)
//     r = np.zeros(n)
//     for t in range(n):
//         r[t] = h[t:] @ h[:n - t]
//     return r


// if __name__ == "__main__":
//     r = np.random.rand(20)
//     h = spectral_fact(r)
//     print(h)
