#include <cmath> // ceil, fabs, log
#include <string>
#include <string_view>

using std::ceil;
using std::fabs;
using std::log2;
using std::pow;
using std::string;

/**
 * @brief Convert to CSD (Canonical Signed Digit) string representation
 *
 * @param num
 * @param places
 * @return string
 */
auto to_csd(double num, int places = 0) noexcept -> string
{
    if (num == 0.)
    {
        return "0";
    }
    auto absnum = fabs(num);
    auto n = absnum < 1. ? 0 : int(ceil(log2(absnum * 1.5)));
    auto csd_str = string {absnum < 1. ? "0" : ""};
    auto limit = pow(2., n) / 3.;
    while (n > -places)
    {
        if (n == 0)
        {
            csd_str += '.';
        }
        n -= 1;
        if (num > limit)
        {
            csd_str += '+';
            num -= 1.5 * limit;
        }
        else if (num < -limit)
        {
            csd_str += '-';
            num += 1.5 * limit;
        }
        else
        {
            csd_str += '0';
        }
        limit /= 2.;
    }
    return csd_str;
}


/**
 * @brief Convert the CSD string to a decimal
 *
 * @param csd_str
 * @return double
 */
auto to_decimal(std::string_view csd_str) noexcept -> double
{
    auto num = 0.0;
    auto loc = 0;
    auto i = 0;
    for (auto c : csd_str)
    {
        switch (c)
        {
            case '0':
                num *= 2.;
                break;
            case '+':
                num *= 2.;
                num += 1.;
                break;
            case '-':
                num *= 2.;
                num -= 1.;
                break;
            case '.':
                loc = i + 1;
                break;
            default:
                break; // ignore
        }
        ++i;
    }
    if (loc != 0)
    {
        num /= std::pow(2, csd_str.size() - loc);
    }
    return num;
}


/**
 * @brief Convert to CSD (Canonical Signed Digit) string representation
 *
 * @param[in] num
 * @param[in] nnz number of non-zero
 * @return string
 */
auto to_csdfixed(double num, unsigned int nnz = 4) noexcept -> string
{
    if (num == 0.)
    {
        return "0";
    }
    auto an = fabs(num);
    auto n = an < 1. ? 0 : int(ceil(log2(an * 1.5)));
    auto csd_str = string {an < 1. ? "0" : ""};
    auto limit = pow(2., n) / 3.;
    while (n > 0 || nnz > 0)
    {
        if (n == 0)
        {
            csd_str += '.';
        }
        n -= 1;
        if (num > limit)
        {
            csd_str += '+';
            num -= 1.5 * limit;
            nnz -= 1;
        }
        else if (num < -limit)
        {
            csd_str += '-';
            num += 1.5 * limit;
            nnz -= 1;
        }
        else
        {
            csd_str += '0';
        }
        limit /= 2.;
        if (nnz == 0)
        {
            num = 0;
        }
    }
    return csd_str;
}
