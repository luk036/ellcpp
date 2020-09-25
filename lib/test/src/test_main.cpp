#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

// #include <boost/coroutine2/all.hpp>
// #include <iostream>
// // #include <boost/coroutine2/detail/push_coroutine.hpp>
// // #include <boost/coroutine2/detail/pull_coroutine.hpp>
// using namespace std;

// // Method 1

// void foo(boost::coroutines2::coroutine<int>::pull_type& sink)
// {
//     // using coIter =
//     boost::coroutines2::coroutine<int>::pull_type::iterator; for (auto start
//     = begin(sink); start != end(sink); ++start)
//     {
//         std::cout << "retrieve " << *start << "\n";
//     }
// }

// // Method 2

// void foo2(boost::coroutines2::coroutine<int>::pull_type& sink)
// {
//     for (auto val : sink)
//     {
//         std::cout << "retrieve " << val << "\n";
//     }
// }


// TEST_CASE("Test Coroutine")
// {
//     using coro_t = boost::coroutines2::coroutine<int>;
//     int max = 8;

//     coro_t::pull_type source([&](coro_t::push_type& sink) {
//         int first = 1;
//         int second = 1;
//         sink(first);
//         sink(second);
//         for (int i = 0; i < max; ++i)
//         {
//             int third = first + second;
//             first = second;
//             second = third;
//             sink(third);
//         }
//     });

//     for (auto i : source)
//     {
//         cout << i << " ";
//     }
//     cout << endl;

//     coro_t::push_type sink([&](coro_t::pull_type& source) {
//         while (source)
//         {
//             cout << source.get() << " ";
//             source();
//         }
//     });
//     vector<int> v {1, 1, 2, 3, 5, 8, 13, 21, 34, 55};
//     copy(begin(v), end(v), begin(sink));

//     foo(source);
// }