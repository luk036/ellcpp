#pragma once

#if defined(__clang__) || defined(__GNUC__)
#define ELL_LIKELY(x) __builtin_expect(!!(x), 1)
#define ELL_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define ELL_LIKELY(x) (!!(x))
#define ELL_UNLIKELY(x) (!!(x))
#endif
