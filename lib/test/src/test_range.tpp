#include "generator.hpp"
#include <stdio.h>


template <typename T>
generator<int> range(T first, T last)
{
    while (first != last)
    {
        co_yield first++;
    }
}

template <typename T>
generator<int> range1(T first, T last)
{
    while (first != last)
    {
        throw std::logic_error("BEGIN");

        co_yield first++;
    }
}

template <typename T>
generator<int> range2(T first, T last)
{
    while (first != last)
    {
        co_yield first++;

        throw std::logic_error("ITERATOR");
    }
}

template <typename T>
generator<int> range4(T first, T last)
{
    co_return;
}

int main()
{
    printf("\nrange\n");

    try
    {
        for (int i : range(0, 10))
        {
            printf("%d\n", i);
        }
    }
    catch (std::exception const& e)
    {
        printf("%s\n", e.what());
    }

    printf("\nrange1\n");

    try
    {
        for (int i : range1(0, 10))
        {
            printf("%d\n", i);
        }
    }
    catch (std::exception const& e)
    {
        printf("%s\n", e.what());
    }

    printf("\nrange2\n");

    try
    {
        for (int i : range2(0, 10))
        {
            printf("%d\n", i);
        }
    }
    catch (std::exception const& e)
    {
        printf("%s\n", e.what());
    }

    try
    {
        for (int i : range4(0, 10))
        {
            printf("%d\n", i);
        }
    }
    catch (std::exception const& e)
    {
        printf("%s\n", e.what());
    }
}
