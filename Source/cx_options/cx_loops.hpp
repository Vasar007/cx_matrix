// Copyright (C) 2018 Vasily Vasilyev (vasar007@yandex.ru)

#pragma once

#include <utility>
#include <type_traits>
#include <initializer_list>


/*
 * Get a hint from https://stackoverflow.com/questions/37602057/why-isnt-a-for-loop-a-compile-time-expression
 *
 * Example of usage:
 *
 *  cx_loops::static_for<10>([] (const auto Index)
 *  {
 *      if constexpr (Index > 2u)
 *      {
 *          constexpr auto A = generate_matrix_with_diagonal_predominance<double, Index, Index>();
 *          constexpr auto b = generate_matrix<double, Index, 1>();
 *          constexpr auto x1 = jacobi_solve(A, b, kEps);
 *          constexpr auto x2 = seidel_solve(A, b, kEps);
 *          std::cout << "Matrix A" << Index << ":\n" << A << "\n\n";
 *          std::cout << "Matrix b" << Index << ":\n" << b << "\n\n";
 *          std::cout << "Matrix x" << Index << ":\n" << x1.first << "\n\n";
 *          std::cout << "Number of iterations Jacobi" << Index << ":\n" << x1.second << "\n\n";
 *          std::cout << "Matrix x" << Index << ":\n" << x2.first << "\n\n";
 *          std::cout << "Number of iterations Seidel" << Index << ":\n" << x2.second;
 *          std::cout << "\n\n------------------------------------\n\n";
 *      }
 *      else
 *      {
 *          std::cout << "Skip redundant steps " << (3u - Index) << "...\n";
 *      }
 *  });
 *
 */

namespace vv::cx_loops
{

template <class F, std::size_t... S>
constexpr void static_for(F&& function, std::index_sequence<S...>)
{
    int unpack[] =
    {
        0,
        (function(std::integral_constant<std::size_t, S>{}), void(), 0)...
    };

    (void)unpack;
}

template <std::size_t iterations, class F>
constexpr void static_for(F&& function)
{
    static_for(std::forward<F>(function), std::make_index_sequence<iterations + 1u>());
}

} // namespace vv::cx_loops
