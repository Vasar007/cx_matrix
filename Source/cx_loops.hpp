#pragma once

#include <utility>
#include <type_traits>
#include <initializer_list>


/*
 *
 * Get a hint from https://stackoverflow.com/questions/37602057/why-isnt-a-for-loop-a-compile-time-expression
 *
 */

namespace cx_loops
{

template <typename F, std::size_t... S>
constexpr void static_for_s(F&& function, std::index_sequence<S...>)
{
    int unpack[] =
    {
        0,
        (function(std::integral_constant<std::size_t, S>{}), void(), 0)...
    };

    (void)unpack;
}

template <std::size_t iterations, typename F>
constexpr void static_for_s(F&& function)
{
    static_for_s(std::forward<F>(function), std::make_index_sequence<iterations + 1u>());
}


// Not tested this variant of implementation constexpr for-loop.
template <typename T>
void static_consume(std::initializer_list<T>)
{ }

template<typename F, std::size_t... S>
constexpr void static_for(F&& function, std::index_sequence<S...>)
{
    return static_consume({ (function(std::integral_constant<std::size_t, S>{}), 0)... });
}

} // namespace cx_loops
