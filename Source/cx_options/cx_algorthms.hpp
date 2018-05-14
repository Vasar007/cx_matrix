#pragma once

#include <algorithm>
#include <type_traits>


namespace vv
{

namespace detail
{

/// ===== ADDITION FUNCTIONAL SECTION =====
template <class InputIt, class T>
constexpr T accumulate(InputIt first, InputIt last, T init)
{
    for (; first != last; ++first)
    {
        init = init + *first;
    }
    return init;
}


template <class InputIt, class T, class BinaryOperation>
constexpr T accumulate(InputIt first, InputIt last, T init, BinaryOperation op)
{
    for (; first != last; ++first)
    {
        init = op(init, *first);
    }
    return init;
}


template <class InputIt1, class InputIt2>
constexpr bool equal(InputIt1 first1, InputIt1 last1, InputIt2 first2)
{
    for (; first1 != last1; ++first1, ++first2)
    {
        if (!(*first1 == *first2))
        {
            return false;
        }
    }
    return true;
}


template <class InputIt1, class InputIt2, class BinaryPredicate>
constexpr bool equal(InputIt1 first1, InputIt1 last1, InputIt2 first2, BinaryPredicate p)
{
    for (; first1 != last1; ++first1, ++first2)
    {
        if (!p(*first1, *first2))
        {
            return false;
        }
    }
    return true;
}


template <class ForwardIt, class T>
constexpr void fill(ForwardIt first, ForwardIt last, const T& value)
{
    for (; first != last; ++first)
    {
        *first = value;
    }
}


template <class T, class U = T>
constexpr T exchange(T& obj, U&& new_value)
{
    T old_value = std::move(obj);
    obj = std::forward<U>(new_value);
    return old_value;
}


template <class T>
constexpr void swap(T& a, T& b) noexcept(std::is_nothrow_move_constructible_v<T>
                                         && std::is_nothrow_move_assignable_v<T>)
{
    T temp = std::move(a);
    a = std::move(b);
    b = std::move(temp);
}


template <class InputIt, class T>
constexpr InputIt find(InputIt first, InputIt last, const T& value)
{
    for (; first != last; ++first)
    {
        if (*first == value)
        {
            return first;
        }
    }
    return last;
}


template <class InputIt, class UnaryPredicate>
constexpr InputIt find_if(InputIt first, InputIt last, UnaryPredicate p)
{
    for (; first != last; ++first)
    {
        if (p(*first))
        {
            return first;
        }
    }
    return last;
}


template <class InputIt, class UnaryPredicate>
constexpr InputIt find_if_not(InputIt first, InputIt last, UnaryPredicate q)
{
    for (; first != last; ++first)
    {
        if (!q(*first))
        {
            return first;
        }
    }
    return last;
}


template <class InputIt, class UnaryPredicate>
constexpr bool all_of(InputIt first, InputIt last, UnaryPredicate p)
{
    return std::find_if_not(first, last, p) == last;
}


template <class InputIt, class UnaryPredicate>
constexpr bool any_of(InputIt first, InputIt last, UnaryPredicate p)
{
    return std::find_if(first, last, p) != last;
}


template <class InputIt, class UnaryPredicate>
constexpr bool none_of(InputIt first, InputIt last, UnaryPredicate p)
{
    return std::find_if(first, last, p) == last;
}

} // namespace detail

} // namespace vv
