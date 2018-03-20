#pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>

#include "matrix.hpp"


namespace vv
{

/// ===== FUNCTION SECTION =====
template <class Type, std::size_t Rows, std::size_t Columns>
[[deprecated("Duplicated function with different steps. Use 'lu_decompose'")]]
constexpr std::array<matrix<Type, Rows, Columns>, 2u>
    lu_decomposition(matrix<Type, Rows, Columns> mat)
{
    static_assert(Rows == Columns, "Matrix is not quadrant!");

    constexpr std::size_t N = Rows;
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

    matrix<Type, Rows, Columns> lower{};
    matrix<Type, Rows, Columns> upper{};

    for (size_type i = 0u; i < N; ++i)
    {
        // Diagonal as 1.
        lower(i, i) = static_cast<Type>(1);
    }

    // Decomposing matrix into Upper and Lower triangular matrix.
    for (size_type k = 0u; k < N; ++k)
    {
        upper(k, k) = mat(k, k);
        for (size_type i = k; i < N; ++i)
        {
            // mat(i, k) == Vi
            lower(i, k) = mat(i, k) / mat(k, k);
            // mat(k, i) == Wi
            upper(k, i) = mat(k, i);             
        }
 
        for (size_type i = k; i < N; ++i)
        {
            for (size_type j = k; j < N; ++j)
            {
                mat(i, j) -= lower(i, k) * upper(k, j);
            }
        }
    } // for

    return { lower, upper };
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr std::array<matrix<Type, Rows, Columns>, 2u>
    lu_decompose(const matrix<Type, Rows, Columns>& mat)
{
    static_assert(Rows == Columns, "Matrix is not quadrant!");

    constexpr std::size_t N = Rows;
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

    matrix<Type, Rows, Columns> lower{};
    matrix<Type, Rows, Columns> upper{};

    // Decomposing matrix into Upper and Lower triangular matrix.
    for (size_type i = 0u; i < N; ++i)
    {
        // Upper Triangular.
        for (size_type k = i; k < N; ++k)
        {
            // Summation of L(i, j) * U(j, k).
            Type sum{};
            for (size_type j = 0u; j < i; ++j)
            {
                sum += (lower(i, j) * upper(j, k));
            }
 
            // Evaluating U(i, k).
            upper(i, k) = mat(i, k) - sum;
        } // for
 
        // Lower Triangular.
        for (size_type k = i; k < N; ++k)
        {
            if (i == k)
            {
                // Diagonal as 1.
                lower(i, i) = static_cast<Type>(1);
            }
            else
            {
                // Summation of L(k, j) * U(j, i).
                Type sum{};
                for (size_type j = 0u; j < i; ++j)
                {
                    sum += (lower(k, j) * upper(j, i));
                }
                // Evaluating L(k, i).
                lower(k, i) = (mat(k, i) - sum) / upper(i, i);
            }
        } // for
    } // for (size_type i = 0u; i < N; ++i)

    return { lower, upper };
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr long rank(const matrix<Type, Rows, Columns>& mat)
{
    static_assert(Rows == Columns, "Matrix is not quadrant!");

    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;
    constexpr std::size_t N = Rows;
    
    long rank = 0;
    const auto [L, U] = lu_decompose(mat);
    
    for (size_type i = 0u; i < N; ++i)
    {
        if (U(i, N - 1u) != static_cast<Type>(0))
        {
            ++rank;
        }
    }
    return rank;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr Type lu_determenant(const matrix<Type, Rows, Columns>& mat)
{
    static_assert(Rows == Columns, "Matrix is not quadrant!");

    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;
    constexpr std::size_t N = Rows;

    const auto [L, U] = lu_decompose(mat);

    Type detL = L(0u, 0u);
    Type detU = U(0u, 0u);
    for (size_type i = 1u; i < N; ++i)
    {
        detL *= L(i, i);
        detU *= U(i, i);
    }

    return detL * detU;
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_LU

TEST_METHOD(lu_test_major_ver)
{
    constexpr auto lu_value = lu_decompose(mat_A);
    std::cout << "LU-decomposition:\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix L:\n" << lu_value.at(0u) << "\n\n";
    std::cout << "Matrix U:\n" << lu_value.at(1u);
    std::cout << "\n\n------------------------------------\n\n";

    constexpr auto lu_value_check = (lu_value.at(0u) * lu_value.at(1u)) - mat_A;
    std::cout << "Checking LU-decomposition:\n\n";
    std::cout << "Matrix (L * U) - A:\n" << lu_value_check;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lu_test_minor_ver)
{
    constexpr auto lu_value = lu_decomposition(mat_A);
    std::cout << "LU-decompose:\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix L:\n" << lu_value.at(0u) << "\n\n";
    std::cout << "Matrix U:\n" << lu_value.at(1u);
    std::cout << "\n\n------------------------------------\n\n";

    constexpr auto lu_value_check = (lu_value.at(0u) * lu_value.at(1u)) - mat_A;
    std::cout << "Checking LU-decompose:\n\n";
    std::cout << "Matrix (L * U) - A:\n" << lu_value_check;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lu_test_determenant)
{
    constexpr auto lu_det = lu_determenant(mat_A);
    std::cout << "Determenant by LU:\n\n";
    std::cout << "Determenant of A:\n" << lu_det;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lu_test_rank)
{
    constexpr auto lu_value = vv::lu_decompose(mat_A);
    constexpr auto rank = vv::rank(mat_A);
    std::cout << "LU-decomposition:\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix L:\n" << lu_value.at(0u) << "\n\n";
    std::cout << "Matrix U:\n" << lu_value.at(1u) << "\n\n";
    std::cout << "Rank A:\n" << rank;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_LU

} // namespace vv
