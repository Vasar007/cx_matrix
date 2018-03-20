#pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>

#include "matrix.hpp"
#include "cx_math.h"


namespace vv
{

/// ===== FUNCTION SECTION =====
template <class Type, std::size_t Rows_A, std::size_t Columns_A>
[[deprecated("Duplicated function with different permutation matrix. Use 'lup_decompose'")]]
constexpr std::tuple<matrix<Type, Rows_A, Columns_A>, matrix<Type, Rows_A, Columns_A>, long>
    lup_decomposition(const matrix<Type, Rows_A, Columns_A>& A)
{
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows_A, Columns_A>::size_type;

    auto C = A;
    // Create identity matrix.
    auto P = matrix<Type, Rows_A, Columns_A>::create_identity();

    long permutations_counter = 0L;
    for (size_type i = 0; i < Columns_A; ++i)
    {
        // Find pivot element.
        Type pivotValue{};
        size_type pivot{};
        for (size_type row = i; row < Rows_A; ++row)
        {
            const auto absValue = static_cast<Type>(cx::abs(C(row, i)));
            if (absValue > pivotValue)
            {
                pivotValue = absValue;
                pivot = row;
            }
        }

        if (pivotValue != static_cast<Type>(0))
        {
            // Swap i row with pivot element.
            P.swap_rows(pivot, i);
            C.swap_rows(pivot, i);
            ++permutations_counter;

            for (size_type j = i + 1; j < Rows_A; ++j)
            {
                C(j, i) /= C(i, i);
                for (size_type k = i + 1; k < Columns_A; ++k)
                {
                    C(j, k) -= C(j, i) * C(i, k);
                }
            }
        }
    } // for (size_type i = 0; i < Columns_A; ++i)
    
    // C = L + U - E
    return { C, P, permutations_counter };
}


template <class Type, std::size_t Rows_P, std::size_t Columns_P>
[[deprecated("Function can contain bugs! Don't use it.")]]
constexpr auto convert_to_array(const matrix<Type, Rows_P, Columns_P>& P)
{
    static_assert(Rows_P == Columns_P, "Matrix is not quadrant!");
    
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows_P, Columns_P>::size_type;
    constexpr std::size_t N = Rows_P;

    matrix<Type, Rows_P, 1u> result{};
    for (size_type i = 0u; i < N; ++i)
    {
        for (size_type j = 0u; j < N; ++j)
        {
            if (P(i, j) == static_cast<Type>(1))
            {
                result(i, 0u) = static_cast<Type>(j);
                break;
            }
            if (j == N - 1u)
            {
                throw std::domain_error("Invalid matrix P!"); 
            }
        }
    }

    return result;
}


template <class Type, std::size_t Rows, std::size_t Columns = 1u>
constexpr auto convert_to_matrix(const matrix<Type, Rows, 1u>& P)
{
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

	matrix<Type, Rows, Rows> result_P{};
	for (size_type i = 0u; i < Rows; ++i)
	{
		result_P(i, P(i, 0u)) = static_cast<Type>(1);
	}
	return result_P;
}


template <class Type, std::size_t Rows_A, std::size_t Columns_A>
constexpr std::tuple<matrix<Type, Rows_A, Columns_A>, matrix<Type, Rows_A, 1u>, long>
    lup_decompose(matrix<Type, Rows_A, Columns_A> A, const Type eps = kDefault_eps<Type>)
{
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows_A, Columns_A>::size_type;

    long permutations_counter = static_cast<long>(Rows_A);
    matrix<Type, Rows_A, 1u> P{};
    // Unit permutation matrix, P[N] initialized with N.
    for (size_type i = 0u; i < Rows_A; ++i)
    {
        P(i, 0u) = static_cast<Type>(i);
    }
    
    for (size_type i = 0u; i < Columns_A; ++i)
    {
        Type max_A{};
        size_type i_max = i;
        for (size_type k = i; k < Rows_A; ++k)
        {
            const auto absA = cx::abs(A(k, i));
            if (absA > max_A)
            { 
                max_A = absA;
                i_max = k;
            }
        }

        if (max_A < cx::abs(eps))
        { 
            // Failure, matrix is degenerate.
            //throw std::domain_error("Matrix is singular!");
            return { matrix<Type, Rows_A, Columns_A>{}, matrix<Type, Rows_A, 1u>{}, -1L };
        }

        if (i_max != i)
        {
            // Pivoting P.
            P.swap_rows(i, i_max);
            // Pivoting rows of A.
            A.swap_rows(i, i_max);
            // Counting pivots starting from N (for determinant).
            ++permutations_counter;
        }

        for (size_type j = i + 1u; j < Rows_A; ++j)
        {
            A(j, i) /= A(i, i);
            for (size_type k = i + 1u; k < Columns_A; ++k)
            {
                A(j, k) -= A(j, i) * A(i, k);
            }
        }
    } // for (size_type i = 0u; i < N; ++i)
    
    return { A, P, permutations_counter };
}


template <class Type, std::size_t Rows, std::size_t Columns_A, std::size_t Columns_b>
constexpr auto lup_solve(const matrix<Type, Rows, Columns_A>& mat,
                         const matrix<Type, Rows, Columns_b>& b,
                         const Type eps = kDefault_eps<Type>)
{
    static_assert(Columns_b == 1u, "Matrix contains more than one columns in right "
                                   "hand side vector!");
                                   
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns_A>::size_type;

    const auto [A, P, permutations_counter] = lup_decompose(mat, eps);

    matrix<Type, Columns_A, Columns_b> x{}; // Columns_b == 1u
    for (size_type i = 0u; i < Rows && i < Columns_A; ++i)
    {
        x(i, 0u) = b(static_cast<size_type>(P(i, 0u)), 0u);

        for (size_type k = 0u; k < i; ++k)
        {
            x(i, 0u) -= A(i, k) * x(k, 0u);
        }
    }

    for (size_type i = Columns_A - 1u; i < Rows && i < Columns_A; --i)
    {
        for (size_type k = i + 1u; k < Rows && k < Columns_A; ++k)
        {
            x(i, 0u) -= A(i, k) * x(k, 0u);
        }

        x(i, 0u) /= A(i, i);

        if (i == 0u) break;
    }
    
    return x;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr Type lup_determenant(const matrix<Type, Rows, Columns>& mat,
                           const Type eps = kDefault_eps<Type>)
{
    static_assert(Rows == Columns, "Matrix is not quadrant!");

    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;
    constexpr std::size_t N = Rows;

    const auto [A, P, permutations_counter] = lup_decompose(mat, eps); 
    Type det = A(0u, 0u);
    for (size_type i = 1u; i < N; ++i)
    {
        det *= A(i, i);
    }

    return ((permutations_counter - N) % 2 == 0) ? det : -det;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr auto lup_invert(const matrix<Type, Rows, Columns>& mat,
                          const Type eps = kDefault_eps<Type>)
{
    static_assert(Rows == Columns, "Matrix is not quadrant!");

    constexpr std::size_t N = Rows;
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

    const auto [A, P, permutations_counter] = lup_decompose(mat, eps);

    matrix<Type, Rows, Columns> IA{};
    for (size_type j = 0u; j < N; ++j)
    {
        for (size_type i = 0u; i < N; ++i)
        {
            if (P(i, 0u) == j)
            { 
                IA(i, j) = static_cast<Type>(1);
            }
            else
            {
                IA(i, j) = static_cast<Type>(0);
            }

            for (size_type k = 0u; k < i; ++k)
            {
                IA(i, j) -= A(i, k) * IA(k, j);
            }
        } // for (size_type i = 0u; i < N; ++i)

        for (size_type i = N - 1u; ; --i)
        {
            for (size_type k = i + 1u; k < N; ++k)
            {
                IA(i, j) -= A(i, k) * IA(k, j);
            }

            IA(i, j) = IA(i, j) / A(i, i);

            if (i == 0u) break;
        } // for (size_type i = N - 1u; ; --i)
    } //  for (size_type j = 0u; j < N; ++j)

    return IA;
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_LUP

TEST_METHOD(lup_test_major_ver)
{
    constexpr auto lup_tuple = lup_decompose(mat_A);
    constexpr auto mat_C = std::get<0u>(lup_tuple);
    constexpr auto mat_P = std::get<1u>(lup_tuple);
    constexpr auto pivots_number = std::get<2u>(lup_tuple);
    std::cout << "LUP-decompose:\n\n";
    std::cout << "Changed matrix A:\n" << mat_C << "\n\n";
    std::cout << "Matrix P:\n" << mat_P << "\n\n";
    std::cout << "Number of pivots:\n" << pivots_number;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lup_test_solve)
{
    constexpr auto vec_x = lup_solve(mat_A, vec_b);
    std::cout << "LUP-solve:\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix b:\n" << vec_b << "\n\n";
    std::cout << "Matrix X:\n" << vec_x;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking LUP-solve:\n\n";
    constexpr auto lup_solve_check = (mat_A * vec_x) - vec_b;
    std::cout << "Result of (A * x - b):\n" << lup_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}


#endif // ENABLE_TESTS_LUP

#if ENABLE_TESTS_LUP_ONLY_SQUARE

TEST_METHOD(lup_test_minor_ver)
{
    constexpr auto lup_tuple = lup_decomposition(mat_A);
    constexpr auto mat_C = std::get<0u>(lup_tuple);
    constexpr auto mat_P = std::get<1u>(lup_tuple);
    constexpr auto pivots_number = std::get<2u>(lup_tuple);
    std::cout << "LUP-decomposition:\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix (L + U - E):\n" <<  mat_C << "\n\n";
    std::cout << "Matrix P:\n" <<  mat_P << "\n\n";
    std::cout << "Number of pivots:\n" << pivots_number;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lup_test_determenant)
{
    constexpr auto lup_det = lup_determenant(mat_A);
    std::cout << "Determenant by LUP:\n\n";
    std::cout << "Determenant of A:\n" << lup_det;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lup_test_inverse)
{
    constexpr auto lup_invert_mat = lup_invert(mat_A);
    std::cout << "Invert by LUP:\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Inverted matrix to A:\n" << lup_invert_mat;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking inverse:\n\n";
    constexpr auto inverse_check1 = mat_A * lup_invert_mat;
    std::cout << "Result of (A * IA):\n" << inverse_check1 << "\n\n";

    constexpr auto inverse_check2 = lup_invert_mat * mat_A;    
    std::cout << "Result of (IA * A):\n" << inverse_check2;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_LUP_ONLY_SQUARE

} // namespace vv
