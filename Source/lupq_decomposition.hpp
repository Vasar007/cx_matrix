#pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>

#include "matrix.hpp"
#include "cx_math.h"


namespace vv::detail
{

/// ===== ADDITION FUNCTIONAL SECTION =====
template <class Type, std::size_t Rows, std::size_t Columns>
constexpr matrix<Type, Rows, Columns> lupq_get_u_solution_impl(const matrix<Type, Rows, Columns>& C,
															   const std::size_t rank) noexcept
{
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;
	
	matrix<Type, Rows, Columns> result_U{};
	// Singular matrix.
	if (rank < Rows || rank < Columns)
	{
		for (size_type i = 0u; i < rank; ++i)
		{
			for (size_type j = i; j < rank; j++)
			{
				result_U(i, j) = C(i, j);
			}
			for (size_type j = rank; j < Columns; j++)
			{
				result_U(i, j) = static_cast<Type>(0);
			}
		}
		for (size_type i = rank; i < Rows; i++)
		{
			for (size_type j = i; j < Columns; j++)
			{
				result_U(i, j) = static_cast<Type>(0);
			}
		}
	}
	// Regular matrix.
	else
	{
		for (size_type i = 0u; i < Rows; ++i)
		{
			for (size_type j = i; j < Columns; ++j)
			{
				result_U(i, j) = C(i, j);
			}
		}
	}

	return result_U;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr matrix<Type, Rows, Rows> lupq_get_l_solution_impl(const matrix<Type, Rows, Columns>& C,
															const std::size_t rank) noexcept
{
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;
	
	matrix<Type, Rows, Rows> result_L{};
	// Singular matrix.
	if (rank < Rows)
	{
		for (size_type i = 0u; i < rank; ++i)
		{
			result_L(i, i) = static_cast<Type>(1);

			if (i == 0u) continue;

			for (size_type j = i - 1u; ; --j)
			{
				result_L(i, j) = C(i, j);

				if (j == 0u) break;
			}
		}
		for (size_type i = rank; i < Rows; ++i)
		{
			if (i == 0u) continue;

			for (size_type j = i; ; --j)
			{
				result_L(i, j) = static_cast<Type>(0);

				if (j == 0u) break;
			}
		}
	}
	// Regular matrix.
	else
	{
		for (size_type i = 0; i < Rows; ++i)
		{
			result_L(i, i) = static_cast<Type>(1);

			if (i == 0u) continue;

			for (size_type j = i - 1u; ; --j)
			{
				result_L(i, j) = C(i, j);

				if (j == 0u) break;
			}
		}
	}

	return result_L;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr matrix<Type, Rows, Columns> lupq_get_u_impl(const matrix<Type, Rows, Columns>& C) noexcept
{
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

	matrix<Type, Rows, Columns> result_U{};
	for (size_type i = 0; i < Rows; ++i)
	{
		for (size_type j = i; j < Columns; ++j)
		{
			result_U(i, j) = C(i, j);
		}
	}

	return result_U;
}

template <class Type, std::size_t Rows, std::size_t Columns>
constexpr matrix<Type, Rows, Rows> lupq_get_l_impl(const matrix<Type, Rows, Columns>& C) noexcept
{
	// Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

	matrix<Type, Rows, Rows> result_L{};
	for (size_type i = 0u; i < Rows; ++i)
	{
		result_L(i, i) = static_cast<Type>(1);

		if (i == 0u) continue;

		size_type temp_index = i >= Columns ? Columns - 1u : i;
		for (size_type j = temp_index - 1u; ; --j)
		{
			if (Columns < Rows)
			{
				result_L(i, j) = static_cast<Type>(0);
			}
			else
			{
				result_L(i, j) = C(i, j);
			}

			if (j == 0u) break;
		}
	}
	
	return result_L;
}


template <class Type, std::size_t Rows, std::size_t Columns = 1u>
constexpr matrix<Type, Rows, Rows> lupq_get_p_impl(const matrix<Type, Rows, 1u>& P) noexcept
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


template <class Type, std::size_t Columns, std::size_t Rows = 1u>
constexpr matrix<Type, Columns, Columns> lupq_get_q_impl(
												const matrix<Type, Columns, 1u>& Q) noexcept
{
	// Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

	matrix<Type, Columns, Columns> result_Q{};
	for (size_type i = 0u; i < Columns; ++i)
	{
		result_Q(Q(i, 0u), i) = static_cast<Type>(1);
	}
	return result_Q;
}

} // namespace vv::detail


namespace vv
{

/// ===== FUNCTION SECTION =====
template <class Type, std::size_t Rows, std::size_t Columns>
constexpr std::tuple<matrix<Type, Rows, Columns>, matrix<Type, Rows, 1u>, matrix<Type, Columns, 1u>,
					 std::size_t, long>
    lupq_decompose(matrix<Type, Rows, Columns> mat, const Type eps = kDefault_eps<Type>)
{
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

	long permutations_counter = 0L;

	matrix<Type, Rows, 1u> P{};
	matrix<Type, Columns, 1u> Q{};

	for (size_type i = 0; i < Rows; ++i)
    {
        P(i, 0u) = static_cast<Type>(i);
    }

	for (size_type i = 0; i < Columns; ++i)
    {
        Q(i, 0u) = static_cast<Type>(i);
    }

	size_type row_index{};
	size_type column_index{};

	// Calculate matrix C = L + U - E, matrix P, matrix Q, rank and permutations_counter.
	while (!mat.is_empty(row_index, column_index) && row_index < Rows && column_index < Columns)
	{
		const auto indices_pair = mat.find_max_element(row_index, column_index);
		if (row_index != indices_pair.first)
		{
			mat.swap_rows(row_index, indices_pair.first);
			P.swap_rows(row_index, indices_pair.first);
			++permutations_counter;
		}
		
		if (column_index != indices_pair.second)
		{
			mat.swap_columns(column_index, indices_pair.second);
			Q.swap_rows(column_index, indices_pair.second);
			++permutations_counter;
		}

		for (size_type j = row_index + 1u; j < Rows; ++j)
		{
			//if (cx::abs(mat(j, column_index)) < eps) continue;

			mat(j, column_index) /= mat(row_index, column_index);
			for (size_type k = column_index + 1u; k < Columns; ++k)					
			{
				mat(j, k) -= mat(row_index, k) * mat(j, column_index);
			}
		}
		++row_index;
		++column_index;
	}
	size_type rank = std::max(row_index, column_index);

    return { mat, P, Q, rank, permutations_counter };
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr std::tuple<matrix<Type, Rows, Columns>, matrix<Type, Rows, Rows>,
					 matrix<Type, Rows, Rows>, matrix<Type, Columns, Columns>, long>
	lupq_get_components(const matrix<Type, Rows, Columns>& mat,
			  		    const Type eps = kDefault_eps<Type>)
{
    const auto [C, P, Q, rank, permutations_counter] = lupq_decompose(mat, eps);

	matrix<Type, Rows, Columns> result_U = detail::lupq_get_u_impl(C);
	matrix<Type, Rows, Rows> result_L = detail::lupq_get_l_impl(C);
	matrix<Type, Rows, Rows> result_P = detail::lupq_get_p_impl(P);
	matrix<Type, Columns, Columns> result_Q = detail::lupq_get_q_impl(Q);

	return { result_U, result_L, result_P, result_Q, rank };
}


template <class Type, std::size_t Rows, std::size_t Columns_A, std::size_t Columns_b>
constexpr auto lupq_solve(const matrix<Type, Rows, Columns_A>& A,
                          const matrix<Type, Rows, Columns_b>& b,
                          const Type eps = kDefault_eps<Type>)
{
	static_assert(Columns_b == 1u, "Matrix contains more than one columns in right "
                                   "hand side vector!");

	// Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns_A>::size_type;

	const auto [C, _P, _Q, rank, permutations_counter] = lupq_decompose(A, eps);
	matrix<Type, Rows, Rows> P = detail::lupq_get_p_impl(_P);
	matrix<Type, Columns_A, Columns_A> Q = detail::lupq_get_q_impl(_Q);

	matrix<Type, Columns_A, Columns_b> solution{}; // Columns_b == 1u
	solution = P * b;

	// Solve equation Lz = b
	matrix<Type, Rows, Rows> L = detail::lupq_get_l_solution_impl(C, rank);
	for (size_type i = 0u; i < rank - 1u; ++i)
	{
		for (size_type j = i + 1u; j < rank; ++j)
		{
			solution(j, 0u) -= solution(i, 0u) * L(j, i);
		}
	}
	
	// Solve system Uy = b, xQ = y
	matrix<Type, Rows, Columns_A> U = detail::lupq_get_u_solution_impl(C, rank);
	for (size_type i = rank - 1u; i > 0u; --i)
	{
		solution(i, 0u) /= U(i, i);

		for (size_type j = i - 1u; ; --j)
		{
			solution(j, 0u) -= solution(i, 0u) * U(j, i);

			if (j == 0u) break;
		}
	}
	solution(0u, 0u) /= U(0u, 0u); // To make U(0, 0) == Type{1}

	// Checking the compatibility of the system of equations.
	if (rank < Rows)
	{
		matrix<Type, Rows, Columns_A> origin = P * A * Q;
		for (size_type i = rank; i < Rows && i < Columns_A; ++i)
		{
			Type x = solution(i, 0);
			for (size_type j = 0u; j < rank && j < Columns_A; ++j)
			{
				x -= origin(i, j) * solution(j, 0u);
			}

			if (cx::abs(x) > eps)
			{
				return matrix<Type, Rows, Columns_b>::get_error_matrix();
			}
		}
	}

	for (size_type i = rank; i < Columns_A; ++i)
	{
		solution(i, 0) = static_cast<Type>(0);
	}

	return Q * solution;
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_LUPQ

TEST_METHOD(lupq_test_decompose)
{
    constexpr auto lupq_tuple = lupq_get_components(mat_A);
    constexpr auto mat_U = std::get<0u>(lupq_tuple);
    constexpr auto mat_L = std::get<1u>(lupq_tuple);
    constexpr auto mat_P = std::get<2u>(lupq_tuple);
    constexpr auto mat_Q = std::get<3u>(lupq_tuple);
    constexpr auto rank = std::get<4u>(lupq_tuple);
    std::cout << "LUPQ-decompose:\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix U:\n" << mat_U << "\n\n";
    std::cout << "Matrix L:\n" << mat_L << "\n\n";
    std::cout << "Matrix P:\n" << mat_P << "\n\n";
    std::cout << "Matrix Q:\n" << mat_Q << "\n\n";
    std::cout << "Rank A:\n" << rank;
    std::cout << "\n\n------------------------------------\n\n";

	constexpr auto lupq_value_check = (mat_L * mat_U) - (mat_P * mat_A * mat_Q);
    std::cout << "Checking LUPQ-decompose:\n\n";
    std::cout << "Matrix (L * U) - (P * A * Q):\n" << lupq_value_check;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(lupq_test_solve)
{
    constexpr auto vec_x = lupq_solve(mat_A, vec_b);
    std::cout << "LUPQ-solve:\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix b:\n" << vec_b << "\n\n";
    std::cout << "Matrix X:\n" << vec_x;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking LUPQ-solve:\n\n";
    constexpr auto lupq_solve_check = (mat_A * vec_x) - vec_b;
    std::cout << "Result of (A * x - b):\n" << lupq_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_LUPQ

} //namesace vv
