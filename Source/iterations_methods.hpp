#pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>
#include <thread>
#include <chrono>

#include "utils.hpp"
#include "cx_loops.hpp"
#include "cx_random.hpp"
#include "matrix.hpp"
#include "cx_math.h"


namespace vv
{

constexpr long kMax_iterations = 10'000L;


/// ===== FUNCTION SECTION =====
template <class Type, std::size_t Rows, std::size_t Columns_A, std::size_t Columns_b>
constexpr std::pair<matrix<Type, Rows, Columns_b>, long>
    jacobi_solve(const matrix<Type, Rows, Columns_A>& A,
                 const matrix<Type, Rows, Columns_b>& b,
                 const Type eps = kDefault_eps<Type>) noexcept
{
    static_assert(Rows == Columns_A, "Matrix is not quadrant!");
    static_assert(Columns_b == 1u, "Matrix contains more than one columns in right "
                                   "hand side vector!");

    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns_A>::size_type;
    constexpr std::size_t N = Rows;

    long iterations_counter = 0L;

	Type norm{};
    matrix<Type, Rows, Columns_b> x{}; // Columns_b == 1u
    matrix<Type, Rows, Columns_b> previous_x{}; // Columns_b == 1u
	do
    {
		for (size_type i = 0; i < N; ++i)
        {
			previous_x(i, 0u) = b(i, 0u);
			for (size_type j = 0; j < N; ++j)
            {
				if (i != j)
                {
                    if (cx::abs(x(j, 0u)) > eps)
                    {
					    previous_x(i, 0u) -= A(i, j) * x(j, 0u);
                    }
                }
			}
			previous_x(i, 0u) /= A(i, i);
		}

        norm = cx::abs(x(0u, 0u) - previous_x(0u, 0u));
		for (size_type k = 0u; k < N; ++k)
        {
			if (cx::abs(x(k, 0u) - previous_x(k, 0u)) > norm)
            {
				norm = cx::abs(x(k, 0u) - previous_x(k, 0u));
            }
			x(k, 0u) = previous_x(k, 0u);
		}

        ++iterations_counter;

        if (iterations_counter > kMax_iterations)
        {
            break;
        }
	}
    while (norm > eps);

    return { x, iterations_counter };
}


template <class Type, std::size_t Rows, std::size_t Columns_A, std::size_t Columns_b>
constexpr std::pair<matrix<Type, Rows, Columns_b>, long>
    seidel_solve(const matrix<Type, Rows, Columns_A>& A,
                 const matrix<Type, Rows, Columns_b>& b,
                 const Type eps = kDefault_eps<Type>) noexcept
{
    static_assert(Rows == Columns_A, "Matrix is not quadrant!");
    static_assert(Columns_b == 1u, "Matrix contains more than one columns in right "
                                   "hand side vector!");

    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns_A>::size_type;
    constexpr std::size_t N = Rows;

    long iterations_counter = 0L;

    constexpr auto converge = [](const auto& x_k, const auto& x_k_prev, const Type epsilon)
    {
        Type norm{};
        for (size_type i = 0u; i < N; ++i)
        {
            norm += (x_k(i, 0u) - x_k_prev(i, 0u)) * (x_k(i, 0u) - x_k_prev(i, 0u));
        }

        return (cx::sqrt(norm) < epsilon);
    };

    matrix<Type, Rows, Columns_b> x{}; // Columns_b == 1u
    matrix<Type, Rows, Columns_b> previous_x{}; // Columns_b == 1u
	do
    {
        for (size_type i = 0u; i < N; ++i)
        {
            previous_x(i, 0u) = x(i, 0u);
        }

        for (size_type i = 0u; i < N; ++i)
        {
            Type var{};
            for (size_type j = 0u; j < i; ++j)
            {
                var += (A(i, j) * x(j, 0u));
            }

            for (size_type j = i + 1u; j < N; ++j)
            {
                var += (A(i, j) * previous_x(j, 0u));
            }

            x(i, 0u) = (b(i, 0u) - var) / A(i, i);
        }

        ++iterations_counter;

        if (iterations_counter > kMax_iterations)
        {
            break;
        }
    }
    while (!converge(x, previous_x, eps));

    return { x, iterations_counter };
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr matrix<Type, Rows, Columns> generate_matrix() noexcept
{
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

    matrix<Type, Rows, Columns> result{};
    for (size_type i = 0u; i < Rows; ++i)
    {
        for (size_type j = 0u; j < Columns; ++j)
        {
            const Type number = static_cast<Type>(cx_random::get_random(i * Rows + j + 1u));

            result(i, j) += number - static_cast<Type>(i * Rows + j + 1u);
            if (i + 1u < Rows)
            {
                result(i + 1u, j) = number * static_cast<Type>(0.001);
            }
        }
    }

    return result;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr matrix<Type, Rows, Columns> generate_matrix_with_diagonal_predominance() noexcept
{
    static_assert(Rows == Columns, "Matrix is not quadrant!");

    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

    matrix<Type, Rows, Columns> result{};
    for (size_type i = 0u; i < Rows; ++i)
    {
        Type sum{};
        for (size_type j = 0u; j < Columns; ++j)
        {
            const Type number = static_cast<Type>(cx_random::get_random(i * Rows + j));

            result(i, j) += number - static_cast<Type>(i * Rows + j);
            if (i + 1u < Rows)
            {
                result(i + 1u, j) = number * static_cast<Type>(0.001);
            }
            sum += result(i, j);
        }

        result(i, i) += sum;
    }

    return result;
}


template <class Type, std::size_t Rows, std::size_t Columns>
constexpr matrix<Type, Rows, Columns> generate_matrix_simmetrical() noexcept
{
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

    auto result = generate_matrix<Type, Rows, Columns>();
    result = result * result.transpose();

    return result;
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_ITERATIONS

TEST_METHOD(jacobi_test_solve)
{
    constexpr auto vec_x = jacobi_solve(mat_A, vec_b, kEps);
    std::cout << "Matrix solved with Jacobi (include eps = " << kEps << "):\n\n";
    std::cout << "Matrix x:\n" << vec_x.first << "\n\n";
    std::cout << "Number of iterations Jacobi:\n" << vec_x.second;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking Jacobi solve (include eps = " << kEps << "):\n\n";
    constexpr auto jacobi_solve_check = (mat_A * vec_x.first) - vec_b;
    std::cout << "Result of (A * x - b):\n" << jacobi_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}

TEST_METHOD(seidel_test_solve)
{
    constexpr auto vec_x = seidel_solve(mat_A, vec_b, kEps);
    std::cout << "Matrix solved with Seidel (include eps = " << kEps << "):\n\n";
    std::cout << "Matrix x:\n" << vec_x.first << "\n\n";
    std::cout << "Number of iterations Seidel:\n" << vec_x.second;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking Seidel solve (include eps = " << kEps << "):\n\n";
    constexpr auto seidel_solve_check = (mat_A * vec_x.first) - vec_b;
    std::cout << "Result of (A * x - b):\n" << seidel_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_ITERATIONS

#if ENABLE_TESTS_ITERATIONS_DIAGONAL

TEST_METHOD(gen_test_diagonal_predominance)
{
    cx_loops::static_for_s<50u>([] (const auto Index)
    {
        if constexpr (Index > 2u)
        {
            constexpr auto A = generate_matrix_with_diagonal_predominance<double, Index, Index>();
            constexpr auto b = generate_matrix<double, Index, 1u>();
            constexpr auto x1 = jacobi_solve(A, b, kEps);
            constexpr auto x2 = seidel_solve(A, b, kEps);
            std::cout << "Matrix A" << Index << ":\n" << A << "\n\n";
            std::cout << "Matrix b" << Index << ":\n" << b << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x1.first << "\n\n";
            std::cout << "Number of iterations Jacobi" << Index << ":\n" << x1.second << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x2.first << "\n\n";
            std::cout << "Number of iterations Seidel" << Index << ":\n" << x2.second;
            std::cout << "\n\n------------------------------------\n\n";

            if constexpr (Index % 10u == 0u)
            {
                utils::pause_clear("\nPress Enter to continue ELIMINATION...");
            }            
        }
        else
        {
            using namespace std::chrono_literals;
            std::cout << "Initiate ELIMINATION " << (3u - Index) << "...\n";
            std::this_thread::sleep_for(1s);
        }
    });

    cx_loops::static_for_s<75u>([] (const auto Index)
    {
        if constexpr (Index > 50u)
        {
            constexpr auto A = generate_matrix_with_diagonal_predominance<double, Index, Index>();
            constexpr auto b = generate_matrix<double, Index, 1u>();
            constexpr auto x1 = jacobi_solve(A, b, kEps);
            constexpr auto x2 = seidel_solve(A, b, kEps);
            std::cout << "Matrix A" << Index << ":\n" << A << "\n\n";
            std::cout << "Matrix b" << Index << ":\n" << b << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x1.first << "\n\n";
            std::cout << "Number of iterations Jacobi" << Index << ":\n" << x1.second << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x2.first << "\n\n";
            std::cout << "Number of iterations Seidel" << Index << ":\n" << x2.second;
            std::cout << "\n\n------------------------------------\n\n";

            if constexpr (Index % 10u == 0u)
            {
                utils::pause_clear("\nPress Enter to continue ELIMINATION...");
            }
        }
        else
        {
            std::cout << "Skip weaklings " << (50u - Index) << "...\n";
        }
    });

    cx_loops::static_for_s<100u>([] (const auto Index)
    {
        if constexpr (Index > 75u)
        {
            constexpr auto A = generate_matrix_with_diagonal_predominance<double, Index, Index>();
            constexpr auto b = generate_matrix<double, Index, 1u>();
            constexpr auto x1 = jacobi_solve(A, b, kEps);
            constexpr auto x2 = seidel_solve(A, b, kEps);
            std::cout << "Matrix A" << Index << ":\n" << A << "\n\n";
            std::cout << "Matrix b" << Index << ":\n" << b << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x1.first << "\n\n";
            std::cout << "Number of iterations Jacobi" << Index << ":\n" << x1.second << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x2.first << "\n\n";
            std::cout << "Number of iterations Seidel" << Index << ":\n" << x2.second;
            std::cout << "\n\n------------------------------------\n\n";

            if constexpr (Index % 10u == 0u && Index != 100u)
            {
                utils::pause_clear("\nPress Enter to continue ELIMINATION...");
            }
        }
        else
        {
            std::cout << "Kill enemies " << (75u - Index) << "...\n";
        }
    });
}

#endif // ENABLE_TESTS_ITERATIONS_DIAGONAL

#if ENABLE_TESTS_ITERATONS_SIMMETRICAL

TEST_METHOD(gen_test_simmetrical)
{
    cx_loops::static_for_s<50u>([] (const auto Index)
    {
        if constexpr (Index > 2u)
        {
            constexpr auto A = generate_matrix_simmetrical<double, Index, Index>();
            constexpr auto b = generate_matrix<double, Index, 1u>();
            constexpr auto x1 = jacobi_solve(A, b, kEps);
            constexpr auto x2 = seidel_solve(A, b, kEps);
            std::cout << "Matrix A" << Index << ":\n" << A << "\n\n";
            std::cout << "Matrix b" << Index << ":\n" << b << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x1.first << "\n\n";
            std::cout << "Number of iterations Jacobi" << Index << ":\n" << x1.second << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x2.first << "\n\n";
            std::cout << "Number of iterations Seidel" << Index << ":\n" << x2.second;
            std::cout << "\n\n------------------------------------\n\n";

            if constexpr (Index % 10u == 0u)
            {
                utils::pause_clear("\nPress Enter to continue ELIMINATION...");
            }            
        }
        else
        {
            using namespace std::chrono_literals;
            std::cout << "Initiate ELIMINATION " << (3u - Index) << "...\n";
            std::this_thread::sleep_for(1s);
        }
    });

    cx_loops::static_for_s<75u>([] (const auto Index)
    {
        if constexpr (Index > 50u)
        {
            constexpr auto A = generate_matrix_simmetrical<double, Index, Index>();
            constexpr auto b = generate_matrix<double, Index, 1u>();
            constexpr auto x1 = jacobi_solve(A, b, kEps);
            constexpr auto x2 = seidel_solve(A, b, kEps);
            std::cout << "Matrix A" << Index << ":\n" << A << "\n\n";
            std::cout << "Matrix b" << Index << ":\n" << b << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x1.first << "\n\n";
            std::cout << "Number of iterations Jacobi" << Index << ":\n" << x1.second << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x2.first << "\n\n";
            std::cout << "Number of iterations Seidel" << Index << ":\n" << x2.second;
            std::cout << "\n\n------------------------------------\n\n";

            if constexpr (Index % 10u == 0u)
            {
                utils::pause_clear("\nPress Enter to continue ELIMINATION...");
            }
        }
        else
        {
            std::cout << "Skip weaklings " << (50u - Index) << "...\n";
        }
    });

    cx_loops::static_for_s<100u>([] (const auto Index)
    {
        if constexpr (Index > 75u)
        {
            constexpr auto A = generate_matrix_simmetrical<double, Index, Index>();
            constexpr auto b = generate_matrix<double, Index, 1u>();
            constexpr auto x1 = jacobi_solve(A, b, kEps);
            constexpr auto x2 = seidel_solve(A, b, kEps);
            std::cout << "Matrix A" << Index << ":\n" << A << "\n\n";
            std::cout << "Matrix b" << Index << ":\n" << b << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x1.first << "\n\n";
            std::cout << "Number of iterations Jacobi" << Index << ":\n" << x1.second << "\n\n";
            std::cout << "Matrix x" << Index << ":\n" << x2.first << "\n\n";
            std::cout << "Number of iterations Seidel" << Index << ":\n" << x2.second;
            std::cout << "\n\n------------------------------------\n\n";

            if constexpr (Index % 10u == 0u && Index != 100u)
            {
                utils::pause_clear("\nPress Enter to continue ELIMINATION...");
            }
        }
        else
        {
            std::cout << "Kill enemies " << (75u - Index) << "...\n";
        }
    });
}

#endif // ENABLE_TESTS_ITERATONS_SIMMETRICAL

} // namespace vv
