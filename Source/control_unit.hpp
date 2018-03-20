#pragma once

#include <iostream>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>

#include "matrix.hpp"


namespace vv
{

/// ===== MACROS SECTION =====
#define TRY_BLOCK(expr) try { expr } \
                        catch(std::exception& ex) { std::cout << ex.what() << " <==\n\n"; }

#define TEST_METHOD_RUNTIME(name) template <class Type, std::size_t Rows, std::size_t Columns> \
                        void name(const matrix<Type, Rows, Columns>& mat_A, \
                                  const matrix<Type, Rows, 1u>& vec_b = matrix<Type, Rows, 1u>{}, \
                                  const Type eps = kDefault_eps<Type>)

#define TEST_METHOD(name) void name()


// Test methods declarations.
TEST_METHOD(general_test_condition_number);

TEST_METHOD(gaussian_test_for_matrix_N);
TEST_METHOD(gaussian_test_for_matrix_eps_N);
TEST_METHOD(gaussian_test_inverse);

TEST_METHOD(lu_test_major_ver);
TEST_METHOD(lu_test_minor_ver);
TEST_METHOD(lu_test_determenant);
TEST_METHOD(lu_test_rank);

TEST_METHOD(lup_test_major_ver);
TEST_METHOD(lup_test_solve);

TEST_METHOD(lup_test_minor_ver);
TEST_METHOD(lup_test_determenant);
TEST_METHOD(lup_test_inverse);

TEST_METHOD(lupq_test_decompose);
TEST_METHOD(lupq_test_solve);

TEST_METHOD(qr_test_decomposition);

TEST_METHOD(jacobi_test_solve);
TEST_METHOD(seidel_test_solve);

TEST_METHOD(gen_test_diagonal_predominance);
TEST_METHOD(gen_test_simmetrical);


/// ===== CONSTANT SECTION =====
// ===== DATA 0 (base constants) =====
template <class Type>
constexpr Type kDefault_eps = static_cast<Type>(0.000001);

constexpr int N = 12;
constexpr double kEps = 0.0000000000001;

// ===== DATA 1 (base tests) =====
constexpr vv::matrix<double, 3u, 3u> data_N{ N + 2, 1,     1,
                                             1,     N + 4, 1,
                                             1,     1,     N + 6 };

constexpr vv::matrix<double, 3u, 1u> vec_b_N{ N + 4, N + 6, N + 8 };

// ===== DATA 2 (bad case) =====
constexpr vv::matrix<double, 3u, 3u> data{ 1, -1, -1,
                                           0,  1, -1,
                                           0,  0,  1  };

constexpr vv::matrix<double, 3u, 3u> eps_mat{ kEps * N, -kEps * N, -kEps * N,
                                              kEps * N,  kEps * N, -kEps * N,
                                              kEps * N,  kEps * N,  kEps * N  };
constexpr auto data_eps = data + eps_mat;

constexpr vv::matrix<double, 3u, 1u> vec_b_eps{ -1, -1, -1 };

// ===== DATA 3 (singular matrix, no solutions) =====
constexpr vv::matrix<double, 4u, 4u> data_singular{  2,  -1, 4, 7,
                                                     4,  -2, 8, 5,
                                                    -2,   1, 6, 9,
                                                     10, -5, 1, 3  };

constexpr vv::matrix<double, 4u, 1u> vec_b_singular{ 4, 8, 6, 1 };

// ===== DATA 4 (from Valera, no solutions) =====
constexpr vv::matrix<double, 3u, 5u> data_val{ 1, 2, -1, 1, 1,
                                               2, 4, -4, 3, 1,
                                               1, 2,  1, 2, 3 };

constexpr vv::matrix<double, 3u, 1u> vec_b_val{ 1, 0, 2 };

// ===== DATA 5 (for rank test) =====
constexpr vv::matrix<double, 4u, 4u> data_rank{  2.0,  1.0, 11.0,  2.0,
                                                 1.0,  0.0,  4.0, -1.0,
                                                11.0,  4.0, 56.0,  5.0,
                                                 2.0, -1.0,  5.0, -6.0  };

// ===== DATA 6 (not square matrix, no solutions) =====
constexpr vv::matrix<double, 5u, 3u> data_not_square_singular{   2.0,  1.0, 11.0,
                                                                 2.0,  1.0,  1.0,
                                                                 0.0,  4.0, -1.0,
                                                                 3.0, 11.0,  4.0,
                                                                56.0,  5.0,  5.0  };

constexpr vv::matrix<double, 5u, 1u> vec_b_not_square_singular{ 1.0, 0.0, 2.0, 3.0, 7.0 };

// ===== DATA 7 (not square matrix) =====
constexpr vv::matrix<double, 5u, 3u> data_not_square_row{  3.0,  1.0,  1.0,
                                                           2.0,  2.0,  1.0,
                                                           7.0,  5.0,  2.0,
                                                           1.0, 11.0,  7.0,
                                                          19.0, 17.0, 16.0  };

constexpr vv::matrix<double, 5u, 1u> vec_b_not_square_row{ 5.0, 5.0, 14.0, 19.0, 52.0 };

// ===== DATA 8 (not square matrix) =====
constexpr vv::matrix<double, 3u, 5u> data_not_square_col{  3.0, 1.0, 1.0, 4.0, 3.0,
                                                           2.0, 2.0, 1.0, 5.0, 4.0,
                                                           7.0, 5.0, 2.0, 1.0, 9.0  };

constexpr vv::matrix<double, 3u, 1u> vec_b_not_square_col{ 12.0, 14.0, 24.0 };


#define mat_A data_eps
#define vec_b vec_b_eps

#define ENABLE_TESTS_GENERAL               1
#define ENABLE_TESTS_GAUSSIAN              1
#define ENABLE_TESTS_LU                    1
#define ENABLE_TESTS_LUP                   1
#define ENABLE_TESTS_LUP_ONLY_SQUARE       1
#define ENABLE_TESTS_LUPQ                  1
#define ENABLE_TESTS_QR                    1
#define ENABLE_TESTS_ITERATIONS            1
#define ENABLE_TESTS_ITERATIONS_DIAGONAL   1
#define ENABLE_TESTS_ITERATONS_SIMMETRICAL 1


void start_tests()
{
#if ENABLE_TESTS_GENERAL
    general_test_condition_number();
#endif // ENABLE_TESTS_GENERAL


#if ENABLE_TESTS_GAUSSIAN
    gaussian_test_for_matrix_N();
    gaussian_test_for_matrix_eps_N();
    gaussian_test_inverse();
#endif // ENABLE_TESTS_GAUSSIAN


#if ENABLE_TESTS_LU
    lu_test_major_ver();
    lu_test_minor_ver();
    lu_test_determenant();
    lu_test_rank();
#endif // ENABLE_TESTS_LU


#if ENABLE_TESTS_LUP
    lup_test_major_ver();
    lup_test_solve();
#endif // ENABLE_TESTS_LUP

#if ENABLE_TESTS_LUP_ONLY_SQUARE
    lup_test_minor_ver();
    lup_test_determenant();
    lup_test_inverse();
#endif // ENABLE_TESTS_LUP_ONLY_SQUARE


#if ENABLE_TESTS_LUPQ
    lupq_test_decompose();
    lupq_test_solve();
#endif // ENABLE_TESTS_LUPQ


#if ENABLE_TESTS_QR
    qr_test_decomposition();
#endif // ENABLE_TESTS_QR


#if ENABLE_TESTS_ITERATIONS
    jacobi_test_solve();
    seidel_test_solve();
#endif // ENABLE_TESTS_ITERATIONS


#if ENABLE_TESTS_ITERATIONS_DIAGONAL
    gen_test_diagonal_predominance();
#endif // ENABLE_TESTS_ITERATIONS_DIAGONAL


#if ENABLE_TESTS_ITERATONS_SIMMETRICAL
    gen_test_simmetrical();
#endif // ENABLE_TESTS_ITERATONS_SIMMETRICAL
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_GENERAL

TEST_METHOD(general_test_condition_number)
{
    constexpr auto condition_number = mat_A.calculate_condition_number();
    std::cout << "Calculated condition number:\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Condition number of matrix A:\n" << general_condition_number;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_GENERAL

#if ENABLE_TESTS_GAUSSIAN

TEST_METHOD(gaussian_test_for_matrix_N)
{
    constexpr auto vec_x = matrix<double>::solve(mat_A, vec_b);
    std::cout << "Matrix solved with Gaussian:\n\n";
    std::cout << "Matrix x:\n" << vec_x;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking Gaussian solve:\n\n";
    constexpr auto gaussian_check = (mat_A * vec_x) - vec_b;
    std::cout << "Result of (A * x - b):\n" << gaussian_check;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(gaussian_test_for_matrix_eps_N)
{
    constexpr auto vec_x = matrix<double>::solve(mat_A, vec_b, kEps);
    std::cout << "Matrix solved with Gaussian (include kEps = " << kEps <<  ")):\n\n";
    std::cout << "Matrix x:\n" << vec_x;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking Gaussian solve (include kEps = " << kEps <<  ")):\n\n";
    constexpr auto gaussian_solve_check = (mat_A * vec_x) - vec_b;
    std::cout << "Result of (Aeps * x - b):\n" << gaussian_solve_check;
    std::cout << "\n\n------------------------------------\n\n";
}


TEST_METHOD(gaussian_test_inverse)
{
    constexpr auto casual_invert = mat_A.inverse();
    std::cout << "Invert by matrix method (Gaussian):\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Inverted matrix to A:\n" << casual_invert;
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking inverse:\n\n";
    constexpr auto inverse_check1 = mat_A * casual_invert;
    std::cout << "Result of (A * IA):\n" << inverse_check1 << "\n\n";

    constexpr auto inverse_check2 = casual_invert * mat_A;    
    std::cout << "Result of (IA * A):\n" << inverse_check2;
    std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_GAUSSIAN

} // namespace vv
