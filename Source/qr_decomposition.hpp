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
template <class Type, std::size_t Rows, std::size_t Columns>
constexpr std::array<matrix<Type, Rows, Columns>, 2u>
    qr_decompose(const matrix<Type, Rows, Columns>& mat, const Type eps = kDefault_eps<Type>)
{
    // Alias the size_type.
    using size_type = typename matrix<Type, Rows, Columns>::size_type;

    // Create identity matrix.
    auto I = matrix<Type, Rows, Rows>::create_identity();
	matrix<Type, Rows, Rows> P = I;
	matrix<Type, Rows, Rows> Q = I;
    
	matrix<Type, Rows, Columns> R = mat;

    // Using Householder reflections.
	for (size_type i = 0u; i < Columns; ++i)
    {
		matrix<Type, Rows, 1u> u{};
        matrix<Type, Rows, 1u> v{};
		
		Type mag{};
		for (size_type j = i; j < Rows; ++j)
        {
			u(j, 0u) = R(j, i); // [j * Columns + i] 
			mag += u(j, 0u) * u(j, 0u);
		}
		mag = cx::sqrt(mag);
		
		Type alpha = u(i, 0u) < static_cast<Type>(0) ? mag : mag;

		mag = static_cast<Type>(0);
		for (size_type j = i; j < Rows; ++j)
        {
			v(j, 0u) = (j == i ? u(j, 0u) + alpha : u(j, 0u));
			mag += v(j, 0u) * v(j, 0u);
		}
		mag = cx::sqrt(mag);

		if (mag < eps) continue;
        

		for (size_type j = i; j < Rows; ++j)
        {
            v(j, 0u) /= mag;
        }

		P = I - (v * v.transpose()) * static_cast<Type>(2);

		R = P * R;
		Q = Q * P;

	} // for (size_type i = 0u; i < Columns; ++i)
    
    return { Q, R };
}


/// ===== TEST SECTION =====
#if ENABLE_TESTS_QR

TEST_METHOD(qr_test_decomposition)
{
    constexpr auto qr_value = qr_decompose(mat_A);
    std::cout << "QR-decompose:\n\n";
    std::cout << "Matrix A:\n" << mat_A << "\n\n";
    std::cout << "Matrix Q:\n" << qr_value.at(0u) << "\n\n";
    std::cout << "Matrix R:\n" << qr_value.at(1u);
    std::cout << "\n\n------------------------------------\n\n";

    std::cout << "Checking QR-decompose:\n\n";
    constexpr auto qr_product = qr_value.at(0u) * qr_value.at(1u);
    std::cout << "Matrix (Q * R):\n" << qr_product;
    std::cout << "\n\n------------------------------------\n\n";

    // std::cout << "QR-solve:\n\n";
    // constexpr auto rx = qr_value.at(1u);
    // constexpr auto qTb = qr_value.at(0u).transpose() * vec_b;
    // //constexpr auto vecX_qr = vv::matrix<double>::solve(rx, qTb, kEps);
    // constexpr auto vecX_qr = vv::matrix<double>::backward_substitution(rx, qTb, kEps);
    // std::cout << "Matrix solved with Gaussian (include eps):\n\n";
    // std::cout << "Matrix (QT * b):\n" << qTb << "\n\n";
    // std::cout << "Matrix x:\n" << vecX_qr;
    // std::cout << "\n\n------------------------------------\n\n";

    // std::cout << "Checking QR-solve:\n\n";
    // constexpr auto qr_solve_check = (rx * vecX_qr) - qTb;
    // std::cout << "Result of (R * x - QT * b):\n" << qr_solve_check;
    // std::cout << "\n\n------------------------------------\n\n";
}

#endif // ENABLE_TESTS_QR

} // namespace vv
