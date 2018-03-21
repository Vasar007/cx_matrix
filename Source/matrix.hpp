#pragma once

#include <iostream>
#include <iterator>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <type_traits>
#include <initializer_list>
#include <string>

#include "cx_math.h"


/*
 *
 * As base for matrix class using https://github.com/akalicki/matrix
 *
 */

namespace vv::detail
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


template<class InputIt, class T, class BinaryOperation>
constexpr T accumulate(InputIt first, InputIt last, T init, BinaryOperation op)
{
    for (; first != last; ++first)
    {
        init = op(init, *first);
    }
    return init;
}


template<class T, class U = T>
constexpr T exchange(T& obj, U&& new_value)
{
    T old_value = std::move(obj);
    obj = std::forward<U>(new_value);
    return old_value;
}


template<class T>
constexpr void swap(T& a, T& b) noexcept(std::is_nothrow_move_constructible_v<T>
                                         && std::is_nothrow_move_assignable_v<T>)
{
    T temp = std::move(a);
    a = std::move(b);
    b = std::move(temp);
}

} // namespace vv::detail


namespace vv
{

/// ===== MATRIX IMLEMENTATION SECTION =====
template <class _Type = int, std::size_t _Rows = 1u, std::size_t _Columns = 1u>
class matrix
{
public:
    using value_type                    = _Type;
	using size_type                     = std::size_t;
	using difference_type               = std::ptrdiff_t;

    template <class Type = value_type, std::size_t N = _Rows>
    using container                     = std::array<Type, N>;
    using row_container                 = container<value_type, _Columns>;
    using row_container_reference       = container<value_type, _Columns>&;    
    using const_row_container_reference = const container<value_type, _Columns>&;

	using pointer                       = value_type*;
	using const_pointer                 = const value_type*;

	using reference                     = value_type&;
	using const_reference               = const value_type&;

    using std_matrix                    = matrix<value_type, _Rows, _Columns>;


    static constexpr value_type EPS = static_cast<value_type>(0.0000000001);

    static_assert(std::is_arithmetic_v<value_type>, "Matrix elements type has to be arithmetic!");
    static_assert(_Rows > 0u && _Columns > 0u, "Incorrect size parameters!");

    constexpr matrix() = default;

    constexpr matrix(std::initializer_list<value_type> list)
    : _data()
    {
        size_type row_counter = 0u;
        size_type col_counter = 0u;
        for (const auto elem : list)
        {
            _data.at(row_counter).at(col_counter) = elem;
            ++col_counter;
            if (row_counter == _Rows && col_counter == _Columns)
            {
                break;
            }
            else if (col_counter == _Columns)
            {
                col_counter = 0u;
                ++row_counter;
            }
        }
    }

    ~matrix() noexcept = default;
    constexpr matrix(const std_matrix& other) = default;
    constexpr std_matrix& operator=(const std_matrix& other) = default;
    constexpr matrix(std_matrix&& other) noexcept = default;
    constexpr std_matrix& operator=(std_matrix&& other) noexcept = default;


    std::string get_dimension() const noexcept
    {
        return std::to_string(get_rows_number()) + std::string("x")
                + std::to_string(get_columns_number());
    }


    constexpr const container<container<value_type, _Columns>, _Rows>& get_data() const noexcept
    {
        return _data;
    }


    constexpr size_type get_rows_number() const noexcept
    {
        return _Rows;
    }


    constexpr size_type get_columns_number() const noexcept
    {
        return _Columns;
    }


    constexpr const_reference operator()(const size_type i, const size_type j) const
    { 
         return (i >= _Rows || j >= _Columns ? throw std::invalid_argument("Out of range!")
                                             : _data.at(i).at(j));
        //return _data.at(i).at(j);
    }


    constexpr reference operator()(const size_type i, const size_type j)
    { 
        return (i >= _Rows || j >= _Columns ? throw std::invalid_argument("Out of range!")
                                            : _data.at(i).at(j));
        //return _data.at(i).at(j);
    }


    constexpr const_row_container_reference& operator()(const size_type i) const
    { 
         return (i >= _Rows ? throw std::invalid_argument("Out of range!") : _data.at(i));
        //return _data.at(i);
    }


    constexpr row_container_reference& operator()(const size_type i)
    { 
        return (i >= _Rows ? throw std::invalid_argument("Out of range!") : _data.at(i));
        //return _data.at(i);
    }


    constexpr auto& operator+=(const std_matrix& rhs) noexcept
    {
        for (size_type i = 0u; i < _Rows; ++i)
        {
            for (size_type j = 0u; j < _Columns; ++j)
            {
                _data.at(i).at(j) += rhs(i, j);
            }
        }
        return *this;
    }
    

    constexpr auto& operator-=(const std_matrix& rhs) noexcept
    {
        for (size_type i = 0u; i < _Rows; ++i)
        {
            for (size_type j = 0u; j < _Columns; ++j)
            {
                _data.at(i).at(j) -= rhs(i, j);
            }
        }
        return *this;
    }
    

    constexpr auto& operator*=(const std_matrix& rhs) noexcept
    {
        const auto transp = rhs.transpose();
        std_matrix temp{};
        for (size_type i = 0u; i < _Rows; ++i)
        {
            for (size_type j = 0u; j < _Columns; ++j)
            {
                for (size_type k = 0u; k < _Columns; ++k)
                {
                    temp(i, j) += (_data.at(i).at(k) * transp(j, k));
                }
            }
        }
        return (*this = temp);
    }
    

    constexpr auto& operator*=(const value_type num) noexcept
    {
        for (auto& row : _data)
        {
            for (auto& elem : row)
            {
                elem *= num;
            }
        }
        return *this;
    }
    

    constexpr auto& operator/=(const value_type num) noexcept
    {
        for (auto& row : _data)
        {
            for (auto& elem : row)
            {
                elem /= num;
            }
        }
        return *this;
    }


    constexpr auto operator^(const int num) const noexcept
    {
        std_matrix temp(*this);
        return _exp_helper(temp, num);
    }


    constexpr void swap_rows(const size_type row_1, const size_type row_2)
    {
        detail::swap(_data.at(row_1), _data.at(row_2));
    }


    constexpr void swap_columns(const size_type column_1, const size_type column_2)
    {
	    for (size_type i = 0; i < _Rows; i++)
	    {
	    	detail::swap(_data.at(i).at(column_1), _data.at(i).at(column_2));
	    }
    }
    

    constexpr auto transpose() const noexcept
    {
        matrix<value_type, _Columns, _Rows> transp{};
        for (size_type i = 0u; i < _Rows; ++i)
        {
            for (size_type j = 0u; j < _Columns; ++j)
            {
                transp(j, i) = _data.at(i).at(j);
            }
        }
        return transp;
    }


    constexpr bool is_empty(const size_type start_row_index = 0u,
                            const size_type start_column_index = 0u,
                            const value_type eps = EPS) const noexcept
    {
        if (start_row_index > _Rows || start_column_index > _Columns)
        {
            return false;
        }

        for (size_type i = start_row_index; i < _Rows; ++i)
        {
	    	for (size_type j = start_column_index; j < _Columns; ++j)
	    	{
	    		if (cx::abs(_data.at(i).at(j)) >= eps)
	    		{
	    			return false;
	    		}
	    	}
        }

        return true;        
    }


    constexpr std::pair<size_type, size_type>
        find_max_element(const size_type start_row_index = 0u,
                         const size_type start_column_index = 0u) const noexcept
    {
        if (start_row_index > _Rows || start_column_index > _Columns)
        {
            return { start_row_index, start_column_index };
        }

        value_type abs_max{};
	    size_type row = start_row_index;
	    size_type column = start_column_index;

	    for (size_type i = start_row_index; i < _Rows; ++i)
        {
	    	for (size_type j = start_column_index; j < _Columns; ++j)
	    	{
	    		if (cx::abs(_data.at(i).at(j)) > abs_max)
	    		{
	    			abs_max = cx::abs(_data.at(i).at(j));
	    			row = i;
	    			column = j;
	    		}
	    	}
        }
	    return { row, column };
    }


    constexpr value_type calculate_condition_number() const noexcept
    {
        constexpr auto abs_plus = [](const value_type& a, const value_type& b)
        {
            return cx::abs(a) + cx::abs(b);
        };

        value_type condition_number{};
        for (const auto& row: _data)
        {
            const auto summ_in_row = detail::accumulate(std::begin(row), std::end(row),
                                                        value_type{}, abs_plus);
            condition_number = std::max(condition_number, summ_in_row);
        }

        return condition_number;
    }


    template <class T = value_type, std::size_t Size = _Rows>
    static constexpr auto create_identity() noexcept
    {
        matrix<T, Size, Size> temp{};
        for (size_type i = 0u; i < Size; ++i)
        {
            for (size_type j = 0u; j < Size; ++j)
            {
                if (i == j)
                {
                    temp(i, j) = static_cast<T>(1);
                }
                else
                {
                    temp(i, j) = static_cast<T>(0);
                }
            }
        }
        return temp;
    }


    template <std::size_t _Rows_, std::size_t _Columns_A, std::size_t _Columns_b>
    static constexpr std::pair<matrix<value_type, _Rows_, _Columns_A>,
                               matrix<value_type, _Rows_, _Columns_b>>
        forward_substitution(matrix<value_type, _Rows_, _Columns_A> A,
                             matrix<value_type, _Rows_, _Columns_b> b,
                             const value_type eps = EPS)
    {
        static_assert(_Columns_b == 1u, "Matrix contains more than one columns in right "
                                        "hand side vector!");

        // Gaussian elimination.
        for (size_type i = 0u; i < _Rows_; ++i)
        {
            if (cx::abs(A(i, i)) < eps)
            {
                // Pivot 0 - throw error.
                throw std::domain_error("Error: the coefficient matrix has 0 as a pivot.");
            }

            for (size_type j = i + 1u; j < _Rows_; ++j)
            {
                for (size_type k = i + 1u; k < _Columns_A; ++k)
                {
                    A(j, k) -= A(i, k) * (A(j, i) / A(i, i));
                    if (cx::abs(A(j, k)) < eps)
                    {
                        A(j, k) = static_cast<value_type>(0);
                    }
                }
                b(j, 0u) -= b(i, 0u) * (A(j, i) / A(i, i));
                if (cx::abs(A(j, 0u)) < eps)
                {
                    A(j, 0u) = static_cast<value_type>(0);
                }
                A(j, i) = static_cast<value_type>(0);
            }
        } // for (size_type i = 0u; i < _Rows_; ++i)

        return { A, b };
    }


    template <std::size_t _Rows_, std::size_t _Columns_A, std::size_t _Columns_b>
    static constexpr matrix<value_type, _Rows_, _Columns_b>
        backward_substitution(const matrix<value_type, _Rows_, _Columns_A>& A,
                              const matrix<value_type, _Rows_, _Columns_b>& b,
                              const value_type eps = EPS)
    {
        static_assert(_Columns_b == 1u, "Matrix contains more than one columns in right "
                                        "hand side vector!");

        // Back substitution.
        matrix<value_type, _Rows_, _Columns_b> x{}; // _Columns_b == 1u
        for (size_type i = _Rows_ - 1u; ; --i)
        {
            value_type sum{};
            for (size_type j = i + 1u; j < _Rows_; ++j)
            {
                sum += A(i, j) * x(j, 0u);
            }
            x(i, 0u) = (b(i, 0u) - sum) / A(i, i);
            if (cx::abs(x(i, 0u)) < eps)
            {
                x(i, 0u) = static_cast<value_type>(0);
            }

            if (i == 0u)
            {
                break;
            }
        } // for (size_type i = _Rows_ - 1u; ; --i)
        
        return x;
    }

    template <std::size_t _Rows_, std::size_t _Columns_A, std::size_t _Columns_b>
    static constexpr auto solve(const matrix<value_type, _Rows_, _Columns_A>& A,
                                const matrix<value_type, _Rows_, _Columns_b>& b,
                                const value_type eps = EPS)
    {
        static_assert(_Columns_b == 1u, "Matrix contains more than one columns in right "
                                        "hand side vector!");

        const auto result_of_forw_subs = forward_substitution(A, b, eps);

        return backward_substitution(result_of_forw_subs.first, result_of_forw_subs.second);
    }


    template <std::size_t _Rows_A, std::size_t _Columns_A,
              std::size_t _Rows_b, std::size_t _Columns_b>
    static constexpr auto band_solve(matrix<value_type, _Rows_A, _Columns_A> A,
                                     matrix<value_type, _Rows_b, _Columns_b> b,
                                     const value_type coef)
    {
        // Optimized Gaussian elimination.
        value_type bandsBelow = static_cast<value_type>((coef - 1) / 2);
        for (size_type i = 0u; i < _Rows_A; ++i)
        {
            if (A(i, i) == static_cast<value_type>(0))
            {
                // Pivot 0 - throw exception.
                throw std::domain_error("Error: the coefficient matrix has 0 as a "
                                        "pivot. Please fix the input and try again.");
            }
            for (size_type j = i + 1; j < _Rows_A && j <= i + bandsBelow; ++j)
            {
                size_type k = i + 1;
                while (k < _Columns_A && A(j, k))
                {
                    A(j, k) -= A(i, k) * (A(j, i) / A(i, i));
                    ++k;
                }
                b(j, 0u) -= b(i, 0u) * (A(j, i) / A(i, i));
                A(j, i) = static_cast<value_type>(0);
            }
        }

        // Back substitution.
        matrix<value_type, _Rows_b, 1u> x{};
        x(_Rows_b - 1u, 0u) = b(_Rows_b - 1u, 0u) / A(_Rows_b - 1u, _Rows_b - 1u);
        for (size_type i = _Rows_b - 2u; ; --i)
        {
            value_type sum{};
            for (size_type j = i + 1u; j < _Rows_b; ++j)
            {
                sum += A(i, j) * x(j, 0u);
            }
            x(i, 0u) = (b(i, 0u) - sum) / A(i, i);

            if (i == 0u)
            {
                break;
            }
        }

        return x;
    }


    // Functions on vectors.
    template <std::size_t _Rows_, std::size_t _Columns_A, std::size_t _Columns_B>
    static constexpr value_type dot_product(const matrix<value_type, _Rows_, _Columns_A>& A,
                                            const matrix<value_type, _Rows_, _Columns_B>& B,
                                            const size_type column) noexcept
    {
        value_type sum{};
        for (size_type i = 0u; i < _Rows_; ++i)
        {
            sum += (A(i, column) * B(i, column));
        }
        return sum;
    }


    // Functions on augmented matrices.
    template <std::size_t _Rows_A, std::size_t _Columns_A,
              std::size_t _Rows_B, std::size_t _Columns_B>
    static constexpr auto augment(const matrix<value_type, _Rows_A, _Columns_A>& A,
                                  const matrix<value_type, _Rows_B, _Columns_B>& B) noexcept
    {
        matrix<value_type, _Rows_A, _Columns_A + _Columns_B> AB{};
        for (size_type i = 0u; i < _Rows_A; ++i)
        {
            for (size_type j = 0u; j < _Columns_A + _Columns_B; ++j)
            {
                if (j < _Columns_A)
                {
                    AB(i, j) = A(i, j);
                }
                else
                {
                    AB(i, j) = B(i, j - _Columns_B);
                }
            }
        }
        return AB;
    }

    constexpr auto gaussian_eliminate(const value_type eps = EPS) const
    {
        std_matrix Ab(*this);
        int rows = _Rows;
        int cols = _Columns;
        int Acols = cols - 1;

        int i = 0; // Row tracker.
        int j = 0; // Column tracker.

        // Iterate through the rows.
        while (i < rows)
        {
            // Find a pivot for the row.
            bool pivot_found = false;
            while (j < Acols && !pivot_found)
            {
                if (Ab(i, j) != static_cast<value_type>(0))
                { // Pivot not equal to 0.
                    pivot_found = true;
                }
                else
                { // Check for a possible swap.
                    int max_row = i;
                    value_type max_val{};
                    for (int k = i + 1; k < rows; ++k)
                    {
                        value_type cur_abs = cx::abs(Ab(k, j));
                        if (cur_abs > max_val)
                        {
                            max_row = k;
                            max_val = cur_abs;
                        }
                    }
                    if (max_row != i)
                    {
                        Ab.swap_rows(max_row, i);
                        pivot_found = true;
                    }
                    else
                    {
                        ++j;
                    }
                }
            }

            // Perform elimination as normal if pivot was found.
            if (pivot_found)
            {
                for (int t = i + 1; t < rows; ++t)
                {
                    for (int s = j + 1; s < cols; ++s)
                    {
                        Ab(t, s) = Ab(t, s) - Ab(i, s) * (Ab(t, j) / Ab(i, j));
                        if (Ab(t, s) < eps && Ab(t, s) > -1 * eps)
                        {
                            Ab(t, s) = static_cast<value_type>(0);
                        }
                    }
                    Ab(t, j) = static_cast<value_type>(0);
                }
            }

            ++i;
            ++j;
        }

        return Ab;
    }

    constexpr auto row_reduce_from_gaussian(const value_type eps = EPS) const
    {
        std_matrix result(*this);
        int rows = _Rows;
        int cols = _Columns;

        int i = rows - 1; // Row tracker.
        int j = cols - 1; // Column tracker.

        // Iterate through every row.
        while (i >= 0 && j >= 0)
        {
            // Find the pivot column.
            int k = j - 1;
            while (k >= 0)
            {
                if (result(i, k) != static_cast<value_type>(0))
                {
                    j = k;
                }
                --k;
            }

            // Zero out elements above pivots if pivot not 0.
            if (result(i, j) != static_cast<value_type>(0))
            {
            
                for (int t = i - 1; t >= 0; --t)
                {
                    for (int s = 0; s < cols; ++s)
                    {
                        if (s != j)
                        {
                            result(t, s) = result(t, s) - result(i, s) 
                                            * (result(t, j) / result(i, j));
                            if (result(t, s) < eps && result(t, s) > -1 * eps)
                            {
                                result(t, s) = static_cast<value_type>(0);
                            }
                        }
                    }
                    result(t, j) = static_cast<value_type>(0);
                }

                // Divide row by pivot.
                for (int k = j + 1; k < cols; ++k)
                {
                    result(i, k) = result(i, k) / result(i, j);
                    if (result(i, k) < eps && result(i, k) > -1 * eps)
                    {
                        result(i, k) = static_cast<value_type>(0);
                    }
                }
                result(i, j) = static_cast<value_type>(1);

            }

            --i;
            --j;
        }

        return result;
    }


    // Not tested this method.
    void read_solutions_from_RREF(std::ostream& os)
    {
        std_matrix result(*this);

        // Print number of solutions.
        bool hasSolutions = true;
        bool doneSearching = false;
        size_type i = 0u;
        while (!doneSearching && i < _Rows)
        {
            bool allZeros = true;
            for (size_type j = 0u; j < _Columns - 1u; ++j)
            {
                if (result(i, j) != static_cast<value_type>(0))
                {
                    allZeros = false;
                }
            }
            if (allZeros && result(i, _Columns - 1u) != static_cast<value_type>(0))
            {
                hasSolutions = false;
                os << "NO SOLUTIONS\n\n";
                doneSearching = true;
            }
            else if (allZeros && result(i, _Columns - 1u) == static_cast<value_type>(0))
            {
                os << "INFINITE SOLUTIONS\n\n";
                doneSearching = true;
            }
            else if (_Rows < _Columns - 1u)
            {
                os << "INFINITE SOLUTIONS\n\n";
                doneSearching = true;
            }
            i++;
        }
        if (!doneSearching)
        {
            os << "UNIQUE SOLUTION\n\n";
        }

        // Get solutions if they exist.
        if (hasSolutions)
        {
            matrix<value_type, _Columns - 1u, 1u> particular{};
            matrix<value_type, _Columns - 1u, 1u> special{};

            for (size_type i = 0u; i < _Rows; ++i)
            {
                bool pivotFound = false;
                bool specialCreated = false;
                for (size_type j = 0u; j < _Columns - 1u; ++j)
                {
                    if (result(i, j) != static_cast<value_type>(0))
                    {
                        // If pivot variable, add b to particular.
                        if (!pivotFound)
                        {
                            pivotFound = true;
                            particular(j, 0) = result(i, _Columns - 1u);
                        }
                        else
                        { // Otherwise, add to special solution.
                            if (!specialCreated)
                            {
                                special = matrix<value_type, _Columns - 1u, 1u>{};
                                specialCreated = true;
                            }
                            special(j, 0) = -result(i, j);
                        }
                    }
                }
                os << "Special solution\n" << special << '\n';
            }
            os << "Particular solution:\n" << particular << '\n';
        }
    }


    constexpr std_matrix inverse() const
    {
        auto I = create_identity();
        auto AI = augment(*this, I);
        auto U = AI.gaussian_eliminate();
        auto IAInverse = U.row_reduce_from_gaussian();
        std_matrix AInverse{};
        for (size_type i = 0u; i < _Rows; ++i)
        {
            for (size_type j = 0u; j < _Columns; ++j)
            {
                AInverse(i, j) = IAInverse(i, j + _Columns);
            }
        }
        return AInverse;
    }

    
    constexpr value_type calculate_elements_summ() const noexcept
    {
        value_type summ{};
        for (const auto& row: _data)
        {
            const auto summ_in_row = detail::accumulate(std::begin(row), std::end(row),
                                                        value_type{});
            summ += summ_in_row;
        }

        return summ;
    }


    template <class T = value_type, std::size_t Rows = _Rows, std::size_t Columns = _Columns>
    static constexpr matrix<T, Rows, Columns> get_error_matrix() noexcept
    {
        matrix<T, Rows, Columns> err_matrix{};
        for (size_type i = 0u; i < Rows; ++i)
        {
            for (size_type j = 0u; j < Columns; ++j)
            {
                err_matrix(i, j) = std::numeric_limits<T>::quiet_NaN();
            }
        }

        return err_matrix;
    }


    friend constexpr auto operator+(const std_matrix& lhs, const std_matrix& rhs) noexcept
    {
        std_matrix temp(lhs);
        return (temp += rhs);
    }


    friend constexpr auto operator-(const std_matrix& lhs, const std_matrix& rhs) noexcept
    {
        std_matrix temp(lhs);
        return (temp -= rhs);
    }


    template <std::size_t _Rows_lhs, std::size_t _Columns_lhs,
              std::size_t _Rows_rhs, std::size_t _Columns_rhs>
    friend constexpr auto operator*(const matrix<value_type, _Rows_lhs, _Columns_lhs>& lhs,
                                    const matrix<value_type, _Rows_rhs, _Columns_rhs>& rhs) noexcept
    {
        static_assert(_Columns_lhs == _Rows_rhs, "Incorrect matrix for product!");
    
        matrix<value_type, _Rows_lhs, _Columns_rhs> result{};
        container<value_type, _Rows_rhs> thatColumn{};

        for (size_type j = 0u; j < _Columns_rhs; ++j)
        {
        	for (size_type k = 0u; k < _Rows_rhs; ++k)
            {
        		thatColumn.at(k) = rhs(k, j);
            }

        	for (size_type i = 0u; i < _Rows_lhs; ++i)
            {
        		const auto thisRow = lhs(i);
        		value_type summand{};
        		for (size_type k = 0u; k < _Rows_rhs; ++k)
                {
        			summand += thisRow.at(k) * thatColumn.at(k);
                }
        		result(i, j) = summand;
        	}
        }
        return result;
    }


    friend constexpr auto operator*(const std_matrix& mat, const value_type num) noexcept
    {
        matrix temp(mat);
        return (temp *= num);
    }


    friend constexpr auto operator*(const value_type num, const std_matrix& mat) noexcept
    {
        return (mat * num);
    }


    friend constexpr auto operator/(const std_matrix& mat, const value_type num)
    {
        if (num == static_cast<value_type>(0))
        {
            throw std::invalid_argument("Division by 0!");
        }

        std_matrix temp(mat);
        return (temp /= num);
    }


    friend std::ostream& operator<<(std::ostream& os, const std_matrix& mat)
    {
        os << "[" << mat.get_dimension() << "]\n";
        for (const auto& row : mat._data)
        {
            std::copy(std::begin(row), std::end(row),
                      std::ostream_iterator<value_type>(os, " "));
            os << '\n';
        }
        return os;
    }


    friend std::istream& operator>>(std::istream& is, std_matrix& mat)
    {
        for (auto& row : mat._data)
        {
            for (auto& elem : row)
            {
                is >> elem;
            }
        }
        return is;
    }


private:
    container<container<value_type, _Columns>, _Rows> _data;

    constexpr auto _exp_helper(const std_matrix& m, const int num) const noexcept
    {
        if (num == 0)
        { 
            return create_identity<value_type, _Rows>();
        }
        else if (num == 1)
        {
            return m;
        }
        else if (num % 2 == 0)
        {
            // num is even.
            return _exp_helper(m * m, num/2);
        }
        else
        {
            // num is odd.
            return m * _exp_helper(m * m, (num - 1) / 2);
        }
    }
};

} // namespace vv
