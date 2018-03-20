#include <iostream>

#include "utils.hpp"
#include "cx_math.h"

// Need to include control unit first.
#include "control_unit.hpp"
#include "lu_decomposition.hpp"
#include "lup_decomposition.hpp"
#include "lupq_decomposition.hpp"
#include "qr_decomposition.hpp"
#include "iterations_methods.hpp"


// clnag parameter compilation: clang++ -std=c++1z -O2 -Xclang -flto-visibility-public-std -fconstexpr-steps=1000000000 -o main.exe main.cpp
int main()
{
    utils::pause("Press Enter to initiate ELIMINATION");

    std::cout.precision(15u);
    std::cout << std::fixed;

    // Matrix A and matrix b defined in control unit.
    std::cout << "Matrix A:\n" << vv::mat_A << "\n\n";    
    std::cout << "Matrix b:\n" << vv::vec_b;
    std::cout << "\n\n------------------------------------\n\n";

    vv::start_tests();

    std::cout << "\n\n";
    utils::pause();
    return 0;
}
