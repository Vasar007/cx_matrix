// Copyright (C) 2018 Vasily Vasilyev (vasar007@yandex.ru)

#pragma once

#include <cstdint>
#include <limits>


/*
 *
 * Special thanks to Jason Turner for this constexpr random!
 * Origin video: https://www.youtube.com/watch?v=rpn_5Mrrxf8
 * Origin source code: https://godbolt.org/g/zbWvXK
 *
 */

namespace vv::cx_random
{

namespace
{

// Additional constant because of compiler has limited number of constexpr steps.
constexpr int MAX_DEPTH = 1000;

}

constexpr auto seed() noexcept
{
    std::uint64_t shifted = 0;

    for (const auto c : __TIME__)
    {
        shifted <<= 8;
        shifted |= c;
    }

    return shifted;
}


template <class T = std::uint64_t>
struct PCG
{
    using result_type = T;
    
    struct pcg32_random_t { std::uint64_t state=0;  std::uint64_t inc=seed(); };

    pcg32_random_t rng;

    
    constexpr result_type operator()() noexcept
    {
        return pcg32_random_r();
    }

    static result_type constexpr min() noexcept
    {
        return std::numeric_limits<result_type>::min();
    }

    static result_type constexpr max() noexcept
    {
        return std::numeric_limits<result_type>::max();
    }

private:
    constexpr std::uint64_t pcg32_random_r() noexcept
    {
        std::uint64_t oldstate = rng.state;
        // Advance internal state.
        rng.state = oldstate * 6364136225ULL + (rng.inc | 1); // 384679300
        // Calculate output function (XSH RR), uses old state for max ILP.
        std::uint64_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
        std::uint64_t rot = oldstate >> 59u;
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

};


constexpr auto get_random(int count) noexcept
{
    PCG<> pcg{};

    if (count > MAX_DEPTH)
    {
        count = count % MAX_DEPTH;
    }

    while (count > 0)
    {
        pcg();
        --count;
    }

    return pcg();
}

} // vv::namespace cx_random

/*
 * Example oof usage:
 * int main()
 * {
 *     constexpr auto r = get_random(10);
 *     return r;
 * }
 */
