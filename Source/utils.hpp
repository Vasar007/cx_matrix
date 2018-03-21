#pragma once

#include <iostream>
#include <iterator>
#include <numeric>
#include <string_view>


namespace utils
{

template <class OutputStream, class Container>
void print(OutputStream& out, const Container& container) noexcept
{
    std::copy(std::begin(container), std::end(container),
              std::ostream_iterator<typename Container::value_type>(out, " "));
}

void pause(const std::string_view message = "\nPress the Enter key to continue...") noexcept
{
    do
    {
        std::cout << message;
    }
    while (std::cin.get() != '\n');
}

void pause_clear(const std::string_view message = "Press ENTER to continue...") noexcept
{
    std::cout << message << std::flush;
    std::cin.seekg(0u, std::ios::end);
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

} // namespace utils
