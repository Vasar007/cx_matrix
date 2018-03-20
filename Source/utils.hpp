#pragma once

#include <iostream>
#include <iterator>
#include <numeric>
#include <string_view>


namespace utils
{

template <class Container>
void print(const Container& container)
{
    std::copy(std::begin(container), std::end(container),
              std::ostream_iterator<typename Container::value_type>(std::cout, " "));
}

void pause(const std::string_view message = "\nPress the Enter key to continue...")
{
    do
    {
        std::cout << message;
    }
    while (std::cin.get() != '\n');
}

void pause_clear(const std::string_view message = "Press ENTER to continue...")
{
    std::cout << message << std::flush;
    std::cin.seekg(0u, std::ios::end);
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

} // namespace utils
