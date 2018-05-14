# CxMatrix

Constexpr matrices with computational mathematic methods. (#constexpr_ALL_the_things)
Inspired by https://www.youtube.com/watch?v=HMB9oXFobJc

## Screenshots

| My proof                                                         |
|------------------------------------------------------------------|
| ![Screenshot](Media/proof.jpg "Proof")                           |


## Compiling

This project is compiled by Clang v5.0.1 and parameters: *clang++ -std=c++1z -O2 -Wall -Wextra -pedantic -Xclang -flto-visibility-public-std -fconstexpr-steps=1000000000*

Also it would be compiled by GCC v7.2 and higher.

MSVC from the v19.14 (VS v15.7) has some issues with constexpr statements and sometimes "Internal compiler error" occurs. But code could be compiled too.


## License information

This project is licensed under the terms of the [MIT License](LICENSE).
