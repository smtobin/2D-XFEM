#include "common.hpp"

#include <iostream>

int main()
{
    Vec3r test(0.1, 0.2, 0.3);

    std::cout << test.transpose() << std::endl;
}