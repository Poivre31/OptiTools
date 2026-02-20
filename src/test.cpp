#include <iostream>
#include <math/vec3d.h>
#include <solver/euler.h>

int main()
{
    solver euler(2.79);
    euler.set_initial_conditions(3.1, 2042.2);
    euler.solve(3);
}