#include "math/vector3d.hpp"
#include "astro/orbit.hpp"

int main()
{
    Math::Vector3D e1 = {1, 0 ,0},
                   e2 = {0, 1, 0},
                   e3 = Math::Vector3D::cross_product(e1, e2),
                   e4 = Math::Vector3D::cross_product(e2, e1);
    e3.print();
    e4.print();
    return 0;
}
