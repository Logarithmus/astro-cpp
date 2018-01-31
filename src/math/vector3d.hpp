#ifndef ASTROCPP_MATH_VECTOR3D_HPP
#define ASTROCPP_MATH_VECTOR3D_HPP

#include <cmath>
#include <cstdio>

namespace Math
{
    class Vector3D
    {
        public:
            Vector3D(const double x, const double y, const double z):
                x(x),
                y(y),
                z(z)
            {

            }

            double mag() const
            {
                return std::sqrt(x*x + y*y + z*z);
            }
            void print() const
            {
                std::printf("\nx = %lf\n"
                            "y = %lf\n"
                            "z = %lf\n"
                            "magnitude = %lf\n", x, y, z, mag());
            }
            Vector3D& operator+() //prefix +
            {
                return *this;
            }
            Vector3D operator-() const  //prefix -
            {
                return Vector3D(-x, -y, -z);
            }
            Vector3D operator*(const double& c) const
            {
                return Vector3D(c*x, c*y, c*z);
            }
            Vector3D& operator*=(const double& c)
            {
                x *= c;
                y *= c;
                z *= c;
                return *this;
            }
            Vector3D operator+(const Vector3D& v2) const
            {
                return Vector3D(x+v2.x, y+v2.y, z+v2.z);
            }
            Vector3D operator-(const Vector3D& v2) const
            {
                return Vector3D(x-v2.x, y-v2.y, z-v2.z);
            }
            Vector3D& operator+=(const Vector3D& v2)
            {
                x += v2.x;
                y += v2.y;
                z += v2.z;
                return *this;
            }
            Vector3D& operator-=(const Vector3D& v2)
            {
                x -= v2.x;
                y -= v2.y;
                z -= v2.z;
                return *this;
            }
            static double dot(const Vector3D& v1, const Vector3D& v2)
            {
                return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
            }
            static Vector3D cross(const Vector3D& v1,
                                          const Vector3D& v2)
            {
                return Vector3D(v1.y*v2.z - v2.y*v1.z,
                                v1.x*v2.z - v2.x*v1.z,
                                v1.x*v2.y - v2.x*v1.y);
            }

        private:
            double x, y, z;
    };
} //namespace math

#endif // ASTROCPP_MATH_VECTOR3D_HPP
