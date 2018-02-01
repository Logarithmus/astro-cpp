#ifndef ASTROCPP_MATH_VECTOR3D_HPP
#define ASTROCPP_MATH_VECTOR3D_HPP

#include <cmath>
#include <cstdio>

namespace Math
{
    class Vector3D
    {
        public:
            double x, y, z;

            Vector3D(const double& x, const double& y, const double& z);

            double mag() const;
            void print() const;

            Vector3D& operator+(); //prefix +
            Vector3D operator-() const;  //prefix -

            Vector3D operator+(const Vector3D& v2) const;
            Vector3D operator-(const Vector3D& v2) const;

            Vector3D& operator*=(const double& c);
            Vector3D& operator+=(const Vector3D& v2);
            Vector3D& operator-=(const Vector3D& v2);

            Vector3D normalized();
            static double dot_product(const Vector3D& v1, const Vector3D& v2);
            static Vector3D cross_product(const Vector3D& v1,
                                  const Vector3D& v2);
    };

    //Vector3D methods' implementation
    Vector3D::Vector3D(const double& x, const double& y, const double& z):
        x(x),
        y(y),
        z(z)
    {

    }

    inline double Vector3D::mag() const
    {
        return std::sqrt(x*x + y*y + z*z);
    }

    inline void Vector3D::print() const
    {
        std::printf("\nx = %lf\n"
                    "y = %lf\n"
                    "z = %lf\n"
                    "magnitude = %lf\n", x, y, z, mag());
    }

    inline Vector3D& Vector3D::operator+() //prefix +
    {
        return *this;
    }

    inline Vector3D Vector3D::operator-() const  //prefix -
    {
        return Vector3D(-x, -y, -z);
    }

    inline Vector3D& Vector3D::operator*=(const double& c)
    {
        x *= c;
        y *= c;
        z *= c;
        return *this;
    }

    inline Vector3D operator*(const double& c, Vector3D v)
    {
        return v *= c;
    }

    inline Vector3D Vector3D::normalized()
    {
        double c = 1 / mag();
        return c*(*this);
    }

    inline Vector3D Vector3D::operator+(const Vector3D& v2) const
    {
        return Vector3D(x+v2.x, y+v2.y, z+v2.z);
    }

    inline Vector3D Vector3D::operator-(const Vector3D& v2) const
    {
        return Vector3D(x-v2.x, y-v2.y, z-v2.z);
    }

    inline Vector3D& Vector3D::operator+=(const Vector3D& v2)
    {
        x += v2.x;
        y += v2.y;
        z += v2.z;
        return *this;
    }

    inline Vector3D& Vector3D::operator-=(const Vector3D& v2)
    {
        x -= v2.x;
        y -= v2.y;
        z -= v2.z;
        return *this;
    }

    inline double Vector3D::dot_product(const Vector3D& v1,
                                        const Vector3D& v2)
    {
        return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    }

    inline Vector3D Vector3D::cross_product(const Vector3D& v1,
                                            const Vector3D& v2)
    {
        return Vector3D(v1.y*v2.z - v2.y*v1.z,
                        v1.x*v2.z - v2.x*v1.z,
                        v1.x*v2.y - v2.x*v1.y);
    }
} //namespace math

#endif // ASTROCPP_MATH_VECTOR3D_HPP
