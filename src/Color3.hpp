#ifndef __Color3_h__
#define __Color3_h__

#include "Common.hpp"

/*
 Color3
*/
template<typename T>
struct Color3
{
    Color3();
    Color3(const Color3&);
    Color3(const T);
    Color3(const T, const T, const T);

    T& operator[](const uint32_t);
    const T& operator[](const uint32_t) const;

    void operator+=(const Color3&);
    void operator-=(const Color3&);
    void operator*=(const Color3&);
    void operator/=(const Color3&);
    void operator*=(const T);
    void operator/=(const T);

    Color3 operator-() const;

    Color3 operator+(const Color3&) const;
    Color3 operator-(const Color3&) const;
    Color3 operator*(const Color3&) const;
    Color3 operator/(const Color3&) const;
    Color3 operator*(const T) const;
    Color3 operator/(const T) const;

    void Clamp(const T low = T(0), const T high = T(1));

    T* Ptr();
    const T* Ptr() const;

    T r, g, b;
};

typedef Color3<uint8_t>     Color3ub;
typedef Color3<uint8_t>     Color3uc;
typedef Color3<uint16_t>    Color3us;
typedef Color3<uint32_t>    Color3ui;
typedef Color3<float>       Color3f;
typedef Color3<double>      Color3d;

void Convert(Color3f& result, Color3ub& color);
void Convert(Color3ub& result, Color3f& color);

template<typename T>
Color3<T>::Color3()
: r(T(0)), g(T(0)), b(T(0))
{
}

template<typename T>
Color3<T>::Color3(const Color3& c)
: r(c.r), g(c.g), b(c.b)
{
}

template<typename T>
Color3<T>::Color3(const T s)
: r(s), g(s), b(s)
{
}

template<typename T>
Color3<T>::Color3(const T r0, const T g0, const T b0)
: r(r0), g(g0), b(b0)
{
}

template<typename T>
T& Color3<T>::operator[](const uint32_t i)
{
    assert(i < 3);
    return *(&r + i);
}

template<typename T>
const T& Color3<T>::operator[](const uint32_t i) const
{
    assert(i < 3);
    return *(&r + i);
}

template<typename T>
void Color3<T>::operator+=(const Color3& c)
{
    r += c.r;
    g += c.g;
    b += c.b;
}

template<typename T>
void Color3<T>::operator-=(const Color3& c)
{
    r -= c.r;
    g -= c.g;
    b -= c.b;
}

template<typename T>
void Color3<T>::operator*=(const Color3& c)
{
    r *= c.r;
    g *= c.g;
    b *= c.b;
}

template<typename T>
void Color3<T>::operator/=(const Color3& c)
{
    r /= c.r;
    g /= c.g;
    b /= c.b;
}

template<typename T>
void Color3<T>::operator*=(const T s)
{
    r *= s;
    g *= s;
    b *= s;
}

template<typename T>
void Color3<T>::operator/=(const T s)
{
    r /= s;
    g /= s;
    b /= s;
}

template<typename T>
Color3<T> Color3<T>::operator-() const
{
    return Color3(-r, -g, -b);
}

template<typename T>
Color3<T> Color3<T>::operator+(const Color3<T>& c) const
{
    return Color3(r + c.r, g + c.g, b + c.b);
}

template<typename T>
Color3<T> Color3<T>::operator-(const Color3<T>& c) const
{
    return Color3(r - c.r, g - c.g, b - c.b);
}

template<typename T>
Color3<T> Color3<T>::operator*(const Color3<T>& c) const
{
    return Color3(r * c.r, g * c.g, b * c.b);
}

template<typename T>
Color3<T> Color3<T>::operator/(const Color3<T>& c) const
{
    return Color3(r / c.r, g / c.g, b / c.b);
}

template<typename T>
Color3<T> Color3<T>::operator*(const T s) const
{
    return Color3(s * r, s * g, s * b);
}

template<typename T>
Color3<T> Color3<T>::operator/(const T s) const
{
    T sInv = T(1) / s;
    return Color3(r * sInv, g * sInv, b * sInv);
}

template<typename T>
inline Color3<T>
operator*(const T s, const Color3<T>& c)
{
    return Color3<T>(s * c.r, s * c.g, s * c.b);
}

template<typename T>
void Color3<T>::Clamp(const T low, const T high)
{
    if (r < low) {
        r = low;
    }
    if (g < low) {
        g = low;
    }
    if (b < low) {
        b = low;
    }
    if (r > high) {
        r = high;
    }
    if (g > high) {
        g = high;
    }
    if (b > high) {
        b = high;
    }
}

template<typename T>
T* Color3<T>::Ptr()
{
    return &r;
}

template<typename T>
const T* Color3<T>::Ptr() const
{
    return &r;
}

inline void Convert(Color3f& result, Color3ub& color)
{
    result.r = float(color.r) / 255.0f;
    result.g = float(color.g) / 255.0f;
    result.b = float(color.b) / 255.0f;
}

inline void Convert(Color3ub& result, Color3f& color)
{
    result.r = uint8_t(color.r * 255.0f);
    result.g = uint8_t(color.g * 255.0f);
    result.b = uint8_t(color.b * 255.0f);
}

#endif
