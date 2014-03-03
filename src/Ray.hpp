#ifndef __Ray_h__
#define __Ray_h__

#include "Vector3.hpp"

//==============================================================================
// Template Ray
//==============================================================================
template<class Real>
struct Ray
{
    Ray();
    Ray(const Vector3<Real>& origin, const Vector3<Real>& direction);

    Vector3<Real> o;
    Vector3<Real> d;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Ray<float>  Rayf;
typedef Ray<double> Rayd;

//------------------------------------------------------------------------------
// Ray::Ray
//------------------------------------------------------------------------------
template<typename Real>
inline
Ray<Real>::Ray()
{
}

//------------------------------------------------------------------------------
// Ray::Ray
//------------------------------------------------------------------------------
template<typename Real>
inline
Ray<Real>::Ray(const Vector3<Real>& origin, const Vector3<Real>& direction)
: o(origin), d(direction)
{
}

#endif
