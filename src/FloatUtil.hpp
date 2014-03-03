//==============================================================================
// Author: Tim Thirion
//
// This file is: FloatUtil.h
//   Created on: July 18, 2004
//  Modified on: July 18, 2004
//
//  Description:
//				Floating point constants and utility functions are defined in
//				this file. It should be noted that the class FloatUtil should
//				only use floating-point inherent C++ types as its sole template
//				argument. Some routines would work with integers, but don't
//				count on it.
//
//    Revisions:
//				1. July 17, 2004
//					FloatUtil.h created and added to the system. Constants
//					PI, DEG_TO_RAD, RAD_TO_DEG, EPSILON, and INFINITY added.
//					The methods Abs, IsZero, Compare, ToDegrees, ToRadians,
//					Sqrt, and RSqrt also added. The implementation for
//					FloatUtil<float>::RSqrt is a hack taken from the Quake III
//					source code. It is efficient although the code looks quite
//					"magical." There has been some study of it; more information
//					to follow.
//
// Added Min, Max, Pow.
// Changed Lerp and changed IsZero to "NearlyZero".
//
//==============================================================================
#ifndef __FloatUtil_h__
#define __FloatUtil_h__

#include <cmath>

//==============================================================================
// class template FloatUtil
//==============================================================================
template<typename Real>
struct FloatUtil
{
	//==========================================================================
	// Public Constants
	//==========================================================================
	static const Real PI;
	static const Real DEG_TO_RAD;
	static const Real RAD_TO_DEG;
	static const Real EPSILON;
	static const Real INFINITY;

	//==========================================================================
	// Public Methods
	//==========================================================================
	static Real Abs(const Real);
	static Real Min(const Real, const Real);
	static Real Max(const Real, const Real);
	static bool NearlyZero(const Real);
	static bool NearlyEqual(const Real, const Real);
	static Real ToDegrees(const Real);
	static Real ToRadians(const Real);
	static Real Pow(const Real, const Real);
	static Real Sqrt(const Real);
	static Real InverseSqrt(Real);
	static void Lerp(Real&, const Real, const Real, const Real);

	static Real Exp(const Real);
	static Real Ln(const Real);

	static Real Sin(const Real);
	static Real Cos(const Real);
	static Real Tan(const Real);
	static Real ArcSin(const Real);
	static Real ArcCos(const Real);
	static Real ArcTan(const Real);
	static Real ArcTan2(const Real, const Real);
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef FloatUtil<float>	FloatUtilf;
typedef FloatUtil<double>	FloatUtild;

//==============================================================================
// Constants
//==============================================================================
const float FloatUtil<float>::PI			= 3.141592653589793238462643383279f;
const float FloatUtil<float>::DEG_TO_RAD	= 0.017453292519943295769236907685f;
const float FloatUtil<float>::RAD_TO_DEG	= 57.29577951308232087679815481411f;
const float FloatUtil<float>::EPSILON		= 1e-3f;
const float FloatUtil<float>::INFINITY		= 1e10f;

const double FloatUtil<double>::PI			= 3.141592653589793238462643383279;
const double FloatUtil<double>::DEG_TO_RAD	= 0.017453292519943295769236907685;
const double FloatUtil<double>::RAD_TO_DEG	= 57.29577951308232087679815481411;
const double FloatUtil<double>::EPSILON		= 1e-3;
const double FloatUtil<double>::INFINITY	= 1e10;

//------------------------------------------------------------------------------
// FloatUtil::Abs
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::Abs(const Real x)
{
	return x < Real(0) ? -x : x;
}

//------------------------------------------------------------------------------
// FloatUtil::Min
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::Min(const Real a, const Real b)
{
	return a < b ? a : b;
}

//------------------------------------------------------------------------------
// FloatUtil::Max
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::Max(const Real a, const Real b)
{
	return a < b ? b : a;
}

//------------------------------------------------------------------------------
// FloatUtil::NearlyZero
//------------------------------------------------------------------------------
template<typename Real>
inline bool
FloatUtil<Real>::NearlyZero(const Real x)
{
	return Abs(x) < EPSILON;
}

//------------------------------------------------------------------------------
// FloatUtil::NearlyEqual
//------------------------------------------------------------------------------
template<typename Real>
inline bool
FloatUtil<Real>::NearlyEqual(const Real x, const Real y)
{
	return Abs(x - y) < EPSILON;
}

//------------------------------------------------------------------------------
// FloatUtil::ToDegrees
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::ToDegrees(const Real radians)
{
	return RAD_TO_DEG * radians;
}

//------------------------------------------------------------------------------
// FloatUtil::ToRadians
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::ToRadians(const Real degrees)
{
	return DEG_TO_RAD * degrees;
}

//------------------------------------------------------------------------------
// FloatUtil::Pow
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::Pow(const float x, const float y)
{
	return ::powf(x, y);
}

//------------------------------------------------------------------------------
// FloatUtil::Pow
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::Pow(const double x, const double y)
{
	return ::pow(x, y);
}

//------------------------------------------------------------------------------
// FloatUtil::Sqrt
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::Sqrt(const float x)
{
	return ::sqrtf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Sqrt
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::Sqrt(const double x)
{
	return ::sqrt(x);
}

//------------------------------------------------------------------------------
// FloatUtil::InverseSqrt
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::InverseSqrt(float x)
{
	float halfX = 0.5f * x;
	int i = *(int*)&x;
	i = 0x5f3759df - (i >> 1);
	x = *(float*)&i;
	x = x * (1.5f - halfX * x * x);
	return x;
}

//------------------------------------------------------------------------------
// FloatUtil::InverseSqrt
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::InverseSqrt(double x)
{
	return 1.0 / ::sqrt(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Lerp
//------------------------------------------------------------------------------
template<typename Real>
inline void
FloatUtil<Real>::Lerp(Real& r, const Real a, const Real b, const Real s)
{
	r = a + s * (b - a);
}

//------------------------------------------------------------------------------
// FloatUtil::Exp
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::Exp(const float x)
{
	return ::expf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Exp
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::Exp(const double x)
{
	return ::exp(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Ln
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::Ln(const float x)
{
	return ::logf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Ln
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::Ln(const double x)
{
	return ::log(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Sin
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::Sin(const float radians)
{
	return ::sinf(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::Sin
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::Sin(const double radians)
{
	return ::sin(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::Cos
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::Cos(const float radians)
{
	return ::cosf(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::Cos
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::Cos(const double radians)
{
	return ::cos(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::Tan
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::Tan(const float radians)
{
	return ::tanf(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::Tan
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::Tan(const double radians)
{
	return ::tan(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcSin
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::ArcSin(const float x)
{
	return ::asinf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcSin
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::ArcSin(const double x)
{
	return ::asin(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcCos
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::ArcCos(const float x)
{
	return ::acosf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcCos
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::ArcCos(const double x)
{
	return ::acos(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcTan
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::ArcTan(const float x)
{
	return ::atanf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcTan
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::ArcTan(const double x)
{
	return ::atan(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcTan2
//------------------------------------------------------------------------------
inline float
FloatUtil<float>::ArcTan2(const float x, const float y)
{
	return ::atan2f(x, y);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcTan2
//------------------------------------------------------------------------------
inline double
FloatUtil<double>::ArcTan2(const double x, const double y)
{
	return ::atan2(x, y);
}

#endif