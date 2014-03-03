#ifndef __Function_h__
#define __Function_h__

#include "Common.hpp"

//==============================================================================
// Template Bivariate
//==============================================================================
template<typename Real>
class Bivariate
{
public:
    typedef Real (*Function)(const Real, const Real);

    Bivariate() : _f(0) {}
    Bivariate(const Function f) : _f(f) {}
    void SetFunction(const Function f) { _f = f; }

    Real operator()(const Real x, const Real y) const
    {
        assert(_f);
        return (*_f)(x, y);
    }

private:
    Function _f;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Bivariate<float>    Bivariatef;
typedef Bivariate<double>   Bivariated;

//==============================================================================
// Template Trivariate
//==============================================================================
template<typename Real>
class Trivariate
{
public:
    typedef Real (*Function)(const Real, const Real, const Real);

    Trivariate() : _f(0) {}
    Trivariate(const Function f) : _f(f) {}
    void SetFunction(const Function f) { _f = f; }

    Real operator()(const Real x, const Real y, const Real z) const
    {
        assert(_f);
        return (*_f)(x, y, z);
    }

private:
    Function _f;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Trivariate<float>   Trivariatef;
typedef Trivariate<double>  Trivariated;

//==============================================================================
// Template Trivariate
//==============================================================================
template<typename Real>
class Quadrivariate
{
public:
    typedef Real (*Function)(const Real, const Real, const Real, const Real);

    Quadrivariate() : _f(0) {}
    Quadrivariate(const Function f) : _f(f) {}
    void SetFunction(const Function f) { _f = f; }

    Real operator()(const Real x,
                    const Real y,
                    const Real z,
                    const Real w) const
    {
        assert(_f);
        return (*_f)(x, y, z);
    }

private:
    Function _f;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Quadrivariate<float>    Quadrivariatef;
typedef Quadrivariate<double>   Quadrivariated;

#endif
