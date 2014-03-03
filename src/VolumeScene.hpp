#ifndef __VolumeScene_h__
#define __VolumeScene_h__

#include <string>
#include <vector>
#include <GLUT/glut.h>
#include <fstream>

#include "Common.hpp"
#include "Function.hpp"
#include "Vector2.hpp"
#include "Vector3.hpp"
#include "Box3.hpp"
#include "Color3.hpp"

//==============================================================================
// Class VolumeScene
//==============================================================================
class VolumeScene
{
public:
    //==========================================================================
    // Public Methods
    //==========================================================================

    //--------------------------------------------------------------------------
    // VolumeScene
    //--------------------------------------------------------------------------
    VolumeScene();

    //--------------------------------------------------------------------------
    // ~MarchingCubes
    //--------------------------------------------------------------------------
    ~VolumeScene();

    //--------------------------------------------------------------------------
    // SetDimensions
    //--------------------------------------------------------------------------
    void SetDimensions(const Vector3ui& dimensions);

    //--------------------------------------------------------------------------
    // SetImageResolution
    //--------------------------------------------------------------------------
    void SetImageResolution(const uint32_t width, const uint32_t height);

    //--------------------------------------------------------------------------
    // SetImageWorldBounds
    //--------------------------------------------------------------------------
    void SetImageWorldBounds(const Vector2f& mins, const Vector2f& maxes);

    //--------------------------------------------------------------------------
    // SetImageBasis
    //--------------------------------------------------------------------------
    void SetImageBasis(const Vector3f& x, const Vector3f& y);

    //--------------------------------------------------------------------------
    // SetEyePosition
    //--------------------------------------------------------------------------
    void SetEyePosition(const Vector3f& eye);

    //--------------------------------------------------------------------------
    // LookAt
    //--------------------------------------------------------------------------
    void LookAt(const Vector3f& at);

    //--------------------------------------------------------------------------
    // SetNumResamples
    //--------------------------------------------------------------------------
    void SetNumResamples(const uint32_t resamples);

    //--------------------------------------------------------------------------
    // SetRayStepSize
    //--------------------------------------------------------------------------
    void SetRayStepSize(const float stepSize);

    //--------------------------------------------------------------------------
    // DataFromFile
    //--------------------------------------------------------------------------
    bool DataFromMRIFile(const std::string& filename, const bool isUnix = true);

    //--------------------------------------------------------------------------
    // DataFromFunction
    //--------------------------------------------------------------------------
    void DataFromFunction(const Trivariatef& t);

    //--------------------------------------------------------------------------
    // SetIsovalueOne
    //--------------------------------------------------------------------------
    void SetIsovalues(const Vector2f& isovalues);

    //--------------------------------------------------------------------------
    // SetThicknesses
    //--------------------------------------------------------------------------
    void SetThicknesses(const Vector2f& thicknesses);

    //--------------------------------------------------------------------------
    // SetPrefactors
    //--------------------------------------------------------------------------
    void SetPrefactors(const Vector2f& prefactors);

    //--------------------------------------------------------------------------
    // SetGradientPower
    //--------------------------------------------------------------------------
    void SetGradientPowers(const Vector2f& powers);

    //--------------------------------------------------------------------------
    // SetAmbientLight
    //--------------------------------------------------------------------------
    void SetAmbientLight(const Color3f& color);

    //--------------------------------------------------------------------------
    // SetAmbientContribution
    //--------------------------------------------------------------------------
    void SetAmbientContribution(const float Ka);

    //--------------------------------------------------------------------------
    // SetDiffuseContribution
    //--------------------------------------------------------------------------
    void SetDiffuseContribution(const float Kd);

    //--------------------------------------------------------------------------
    // SetSpecularContribution
    //--------------------------------------------------------------------------
    void SetSpecularContribution(const float Ks);

    //--------------------------------------------------------------------------
    // SetPhongPower
    //--------------------------------------------------------------------------
    void SetPhongPower(const float power);

    //--------------------------------------------------------------------------
    // SetPointLight
    //--------------------------------------------------------------------------
    void SetPointLight(const Vector3f& position, const Color3f& color);

    //--------------------------------------------------------------------------
    // SetConstantAttenuation
    //--------------------------------------------------------------------------
    void SetConstantAttenuation(const float constant);

    //--------------------------------------------------------------------------
    // SetLinearAttenuation
    //--------------------------------------------------------------------------
    void SetLinearAttenuation(const float linear);

    //--------------------------------------------------------------------------
    // SetFirstObjectColor
    //--------------------------------------------------------------------------
    void SetFirstObjectColor(const Color3f& color);

    //--------------------------------------------------------------------------
    // SetSecondObjectMaterial
    //--------------------------------------------------------------------------
    void SetSecondObjectColor(const Color3f& color);

    //--------------------------------------------------------------------------
    // SetThresholdIsovalue
    //--------------------------------------------------------------------------
    void SetThresholdIsovalue(const float threshold);

    //--------------------------------------------------------------------------
    // HasErrors
    //--------------------------------------------------------------------------
    bool HasErrors() const;

    //--------------------------------------------------------------------------
    // GetErrorString
    //--------------------------------------------------------------------------
    std::string GetErrorString() const;

    //--------------------------------------------------------------------------
    // BuildImage
    //--------------------------------------------------------------------------
    void BuildImage();

    //--------------------------------------------------------------------------
    // Render
    //--------------------------------------------------------------------------
    void Render();

private:
    //==========================================================================
    // Private Methods
    //==========================================================================

    //--------------------------------------------------------------------------
    // Copy Constructor
    //
    // The copy constructor is not implemented. Prevent its use.
    //--------------------------------------------------------------------------
    VolumeScene(const VolumeScene&);

    //--------------------------------------------------------------------------
    // Assignment Operator
    //
    // The assignment operator is not implemented. Prevent its use.
    //--------------------------------------------------------------------------
    const VolumeScene& operator=(const VolumeScene&);

    //--------------------------------------------------------------------------
    // Index
    //
    // Get the index of the one-dimensional array from the three-dimensional
    // coordinates.
    //--------------------------------------------------------------------------
    uint32_t Index(const uint32_t i, const uint32_t j, const uint32_t k) const;

    //--------------------------------------------------------------------------
    // Sample
    //
    // Get the density value at the voxel grid position (i, j, k).
    //--------------------------------------------------------------------------
    float Sample(const uint32_t i, const uint32_t j, const uint32_t k) const;

    //--------------------------------------------------------------------------
    // SampleTrilinear
    //
    // Use trilinear interpolation to compute the density value at the position
    // (x, y, z). This method returns false if the point is outside the volume,
    // returns true otherwise.
    //--------------------------------------------------------------------------
    bool SampleTrilinear(float& value,
                         const float x,
                         const float y,
                         const float z) const;

    //--------------------------------------------------------------------------
    // SampleTrilinear
    //
    // Use trilinear interpolation to compute the density value at the position
    // p = (x, y, z).
    //--------------------------------------------------------------------------
    bool SampleTrilinear(float& value, const Vector3f& p) const;

    //--------------------------------------------------------------------------
    // GetGradient
    //
    // Compute the gradient at voxel grid position (i, j, k).
    //--------------------------------------------------------------------------
    void GetGradient(Vector3f& g,
                     const uint32_t i,
                     const uint32_t j,
                     const uint32_t k) const;

    //--------------------------------------------------------------------------
    // GetGradientTrilinear
    //
    //  Use trilinear interpolation to compute the gradient at the position
    //  (x, y, z).
    //--------------------------------------------------------------------------
    bool GetGradientTrilinear(Vector3f& g,
                              const float x,
                              const float y,
                              const float z) const;

    //--------------------------------------------------------------------------
    // GetGradientTrilinear
    //
    // Use trilinear interpolation to compute the gradient at the position
    // p = (x, y, z).
    //--------------------------------------------------------------------------
    bool GetGradientTrilinear(Vector3f& g, const Vector3f& p) const;

    //--------------------------------------------------------------------------
    // Alpha
    //
    // The opacity transfer function.
    //--------------------------------------------------------------------------
    float Alpha(const Vector3f& r);

    //--------------------------------------------------------------------------
    // ComputeColorByIllumination
    //--------------------------------------------------------------------------
    void ComputeColorByIllumination(Color3f& color,
                                    const Vector3f& position,
                                    const Vector3f& rayOrigin);

    //==========================================================================
    // Private Data
    //==========================================================================

    //
    // Voxel grid information
    //
    Vector3ui                   _dimensions;
    uint32_t                    _dimensionsProduct;
    Vector3f                    _scale;
    Vector3f                    _widths;
    Vector3f                    _rotationAngles;
    Box3f                       _boundingBox;

    //
    // Voxel data
    //
    bool                        _renderingMRI;
    float                       _headData[256][256][109];
    float*                      _data;

    //
    // Image information
    //
    Vector2ui                   _imageResolution;
    Vector2f                    _imageMins;
    Vector2f                    _imageMaxes;
    Vector3f                    _x;
    Vector3f                    _y;
    Color3ub*                   _pixels;

    //
    // Raycasting information
    //
    uint32_t                    _numResamples;
    float                       _rayStepSize;

    //
    // Camera information
    //
    Vector3f                    _eyePosition;
    Vector3f                    _eyeDirection;

    //
    // Lighting information
    //
    Vector3f                    _lightPosition;
    Color3f                     _lightColor;
    float                       _Ka;
    float                       _Kd;
    float                       _Ks;
    float                       _linearAtten;
    float                       _constantAtten;
    float                       _phongPower;
    Color3f                     _ambientLight;
    Color3f                     _firstColor;
    Color3f                     _secondColor;

    //
    // Visualization parameters
    //
    Vector2f                    _isovalues;
    float                       _thresholdIsovalue;
    Vector2f                    _thicknesses;
    Vector2f                    _prefactors;
    Vector2f                    _gradientPowers;

    //
    // Error information
    //
    bool                        _hasErrors;
    std::string                 _errorString;
};

//------------------------------------------------------------------------------
// VolumeScene::VolumeScene
//------------------------------------------------------------------------------
inline
VolumeScene::VolumeScene()
:
_dimensionsProduct(0),
//_data(0),
_pixels(0),
_numResamples(0),
_rayStepSize(0),
_hasErrors(false)
{
}

//------------------------------------------------------------------------------
// VolumeScene::~VolumeScene
//------------------------------------------------------------------------------
inline
VolumeScene::~VolumeScene()
{
    //if (_data) {
    //  delete [] _data;
    //}
    //_data = 0;
    if (_pixels) {
        delete [] _pixels;
    }
    _pixels = 0;
}

//------------------------------------------------------------------------------
// VolumeScene::SetDimensions
//------------------------------------------------------------------------------
inline void
VolumeScene::SetDimensions(const Vector3ui& dimensions)
{
    _dimensionsProduct = dimensions.x * dimensions.y * dimensions.z;
    _dimensions = dimensions;

    Vector3f maxes = Vector3f(dimensions.x - 1,
                              dimensions.y - 1,
                              dimensions.z - 1);
    //
    // The bounding box is an axis-aligned bounding box with minimums
    // (0, 0, 0) and maximums (dimensions[0], dimensions[1], dimensions[2]).
    //
    _boundingBox.SetCenter(maxes * 0.5f);

    _boundingBox.SetAxis(0, Vector3f(1.0f, 0.0f, 0.0f));
    _boundingBox.SetAxis(1, Vector3f(0.0f, 1.0f, 0.0f));
    _boundingBox.SetAxis(2, Vector3f(0.0f, 0.0f, 1.0f));
    _boundingBox.SetExtent(0, maxes.x);
    _boundingBox.SetExtent(1, maxes.y);
    _boundingBox.SetExtent(2, maxes.z);

    _widths = Vector3f(1.0f, 1.0f, 1.0f);
}

//------------------------------------------------------------------------------
// VolumeScene::SetImageResolution
//------------------------------------------------------------------------------
inline void
VolumeScene::SetImageResolution(const uint32_t width, const uint32_t height)
{
    if (_pixels) {
        delete [] _pixels;
    }
    _pixels = new Color3ub[width * height];
    assert(_pixels);
    _imageResolution = Vector2ui(width, height);
}

//------------------------------------------------------------------------------
// VolumeScene::SetImageBasis
//------------------------------------------------------------------------------
inline void
VolumeScene::SetImageBasis(const Vector3f& x, const Vector3f& y)
{
    _x = x;
    Normalize(_x);
    _y = y;
    Normalize(_y);
}

//------------------------------------------------------------------------------
// VolumeScene::SetImageBounds
//------------------------------------------------------------------------------
inline void
VolumeScene::SetImageWorldBounds(const Vector2f& mins, const Vector2f& maxes)
{
    _imageMins = mins;
    _imageMaxes = maxes;
}

//------------------------------------------------------------------------------
// VolumeScene::SetEyePosition
//------------------------------------------------------------------------------
inline void
VolumeScene::SetEyePosition(const Vector3f& position)
{
    _eyePosition = position;
}

//------------------------------------------------------------------------------
// VolumeScene::LookAt
//------------------------------------------------------------------------------
inline void
VolumeScene::LookAt(const Vector3f& direction)
{
    _eyeDirection = direction;
    Normalize(_eyeDirection);
}

//------------------------------------------------------------------------------
// VolumeScene::SetNumResamples
//------------------------------------------------------------------------------
inline void
VolumeScene::SetNumResamples(const uint32_t numResamples)
{
    _numResamples = numResamples;
}

//------------------------------------------------------------------------------
// VolumeScene::SetRayStepSize
//------------------------------------------------------------------------------
inline void
VolumeScene::SetRayStepSize(const float stepSize)
{
    _rayStepSize = stepSize;
}

//------------------------------------------------------------------------------
// VolumeScene::SetIsovalueOne
//------------------------------------------------------------------------------
inline void
VolumeScene::SetIsovalues(const Vector2f& isovalues)
{
    _isovalues = isovalues;
}

//------------------------------------------------------------------------------
// VolumeScene::SetThicknesses
//------------------------------------------------------------------------------
inline void
VolumeScene::SetThicknesses(const Vector2f& thicknesses)
{
    _thicknesses = thicknesses;
}

//------------------------------------------------------------------------------
// VolumeScene::SetPrefactors
//------------------------------------------------------------------------------
inline void
VolumeScene::SetPrefactors(const Vector2f& prefactors)
{
    _prefactors = prefactors;
}

//------------------------------------------------------------------------------
// VolumeScene::SetGradientPower
//------------------------------------------------------------------------------
inline void
VolumeScene::SetGradientPowers(const Vector2f& powers)
{
    _gradientPowers = powers;
}

//------------------------------------------------------------------------------
// VolumeScene::SetAmbientLight
//------------------------------------------------------------------------------
inline void
VolumeScene::SetAmbientLight(const Color3f& color)
{
    _ambientLight = color;
}

//------------------------------------------------------------------------------
// VolumeScene::SetAmbientContribution
//------------------------------------------------------------------------------
inline void
VolumeScene::SetAmbientContribution(const float Ka)
{
    _Ka = Ka;
}

//------------------------------------------------------------------------------
// VolumeScene::SetDiffuseContribution
//------------------------------------------------------------------------------
inline void
VolumeScene::SetDiffuseContribution(const float Kd)
{
    _Kd = Kd;
}

//------------------------------------------------------------------------------
// VolumeScene::SetSpecularContribution
//------------------------------------------------------------------------------
inline void
VolumeScene::SetSpecularContribution(const float Ks)
{
    _Ks = Ks;
}

//------------------------------------------------------------------------------
// SetPhongPower
//------------------------------------------------------------------------------
inline void
VolumeScene::SetPhongPower(const float power)
{
    _phongPower = power;
}

//------------------------------------------------------------------------------
// VolumeScene::SetPointLight
//------------------------------------------------------------------------------
inline void
VolumeScene::SetPointLight(const Vector3f& position, const Color3f& color)
{
    _lightPosition = position;
    _lightColor = color;
}

//------------------------------------------------------------------------------
// VolumeScene::SetConstantAttenuation
//------------------------------------------------------------------------------
inline void
VolumeScene::SetConstantAttenuation(const float constant)
{
    _constantAtten = constant;
}

//------------------------------------------------------------------------------
// VolumeScene::SetLinearAttenuation
//------------------------------------------------------------------------------
inline void
VolumeScene::SetLinearAttenuation(const float linear)
{
    _linearAtten = linear;
}

//------------------------------------------------------------------------------
// VolumeScene::SetFirstObjectColor
//------------------------------------------------------------------------------
inline void
VolumeScene::SetFirstObjectColor(const Color3f& color)
{
    _firstColor = color;
}

//------------------------------------------------------------------------------
// VolumeScene::SetSecondObjectMaterial
//------------------------------------------------------------------------------
inline void
VolumeScene::SetSecondObjectColor(const Color3f& color)
{
    _secondColor = color;
}

//------------------------------------------------------------------------------
// VolumeScene::SetThresholdIsovalue
//------------------------------------------------------------------------------
inline void
VolumeScene::SetThresholdIsovalue(const float thresholdIsovalue)
{
    _thresholdIsovalue = thresholdIsovalue;
}

//------------------------------------------------------------------------------
// VolumeScene::HasErrors
//------------------------------------------------------------------------------
inline bool
VolumeScene::HasErrors() const
{
    return _hasErrors;
}

//------------------------------------------------------------------------------
// VolumeScene::GetErrorString
//------------------------------------------------------------------------------
inline std::string
VolumeScene::GetErrorString() const
{
    return _errorString;
}

//------------------------------------------------------------------------------
// VolumeScene::Render
//------------------------------------------------------------------------------
inline void
VolumeScene::Render()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, (GLdouble)_imageResolution.x,
               0.0, (GLdouble)_imageResolution.y);
    glRasterPos2i(0, 0);
    glDrawPixels(_imageResolution.x, _imageResolution.y,
                 GL_RGB, GL_UNSIGNED_BYTE, _pixels);
    glFlush();
}

//------------------------------------------------------------------------------
// VolumeScene::Index
//------------------------------------------------------------------------------
inline uint32_t
VolumeScene::Index(const uint32_t i,
                   const uint32_t j,
                   const uint32_t k) const
{
    return i + _dimensions.x * j + _dimensions.y * _dimensions.x * k;
}

//------------------------------------------------------------------------------
// VolumeScene::Sample
//------------------------------------------------------------------------------
inline float
VolumeScene::Sample(const uint32_t i,
                    const uint32_t j,
                    const uint32_t k) const
{
    assert(i < _dimensions.x);
    assert(j < _dimensions.y);
    assert(k < _dimensions.z);
    if (_renderingMRI) {
        return _headData[i][j][k];
    }
    return _data[Index(i, j, k)];
}

//------------------------------------------------------------------------------
// VolumeScene::SampleTrilinear
//------------------------------------------------------------------------------
inline bool
VolumeScene::SampleTrilinear(float& value, const Vector3f& p) const
{
    return SampleTrilinear(value, p.x, p.y, p.z);
}

//------------------------------------------------------------------------------
// VolumeScene::GetGradient
//------------------------------------------------------------------------------
inline void
VolumeScene::GetGradient(Vector3f& g,
                         const uint32_t i,
                         const uint32_t j,
                         const uint32_t k) const
{
    assert(i < _dimensions.x);
    assert(j < _dimensions.y);
    assert(k < _dimensions.z);

    if (i == 0) {
        g.x = Sample(i + 1, j, k) - Sample(i, j, k);
    } else if (i == _dimensions.x - 1) {
        g.x = Sample(i, j, k) - Sample(i - 1, j, k);
    } else {
        g.x = Sample(i + 1, j, k) - Sample(i - 1, j, k);
        g.x *= 0.5f;
    }

    if (j == 0) {
        g.y = Sample(i, j + 1, k) - Sample(i, j, k);
    } else if (j == _dimensions.y - 1) {
        g.y = Sample(i, j, k) - Sample(i, j - 1, k);
    } else {
        g.y = Sample(i, j + 1, k) - Sample(i, j - 1, k);
        g.y *= 0.5f;
    }

    if (k == 0) {
        g.z = Sample(i, j, k + 1) - Sample(i, j, k);
    } else if (k == _dimensions.z - 1) {
        g.z = Sample(i, j, k) - Sample(i, j, k - 1);
    } else {
        g.z = Sample(i, j, k + 1) - Sample(i, j, k - 1);
        g.z *= 0.5f;
    }
}

//------------------------------------------------------------------------------
// VolumeScene::GetGradientTrilinear
//------------------------------------------------------------------------------
inline bool
VolumeScene::GetGradientTrilinear(Vector3f& g, const Vector3f& p) const
{
    return GetGradientTrilinear(g, p.x, p.y, p.z);
}

#endif
