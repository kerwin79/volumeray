#include <fstream>
#include <iostream>
#include "Ray.hpp"
#include "VolumeScene.hpp"

//------------------------------------------------------------------------------
// VolumeScene::DataFromFile
//------------------------------------------------------------------------------
bool
VolumeScene::DataFromMRIFile(const std::string& fileName, const bool isUnix)
{
    _renderingMRI = true;
    std::fstream fin(fileName.c_str(), std::ios::in | std::ios::binary);

    if (!fin.good()) {
        _errorString = "Could not open specified file: " + fileName + "\n";
        _hasErrors = true;
        return false;
    }

    int magicNumber     = 0;
    int headerLength    = 0;
    int sliceWidth      = 0;
    int sliceHeight     = 0;
    int numberOfSlices  = 0;
    int bitsPerPixel    = 0;
    int indexBits       = 0;

    fin.read(reinterpret_cast<char*>(&magicNumber), sizeof(int));
    fin.read(reinterpret_cast<char*>(&sliceWidth), sizeof(int));
    fin.read(reinterpret_cast<char*>(&sliceHeight), sizeof(int));
    fin.read(reinterpret_cast<char*>(&numberOfSlices), sizeof(int));
    fin.read(reinterpret_cast<char*>(&bitsPerPixel), sizeof(int));
    fin.read(reinterpret_cast<char*>(&indexBits), sizeof(int));

    if (isUnix) {
        magicNumber     = LongSwap(magicNumber);
        headerLength    = LongSwap(headerLength);
        sliceWidth      = LongSwap(sliceWidth);
        sliceHeight     = LongSwap(sliceHeight);
        numberOfSlices  = LongSwap(numberOfSlices);
        bitsPerPixel    = LongSwap(bitsPerPixel);
        indexBits       = LongSwap(indexBits);
    }

    SetDimensions(Vector3ui(sliceHeight, numberOfSlices, bitsPerPixel));
    assert(_dimensionsProduct > 0);

    _numResamples = sliceWidth;

    fin.read(reinterpret_cast<char*>(&_scale.x), 3*sizeof(float));
    fin.read(reinterpret_cast<char*>(&_rotationAngles.x), 3*sizeof(float));

    if (isUnix) {
        _scale.x = FloatSwap(_scale.x);
        _scale.y = FloatSwap(_scale.y);
        _scale.z = FloatSwap(_scale.z);
        _rotationAngles.x = FloatSwap(_rotationAngles.x);
        _rotationAngles.y = FloatSwap(_rotationAngles.y);
        _rotationAngles.z = FloatSwap(_rotationAngles.z);
    }

    std::cout << "Populating data from the file";

    _widths.x = _boundingBox.GetExtent(0) / (_dimensions.x - 1);
    _widths.y = _boundingBox.GetExtent(1) / (_dimensions.y - 1);
    _widths.z = _boundingBox.GetExtent(2) / (_dimensions.z - 1);

    float maxValue(0.0f);
    float minValue(INFINITY);
    float averageValue(0.0f);
    float sum(0.0f);
    float temp(0.0f);

    uint8_t c = 0;
    uint32_t i, j, k;
    for (k = 0; k < _dimensions.z; ++k) {
        for (j = 0; j < _dimensions.y; ++j) {
            for (i = 0; i < _dimensions.x; ++i) {
                fin.read(reinterpret_cast<char*>(&c), sizeof c);
                _headData[i][j][k] = float(c);

                temp = _headData[i][j][k];
                if (temp > maxValue) {
                    maxValue = temp;
                }
                if (temp < minValue) {
                    minValue = temp;
                }
                sum += temp;
            }
        }
    }

    averageValue = sum / _dimensionsProduct;

    std::cout << "\nMaximum: " << maxValue << "\n";
    std::cout << "Minimum: " << minValue << "\n";
    std::cout << "Average: " << averageValue << "\n";
    std::cout << "Done.\n" << std::endl;
    return true;
}

//------------------------------------------------------------------------------
// VolumeScene::DataFromFunction
//------------------------------------------------------------------------------
void
VolumeScene::DataFromFunction(const Trivariatef& f)
{
    _renderingMRI = false;
    std::cout << "Populating data from the function.";

    assert(_dimensionsProduct > 0);
    _data = new float[_dimensionsProduct];

    Vector3f d;
    d.x = 2.0f / _dimensions.x;
    d.y = 2.0f / _dimensions.y;
    d.z = 2.0f / _dimensions.z;

    float maxValue(0.0f);
    float minValue(INFINITY);
    float averageValue(0.0f);
    float sum(0.0f);
    float temp(0.0f);

    Vector3f  p;
    uint32_t index = 0;
    for (uint32_t i = 0; i < _dimensions.x; ++i) {
        for (uint32_t j = 0; j < _dimensions.y; ++j) {
            for (uint32_t k = 0; k < _dimensions.z; ++k) {
                index = Index(i, j, k);
                p.x = -1.0f + d.x * i;
                p.y = -1.0f + d.y * j;
                p.z = -1.0f + d.z * k;

                _data[index] = f(p.x, p.y, p.z);

                temp = _data[index];
                if (temp > maxValue) {
                    maxValue = temp;
                }
                if (temp < minValue) {
                    minValue = temp;
                }
                sum += temp;
            }
        }
    }

    averageValue = sum / _dimensionsProduct;

    std::cout << "\nMaximum: " << maxValue;
    std::cout << "\nMinimum: " << minValue;
    std::cout << "\nAverage: " << averageValue;
    std::cout << "\nDone.\n" << std::endl;
}

//------------------------------------------------------------------------------
// VolumeScene::BuildImage
//------------------------------------------------------------------------------
void
VolumeScene::BuildImage()
{
    std::cout << "Building image, please wait.\n";

    //
    // Associate each pixel with a ray.
    //
    Rayf r;
    r.d = _eyeDirection;
    Vector3f position;

    //
    // Color information
    //
    Color3f     illuminatedColor;
    Color3f     accumulatedColor;
    Color3ub    pixelColor;

    //
    // Opacity information
    //
    float       currentOpacity = 0.0f;
    float       accumulatedOpacity = 0.0f;

    //
    // Determine the color of each pixel via raycasting.
    //
    float maxZ = 0.0f;

    Vector2f pixelWidths = _imageMaxes - _imageMins;
    pixelWidths.x /= float(_imageResolution.x);
    pixelWidths.y /= float(_imageResolution.y);

    Vector2f worldImageWidths;
    worldImageWidths.x = Abs(_imageMaxes.x - _imageMins.x);
    worldImageWidths.y = Abs(_imageMaxes.y - _imageMins.y);

    Vector3f imageWorldMins = _eyePosition
                              - _x * (worldImageWidths.x * 0.5f)
                              - _y * (worldImageWidths.y * 0.5f);

    const uint32_t numPixels = _imageResolution.x * _imageResolution.y;
    for (uint32_t i = 0; i < numPixels; ++i) {
        //
        // Compute the origin of the ray.
        //
        Vector2ui pixelIndex(i % _imageResolution.x,
                             i / _imageResolution.y);

        if (pixelIndex.x == 0) {
            std::cout << ".";
        }

        float x = (pixelIndex.x + 0.5f) * pixelWidths.x + imageWorldMins.x;
        float y = (pixelIndex.y + 0.5f) * pixelWidths.y + imageWorldMins.y;

        r.o = _x * x + _y * y;
        position = r.o;

        //
        // For each ray, calculate color and opacity.
        //
        for (uint32_t j = 0; j < _numResamples; ++j) {
            if (accumulatedOpacity >= 1.0f)
                break;

            //
            // Compute the current opacity.
            //
            currentOpacity = Alpha(position);

            //
            // Compute the new color (skip if the opacity is negligible).
            //
            if (currentOpacity > 0.01f) {
                ComputeColorByIllumination(illuminatedColor, position, r.o);
                accumulatedColor = illuminatedColor * currentOpacity *
                                   (1 - accumulatedOpacity) + accumulatedColor;
            }

            //
            // Accumulate the opacity.
            //
            accumulatedOpacity = (1 - accumulatedOpacity) * currentOpacity +
                                 accumulatedOpacity;

            //
            // Advance the ray.
            //
            position += r.d * _rayStepSize;

            if (position.z > maxZ) {
                maxZ = position.z;
            }
        }
        accumulatedColor.Clamp();
        Convert(pixelColor, accumulatedColor);
        _pixels[i] = pixelColor;

        accumulatedColor = Color3f(0.0f, 0.0f, 0.0f);
        accumulatedOpacity = 0.0f;

        if (_hasErrors) {
            std::cout << _errorString << std::endl;
            return;
        }
    }

    std::cout << "\n\nMaximum depth: " << maxZ << std::endl;
}

//------------------------------------------------------------------------------
// VolumeScene::SampleTrilinear
//------------------------------------------------------------------------------
bool
VolumeScene::SampleTrilinear(float& value,
                             const float x,
                             const float y,
                             const float z) const
{
    if (x < 0.0f || x >= _dimensions.x - 1)
        return false;
    if (y < 0.0f || y >= _dimensions.y - 1)
        return false;
    if (z < 0.0f || z >= _dimensions.z - 1)
        return false;

    int i0 = int(x / _widths.x);
    float xLerp = (x - float(i0)) / _widths.x;

    int j0 = int(y / _widths.y);
    float yLerp = (y - float(j0)) / _widths.y;

    int k0 = int(z / _widths.z);
    float zLerp = (z - float(k0)) / _widths.z;

    uint32_t i = uint32_t(i0);
    uint32_t j = uint32_t(j0);
    uint32_t k = uint32_t(k0);

    float f[8];
    f[0] = Sample(    i,     j,     k);
    f[1] = Sample(i + 1,     j,     k);
    f[2] = Sample(i + 1,     j, k + 1);
    f[3] = Sample(    i,     j, k + 1);
    f[4] = Sample(    i, j + 1,     k);
    f[5] = Sample(i + 1, j + 1,     k);
    f[6] = Sample(i + 1, j + 1, k + 1);
    f[7] = Sample(    i, j + 1, k + 1);

    Lerp(f[0], f[0], f[1], xLerp);
    Lerp(f[1], f[3], f[2], xLerp);
    Lerp(f[2], f[7], f[6], xLerp);
    Lerp(f[3], f[4], f[5], xLerp);

    Lerp(f[0], f[0], f[1], zLerp);
    Lerp(f[1], f[3], f[2], zLerp);

    Lerp(value, f[0], f[1], yLerp);

    //FloatUtilf::Lerp(f[0], f[0], f[1], xLerp);
    //FloatUtilf::Lerp(f[1], f[3], f[2], xLerp);
    //FloatUtilf::Lerp(f[2], f[7], f[6], xLerp);
    //FloatUtilf::Lerp(f[3], f[4], f[5], xLerp);

    //FloatUtilf::Lerp(f[0], f[0], f[1], zLerp);
    //FloatUtilf::Lerp(f[1], f[3], f[2], zLerp);

    //FloatUtilf::Lerp(value, f[0], f[1], yLerp);
    return true;
}

//------------------------------------------------------------------------------
// VolumeScene::GetGradientTrilinear
//------------------------------------------------------------------------------
bool
VolumeScene::GetGradientTrilinear(Vector3f& gradient,
                                  const float x,
                                  const float y,
                                  const float z) const
{
    if (x < 0.0f || x >= _dimensions.x - 1)
        return false;
    if (y < 0.0f || y >= _dimensions.y - 1)
        return false;
    if (z < 0.0f || z >= _dimensions.z - 1)
        return false;

    int i0 = int(x / _widths.x);    
    float xLerp = (x - float(i0)) / _widths.x;

    int j0 = int(y / _widths.y);
    float yLerp = (y - float(j0)) / _widths.y;

    int k0 = int(z / _widths.z);
    float zLerp = (z - float(k0)) / _widths.z;

    uint32_t i = uint32_t(i0);
    uint32_t j = uint32_t(j0);
    uint32_t k = uint32_t(k0);

    Vector3f g[8];
    GetGradient(g[0],     i,     j,     k);
    GetGradient(g[1], i + 1,     j,     k);
    GetGradient(g[2], i + 1,     j, k + 1);
    GetGradient(g[3],     i,     j, k + 1);
    GetGradient(g[4],     i, j + 1,     k);
    GetGradient(g[5], i + 1, j + 1,     k);
    GetGradient(g[6], i + 1, j + 1, k + 1);
    GetGradient(g[7],     i, j + 1, k + 1);

    Lerp(g[0], g[0], g[1], xLerp);
    Lerp(g[1], g[3], g[2], xLerp);
    Lerp(g[2], g[7], g[6], xLerp);
    Lerp(g[3], g[4], g[5], xLerp);

    Lerp(g[0], g[0], g[1], zLerp);
    Lerp(g[1], g[3], g[2], zLerp);

    Lerp(gradient, g[0], g[1], yLerp);
    return true;
}

//------------------------------------------------------------------------------
//  VolumeScene::Alpha
//------------------------------------------------------------------------------
float
VolumeScene::Alpha(const Vector3f& r)
{
    float f;

    if (SampleTrilinear(f, r)) {
        SampleTrilinear(f, r);
        Vector3f g;
        Vector2f t = _thicknesses;

        GetGradientTrilinear(g, r);
        float gLength = Length(g);
        float alpha1 = _prefactors[0] *
                       powf(gLength, _gradientPowers[0]);

        alpha1 *= expf(-((f - _isovalues[0]) * (f - _isovalues[0]))
                  / (2.0f * (t[0] * t[0])));

        float alpha2 = _prefactors[1] *
                       powf(gLength, _gradientPowers[1]);

        alpha2 *= expf(-((f - _isovalues[1]) * (f - _isovalues[1]))
                  / (2.0f * (t[1] * t[1])));

        return alpha1 + alpha2;
    }
    return 0.0f;
}

//------------------------------------------------------------------------------
// VolumeScene::ComputeColorByIllumination
//------------------------------------------------------------------------------
void
VolumeScene::ComputeColorByIllumination(Color3f& color,
                                        const Vector3f& position,
                                        const Vector3f& rayOrigin)
{
    //
    // Set initial color to zero.
    //
    color = Color3f(0, 0, 0);

    //
    // Compute the unit normal vector.
    //
    Vector3f N;
    if (!GetGradientTrilinear(N, position)) {
        _hasErrors = true;
        _errorString += "\nAttempted to compute gradient "
                        "outside of voxel grid.";
        return;
    }
    N = -N;
    Normalize(N);

    //
    // Compute the view vector.
    //
    Vector3f V = _eyePosition - position;
    Normalize(V);

    //
    // If the value at the current position is greater than the threshold
    // isovalue, use the second object material properties.
    //
    float f(0.0f);
    SampleTrilinear(f, position);

    Color3f objectColor = _firstColor;
    if (f > _thresholdIsovalue) {
        objectColor = _secondColor;
    }

    //
    // Compute the light vector.
    //
    Vector3f L = _lightPosition - position;
    Normalize(L);

    //
    // Compute the ambient contribution.
    //
    color += objectColor * _ambientLight * _Ka;

    //
    // Compute the diffuse contribution.
    //
    float diffuse = _Kd * Max(Dot(N, L), 0.0f);

    //
    // Compute the half vector.
    //
    Vector3f H = L + V;
    Normalize(H);

    //
    // Compute the specular contribution.
    //
    float specular;
    if (diffuse < 0.01f) {
        specular = 0.0f;
    } else {
        specular = _Ks * Max(Dot(N, H), 0.0f);
        specular = powf(specular, _phongPower);
    }

    //
    // Attenuate diffuse and specular.
    //
    float distanceToImage = Length(position - rayOrigin);
    float attenuation = 1.0f;
    attenuation /= (_linearAtten * distanceToImage + _constantAtten);
    color += objectColor * _lightColor * attenuation * (specular + diffuse);
}
