//==============================================================================
// Volume Rendering by Raycasting
// Tim Thirion
// CS 530 - Purdue University
// October 2004
//
//    File: Main.cpp
//  Author: Tim Thirion
// Created: October 7, 2004
// Revised: October 13, 2004
//
// License: You may redistribute this file as long as the author's name and
//          this header appears in its entirety. If this source is used in a
//          binary distribution, this source file must also be included in said
//          distribution. The author is not responsible for any deleterious
//          effects this code may have on your computer.
//
//   Notes: This file is the entry point for this program. The original source
//          file was taken from the first CS 530 project (Marching Cubes).
//==============================================================================

// Library Inclusions
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>

// Project inclusions
#include "Clock.hpp"
#include "VolumeScene.hpp"

//------------------------------------------------------------------------------
// Global Variables
//------------------------------------------------------------------------------
Clock       g_clock;
bool        g_renderFunction(false);
bool        g_renderDiagonalView(false);
int         g_windowWidth(600);
int         g_windowHeight(600);
Trivariatef g_density;
float       g_a(0.0f);
float       g_b(0.0f);
float       g_r0(0.0f);
VolumeScene g_volume;

//------------------------------------------------------------------------------
// Function Prototypes
//------------------------------------------------------------------------------
float Rho(const float x, const float y, const float z);
void RenderFunction();
void RenderMRIData();
void Initialize();
void Display();
void Keyboard(unsigned char key, int x, int y);
void Redisplay();
void Shutdown();

//------------------------------------------------------------------------------
// Rho
//------------------------------------------------------------------------------
float Rho(const float x, const float y, const float z)
{
    const Vector3f A(-0.5f, -0.5f, -0.5f);
    const Vector3f B( 0.5f,  0.5f, -0.5f);
    const Vector3f C( 0.5f, -0.5f,  0.5f);
    const Vector3f D(-0.5f,  0.5f,  0.5f);

    Vector3f r(x, y, z);
    return g_a * expf(-Length(r) / g_r0)
         + g_b * expf(-Length(r - A) / g_r0)
         + g_b * expf(-Length(r - B) / g_r0)
         + g_b * expf(-Length(r - C) / g_r0)
        + g_b * expf(-Length(r - D) / g_r0);
//  return g_a * FloatUtilf::Exp(-Length(r) / g_r0)
//       + g_b * FloatUtilf::Exp(-Length(r - A) / g_r0)
//       + g_b * FloatUtilf::Exp(-Length(r - B) / g_r0)
//       + g_b * FloatUtilf::Exp(-Length(r - C) / g_r0)
//       + g_b * FloatUtilf::Exp(-Length(r - D) / g_r0);
}

//------------------------------------------------------------------------------
// RenderFunction
//------------------------------------------------------------------------------
void RenderFunction()
{
    g_a = 80.0f;
    g_b = 50.0f;
    g_r0 = 0.5f;

    g_volume.SetDimensions(Vector3ui(128, 128, 128));
    g_volume.DataFromFunction(g_density);

    //
    // Set illumination parameters.
    //
    g_volume.SetAmbientLight(Color3f(0.1f, 0.1f, 0.1f));
    g_volume.SetAmbientContribution(1.0f);
    g_volume.SetDiffuseContribution(0.4f);
    g_volume.SetSpecularContribution(0.6f);
    g_volume.SetPhongPower(20.0f);
    g_volume.SetConstantAttenuation(1.0f);
    g_volume.SetLinearAttenuation(0.01f);

    g_volume.SetFirstObjectColor(Color3f(1.0f, 1.0f, 1.0f));
    g_volume.SetSecondObjectColor(Color3f(0.0f, 0.0f, 1.0f));

    //
    // Set ray parameters.
    //
    g_volume.SetNumResamples(60);
    g_volume.SetRayStepSize(9.04);

    //
    // Set image parameters.
    //
    g_volume.SetImageResolution(g_windowWidth, g_windowHeight);

    if (g_renderDiagonalView) {
        g_volume.SetEyePosition(Vector3f(0.0f, 63.5f, 0.0f));
        g_volume.SetPointLight(Vector3f(0.0f, 63.5f, 0.0f),
                               Color3f(3.0f, 3.0f, 3.0f));
        g_volume.LookAt(Vector3f(1.0f, 0.0f, 1.0f));

        g_volume.SetImageWorldBounds(Vector2f(0.0f, 0.0f),
                                     Vector2f(180.0f, 180.0f));
        g_volume.SetImageBasis(Vector3f(1.0f, 0.0f, -1.0f),
                               Vector3f(0.0f, 1.0f, 0.0f));
    } else {
        g_volume.SetEyePosition(Vector3f(63.5f, 63.5f, -10.0f));
        g_volume.SetPointLight(Vector3f(63.5f, 63.5f, -10.0f),
                               Color3f(3.0f, 3.0f, 3.0f));
        g_volume.LookAt(Vector3f(0.0f, 0.0f, 1.0f));

        g_volume.SetImageWorldBounds(Vector2f(0.0f, 0.0f),
                                     Vector2f(127.0f, 127.0f));

        g_volume.SetImageBasis(Vector3f(1.0f, 0.0f, 0.0f),
                               Vector3f(0.0f, 1.0f, 0.0f));
    }

    //
    // Set transfer function parameters and isovalues.
    //
    g_volume.SetIsovalues(Vector2f(35.0f, 80.0f));
    g_volume.SetThresholdIsovalue(60.0f);

    g_volume.SetThicknesses(Vector2f(5.0f, 10.0f));
    g_volume.SetPrefactors(Vector2f(0.05f, 0.8f));
    g_volume.SetGradientPowers(Vector2f(0.0f, 0.0f));
}

//------------------------------------------------------------------------------
// RenderMRIData
//------------------------------------------------------------------------------
void RenderMRIData()
{
    if (!g_volume.DataFromMRIFile("3DHEAD.vol")) {
        std::cout << g_volume.GetErrorString() << std::endl;
        std::cout << "Press any key to continue." << std::endl;
        char c;
        std::cin >> c;
        exit(-1);
    }

    //
    // Set illumination parameters.
    //
    g_volume.SetAmbientLight(Color3f(0.4f, 0.4f, 0.4f));
    g_volume.SetAmbientContribution(1.0f);
    g_volume.SetDiffuseContribution(0.5f);
    g_volume.SetSpecularContribution(0.5f);
    g_volume.SetPhongPower(20.0f);
    g_volume.SetConstantAttenuation(1.0f);
    g_volume.SetLinearAttenuation(0.01f);

    g_volume.SetFirstObjectColor(Color3f(0.0f, 1.0f, 0.0f));
    g_volume.SetSecondObjectColor(Color3f(1.0f, 1.0f, 1.0f));

    //
    // Set ray parameters.
    //
    g_volume.SetRayStepSize(9.04f);

    //
    // Set image parameters.
    //
    g_volume.SetImageResolution(g_windowWidth, g_windowHeight);

    //
    // For now assume z = eye.z. Solve the problem of having
    // an arbitrarily oriented image plane.
    //
    if (g_renderDiagonalView) {
        g_volume.SetEyePosition(Vector3f(0.0f, 127.5f, 0.0f));
        g_volume.SetPointLight(Vector3f(0.0f, 127.5f, 0.0f),
                               Color3f(3.0f, 3.0f, 3.0f));
        g_volume.LookAt(Vector3f(1.0f, 0.0f, 1.0f));

        g_volume.SetImageWorldBounds(Vector2f(0.0f, 0.0f),
                                     Vector2f(260.0f, 260.0f));
        g_volume.SetImageBasis(Vector3f(1.0f, 0.0f, -1.0f),
                               Vector3f(0.0f, 1.0f, 0.0f));
    } else {
        g_volume.SetEyePosition(Vector3f(127.5f, 127.5f, -10.0f));
        g_volume.SetPointLight(Vector3f(127.0f, 127.0f, -10.0),
                               Color3f(3.0f, 3.0f, 3.0f));
        g_volume.LookAt(Vector3f(0.0f, 0.0f, 1.0f));

        g_volume.SetImageWorldBounds(Vector2f(0.0f, 0.0f),
                                     Vector2f(255.0f, 255.0f));

        g_volume.SetImageBasis(Vector3f(1.0f, 0.0f, 0.0f),
                               Vector3f(0.0f, 1.0f, 0.0f));
    }


    //
    // Set transfer function parameters and isovalues.
    //
    g_volume.SetIsovalues(Vector2f(50.0f, 90.0f));
    g_volume.SetThresholdIsovalue(70.0f);

    g_volume.SetThicknesses(Vector2f(10.0f, 20.0f));
    g_volume.SetPrefactors(Vector2f(0.1f, 0.1f));
    g_volume.SetGradientPowers(Vector2f(0.0f, 0.0f));
}

//------------------------------------------------------------------------------
// Initialize
//------------------------------------------------------------------------------
void Initialize()
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glDisable(GL_DEPTH_TEST);

    //std::cout << "The Density Function is:\n\n"
    //            << "p(r) =   a * f(|r|)\n"
    //               "       + b * f(|r - A|) + b * f(|r - B|)\n"
    //               "       + b * f(|r - C|) + b * f(|r - D|)\n"
    //            << "Where f(x) = exp(-x/r0).\n";

    g_density.SetFunction(Rho);

    if (g_renderFunction) {
        RenderFunction();
    } else {
        RenderMRIData();
    }

    std::cout << "The image is "
              << g_windowWidth << "x" << g_windowHeight
              << " pixels.\n";
    std::cout << "One dot will be printed for each "
                 "row of pixels determined.\n";

    //
    // Build the image; record the time taken.
    //
    g_clock.Tick();
    g_volume.BuildImage();
    g_clock.Tick();
    double renderTime = g_clock.GetDeltaTime();
    std::cout << "\nImage rendered in "
              << renderTime << " seconds."
              << std::endl;
}

//------------------------------------------------------------------------------
// Display
//------------------------------------------------------------------------------
void Display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    g_volume.Render();
    glutSwapBuffers();
}

//------------------------------------------------------------------------------
// Keyboard
//------------------------------------------------------------------------------
void Keyboard(unsigned char key, int, int)
{
    switch (key) {
        case 27: // Escape
        {
            exit(0);
        } break;
        default:
            break;
    }
}

//------------------------------------------------------------------------------
// Redisplay
//------------------------------------------------------------------------------
void Redisplay()
{
    glutPostRedisplay();
}

//------------------------------------------------------------------------------
// Shutdown
//------------------------------------------------------------------------------
void Shutdown()
{
}

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(g_windowWidth, g_windowHeight);
    glutCreateWindow("volumeray");
    Initialize();
    glutDisplayFunc(Display);
    glutKeyboardFunc(Keyboard);
    glutIdleFunc(Redisplay);
    glutMainLoop();
    Shutdown();
    return 0;
}
