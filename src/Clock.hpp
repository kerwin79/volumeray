#ifndef __Clock_h__
#define __Clock_h__

#if defined(WIN32)
#define WIN32_LEANANDMEAN
#include <windows.h>
#endif

#include "Common.hpp"

class Clock
{
public:
    Clock(const float interval = 1.0f);

    void SetInterval(const float interval);

    void Tick();
    double GetTime();

    double GetDeltaTime();
    double GetFPS();

private:
    double _start;
    double _previous;
    double _current;
    double _dt;
    double _fps;
    uint32_t _count;
    uint32_t _interval;

#if defined(WIN32)
    LARGE_INTEGER _frequency, _counter;
#endif
};

#endif
