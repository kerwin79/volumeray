#include "Clock.hpp"

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

//------------------------------------------------------------------------------
// Clock::Clock
//------------------------------------------------------------------------------
Clock::Clock(const float interval)
: _start(0.0)
, _previous(0.0)
, _current(0.0)
, _dt(0.0)
, _fps(0.0)
, _count(0)
, _interval(0)
#if defined(WIN32)
, _frequency()
, _counter()
#endif
{
    _start = _previous = GetTime();
    _interval = uint32_t(interval);
}

//------------------------------------------------------------------------------
// Clock::Update
//------------------------------------------------------------------------------
void Clock::Tick()
{
    _current = GetTime();
    _dt = _current - _previous;
    _previous = _current;
    _count++;

    if (_current - _start > _interval)
    {
        _fps = _count / double(_current - _start);
        _start = _current;
        _count = 0;
    }
}

//------------------------------------------------------------------------------
// Clock::DeltaTime
//------------------------------------------------------------------------------
double Clock::GetDeltaTime()
{
    return _dt;
}

//------------------------------------------------------------------------------
// Clock::GetFPS
//------------------------------------------------------------------------------
double Clock::GetFPS()
{
    return _fps;
}

//------------------------------------------------------------------------------
// Clock::GetTime
//------------------------------------------------------------------------------
double Clock::GetTime()
{
#if defined(WIN32)
    QueryPerformanceFrequency(&_frequency);
    QueryPerformanceCounter(&_counter);
    return (double(_counter.QuadPart) / double(_frequency.QuadPart));
#elif __MACH__
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    return mts.tv_sec + 1e-9 * mts.tv_nsec;
#endif
}

//------------------------------------------------------------------------------
// Clock::SetInterval
//------------------------------------------------------------------------------
void Clock::SetInterval(const float interval)
{
    _interval = uint32_t(interval);
}
