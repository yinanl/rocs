/**
 *  timer.hpp
 *
 *  A timer class.
 *
 *  Created by Yinan Li on Mar. 28, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _timer_h
#define _timer_h


#include <chrono>
#include <iostream>


class timer
{
private:
    // Type aliases to make accessing nested type easier
    using clock_t = std::chrono::high_resolution_clock;
    using second_t = std::chrono::duration<double, std::ratio<1> >;
	
    std::chrono::time_point<clock_t> m_beg;
 
public:
    timer() : m_beg(clock_t::now())
    {
    }
	
    void reset()
    {
	m_beg = clock_t::now();
    }
	
    double elapsed() const
    {
	return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
    }
};

#endif
