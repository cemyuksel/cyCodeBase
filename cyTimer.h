// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyTimer.h 
//! \author Cem Yuksel
//! 
//! \brief  Timer classes to measure performance
//! 
//! This file includes timer classes to measure performance in windows systems.
//! Timer class uses Windows specific calls to measure the time,
//! therefore this file can only be used in windows platforms.
//!
//-------------------------------------------------------------------------------
//
// Copyright (c) 2016, Cem Yuksel <cem@cemyuksel.com>
// All rights reserved.
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy 
// of this software and associated documentation files (the "Software"), to deal 
// in the Software without restriction, including without limitation the rights 
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
// copies of the Software, and to permit persons to whom the Software is 
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
// SOFTWARE.
// 
//-------------------------------------------------------------------------------

#ifndef _CY_TIMER_H_INCLUDED_
#define _CY_TIMER_H_INCLUDED_

//-------------------------------------------------------------------------------

#include <time.h>
#include <math.h>

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

//! Simple stopwatch class.
//!
//! Use this class to measure the time between Start and Stop calls.

class Timer
{
public:
	//! Starts the timer
	void Start() {
		startTime = clock();
	}

	//! Returns the time passed since Start call in seconds.
	//! Note that this method does not actually stop the timer,
	//! it only returns the time passed since Start call.
	//! Therefore, you can call this method as many times as you like
	//! once you call Start method once.
	double Stop() const {
		clock_t endTime = clock();
		return float(endTime-startTime)/CLOCKS_PER_SEC;;
	}

protected:
	clock_t startTime;	//!< Keeps the starting time
};

//-------------------------------------------------------------------------------

//! Stopwatch class with statistics.
//!
//! Use this class to measure the time between Start and Stop calls.
//! Unlike Timer class, this class also provides statistical information.

class TimerStats
{
public:
	TimerStats() { Clear(); }

	//!@name Timer Methods

	//! Starts the timer
	void Start() {
		timer.Start();
	}

	//! Stops the timer and records the current measurement.
	//! Returns the time passed since Start call in seconds.
	double Stop() {
		double t = timer.Stop();
		unsigned char p = pos & 0x7F;
		totalTime += t - times[ p ];
		times[ p ] = t;
		if ( minTime > t ) minTime = t;
		if ( maxTime < t ) maxTime = t;
		// Increment pos such that the first bit shows if times array is full.
		p++;
		pos = ( p & 0x7F ) + ( ( pos | p ) & 0x80 );
		return t;
	}


	//!@name Statistics Methods

	//! Clears all the time records
	void Clear() {
		for ( int i=0; i<128; i++ ) times[i] = 0;
		minTime = 1.0e30;
		maxTime = 0;
		totalTime = 0;
		pos = 0;
	}

	//! Returns the last measured time.
	//! If no time is measured before, returns zero.
	double GetLastTime() const { return times[ (pos-1) & 0x7F ]; }

	//! Returns the minimum measured time
	double GetMin() const { return minTime; }

	//! Returns the maximum measured time
	double GetMax() const { return maxTime; }

	//! Returns the average of all recorded times (max previous 128 records)
	double GetAverage() const { return totalTime / (double) GetRecordCount(); }

	//! Returns the variance of the time records.
	//! Note that this method goes over all time records,
	//! so it may take a little time to compute the variance.
	double GetVariance() const {
		unsigned char count = GetRecordCount();
		double avrg = totalTime / (double) count;
		double var = 0;
		for ( int i=0; i<count; i++ ) var += (times[i]-avrg) * (times[i]-avrg);
		return var / (double) count;
	}

	//! Returns the standard deviation of the time records.
	//! Note that this method goes over all time records,
	//! so it may take a little time to compute the standard deviation.
	double GetStdev() const { return sqrt( GetVariance() ); }


	//!@name Access time records

	//! Returns the number of time records.
	unsigned char GetRecordCount() const { return ( pos & 0x80 ) ? 128 : pos; }

	//! Returns the array of time records.
	const double* GetRecords() const { return times; }


protected:
	//! \internal
	Timer   timer;		// The timer
	double	times[128];	// List of measured times
	double	minTime;	// Minimum time recorded
	double	maxTime;	// Maximum time recorded
	double	totalTime;	// Total of recorded times
	unsigned char pos;	// Position of the next time record
};

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Timer      cyTimer;			//!< Simple stopwatch class
typedef cy::TimerStats cyTimerStats;	//!< Stopwatch class with statistics

//-------------------------------------------------------------------------------

#endif
