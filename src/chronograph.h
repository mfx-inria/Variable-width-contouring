/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2015
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#ifndef CHRONOGRAPH_H
#define CHRONOGRAPH_H
#endif // CHRONOGRAPH_H

#include <ctime>
#include <string>
#include <iostream>

#ifdef __APPLE__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#elif defined(__linux__)
  #include <time.h>
#endif

class Chronograph {
  public:
    Chronograph(const std::string & s) : acc_time_(0.0), name_(s), running(false) {
#ifdef __APPLE__
      mach_timebase_info(&timebase_info_);
#elif defined(__linux__)
      clock_getres(CLOCK_MONOTONIC, &time_res);
#endif
    }
    Chronograph() : acc_time_(0.0), name_(""), running(false) {
#ifdef __APPLE__
      mach_timebase_info(&timebase_info_);
#elif defined(__linux__)
      clock_getres(CLOCK_MONOTONIC, &time_res);
#endif
    }

#ifdef __APPLE__
    inline uint64_t now() const {
        return mach_absolute_time();
    }
#elif defined(__linux__)
    inline uint64_t now() const {
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return static_cast<uint64_t>(ts.tv_sec)*1000000000 + ts.tv_nsec;
    }
#endif

    void start() {
      if( running ) return;
      running = true;
      absolute_ = now();
    }

    void stop() {
      if( ! running ) return;
      acc_time_ += now() - absolute_;
      running = false;
    }

    double elapsed_time() const {
      uint64_t t;
      if( running ) {
        t = acc_time_ + now() - absolute_;
      } else {
        t = acc_time_;
      }
#ifdef __APPLE__
      return 1e-9 * (t * timebase_info_.numer / timebase_info_.denom);
#elif defined(__linux__)
      return 1e-9 * t;
#endif
    }

    template< class Output >
      void print(Output & o) {
        o << name_ << " [" << elapsed_time() << " s.]";
      }

    void reset() {
      acc_time_ = 0;
      absolute_ = 0;
      running = false;
    }

    void restart() {
      reset();
      start();
    }

  private:
    uint64_t acc_time_;
    std::string name_; // not a const reference otherwise scoped chronograph crashes
    uint64_t absolute_;
    bool running;
#ifdef __APPLE__
    struct mach_timebase_info timebase_info_;
#elif defined(__linux__)
    mutable struct timespec ts;
    struct timespec time_res;
#endif
};

template< typename OutputStream >
class ScopedChronograph : public Chronograph {
  public:
    ScopedChronograph(OutputStream & o) : Chronograph(), out_(o) { Chronograph::start(); }
    ScopedChronograph(OutputStream & o, const std::string & s) : Chronograph(s), out_(o) { Chronograph::start(); }
    ~ScopedChronograph() {
      Chronograph::stop();
      Chronograph::print(out_);
    }
  private:
    OutputStream & out_;
};
