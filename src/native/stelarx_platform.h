#pragma once

// Small host-side portability layer shared by the CUDA JNI libraries. Device
// code remains unchanged; these helpers only cover timing, sleeping, and TTY
// detection used by progress reporting.

#include <chrono>
#include <thread>

#ifdef _WIN32
#  include <io.h>
#else
#  include <unistd.h>
#endif

static inline double stelarx_now_sec() {
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

static inline void stelarx_sleep_millis(int millis) {
    std::this_thread::sleep_for(std::chrono::milliseconds(millis));
}

static inline int stelarx_stderr_isatty() {
#ifdef _WIN32
    return _isatty(_fileno(stderr));
#else
    return isatty(STDERR_FILENO);
#endif
}
