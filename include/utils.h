#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

vector<int> get_ticks(int start, int end, int num_ticks) {
    vector<int> ticks;
    int step = (end - start) / (num_ticks - 1);
    for (int i = start; i <= end; i += step)
        ticks.push_back(i);
    return ticks;
}

vector<std::string> get_ticklabels(double start, double end, int num_ticks) {
    vector<std::string> ticklabels;
    double step = (end - start) / (num_ticks - 1);
    for (double i = start; i <= end; i += step) {
        stringstream ss;
        ss << i;
        ticklabels.push_back(ss.str());
    }
    return ticklabels;
}

class Timer {
   private:
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point _start = Clock::now();
    Clock::time_point _end = Clock::now();

   public:
    void tick() { _start = Clock::now(); }

    void tock() { _end = Clock::now(); }

    double duration() const { return std::chrono::duration<double>(_end - _start).count(); }
};

#endif