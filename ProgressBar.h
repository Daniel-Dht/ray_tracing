#ifndef PROGBAR_H
#define PROGBAR_H

#include<iostream>
#include <chrono> 

using namespace std;

class ProgressBar {
public:
    int i=0;
    int N;
    int barWidth  = 60;
    int modulo;     // step between each call of cout_progress()
    int count = -1;  // number of call
    int total_call; // number of time cout_progress() will be called
    std::chrono::time_point<std::chrono::system_clock> start, end;
    float t_estimation = 0 ; // etimation time of the raytracing algorithm

    ProgressBar(int N, int m) : N(N), modulo(m) {        
        cout_progress(0);  
        start = std::chrono::system_clock::now();
        total_call = N/modulo + 1;
    };

    void cout_progress(float progress){

        int pos = barWidth * progress/N;

        cout << "[";
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) cout << "=";
            else if (i == pos) cout << ">";
            else cout << " ";
        }        
        count++;
        // calculate an estimation of total time remaining
        if(count>0){
            end = std::chrono::system_clock::now();
            float elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds> (end-start).count();
            t_estimation = elapsed_seconds / total_call * float(total_call-count);
        }
        cout << "] " << progress/N * 100.0 <<"%\r"; //, remaining:" << (int) t_estimation << " seconds
        cout.flush();
    }
};

#endif // PROGBAR_H