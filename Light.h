#ifndef LIGHT_H
#define LIGHT_H

#include <Eigen/Dense>

using Eigen::Vector3d;

class Light {
public:
    const double power; // light intensity
    const Vector3d source;
    Light(const Vector3d& src, double pow) : power(pow), source(src) {};
    
};

#endif // LIGHT_H









