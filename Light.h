#ifndef LIGHT_H
#define LIGHT_H

#include <Eigen/Dense>

using Eigen::Vector3d;

class Light {
public:
    const double power = 400000.0; // light intensity
    const Vector3d source = Vector3d(8.0, 0.0, 45.0);
    // const Vector3d source = Vector3d(0.0, 0.0, 20.0);
};

#endif // LIGHT_H









