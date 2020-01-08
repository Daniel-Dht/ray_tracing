#ifndef LIGHT_H
#define LIGHT_H

#include <Eigen/Dense>

using Eigen::Vector3d;

class Light {
public:
    const double power = 20000.0; // light intensity
    const Vector3d source = Vector3d(5.0, 5.0, 5.0);
};

#endif // LIGHT_H