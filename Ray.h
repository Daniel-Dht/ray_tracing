
#ifndef RAY_H
#define RAY_H

#include <Eigen/Dense>
using Eigen::Vector3d;

class Ray {
public:
    Vector3d C, u;
    Ray(const Vector3d& C, Vector3d u) : C(C), u(u) {}; // Centre et direction    
};

#endif // RAY_H
