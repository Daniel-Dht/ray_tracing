
#ifndef RAY_H
#define RAY_H

#include "Vector.h"

class Ray {
public:
    Vector3d C, u;
    Ray(const Vector3d& C, Vector3d u) : C(C), u(u) {}; // Origin et direction    
};

#endif // RAY_H
