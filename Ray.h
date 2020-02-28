
#ifndef RAY_H
#define RAY_H

#include "Vector.h"

class Ray {
public:
    Vector3d C, u;
    Ray(const Vector3d& C, Vector3d u) : C(C), u(u) {}; // Origin et direction    

    Vector3d computeIntersection(const double& t) const { 
        // return the intersection point
        
        return C + u*t;           
    }
};

struct InterStruct{
public:
    bool computeP = true;
    bool computeN = true;
    Vector3d N = Vector3d(0);
    Vector3d P = Vector3d(0);

};  
#endif // RAY_H
