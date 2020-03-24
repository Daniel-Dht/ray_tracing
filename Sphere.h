#ifndef SPHERE_H
#define SPHERE_H

#include "Ray.h"
#include "Material.h"

#include "Object.h"
#include "Vector.h"
#include "iostream" 
#include <limits>


using namespace std;

class Sphere : public Object { 
public:
    Vector3d O;
    double R;    
   
    Sphere(const Vector3d O, double R, Material mat) : O(O), R(R), Object(mat) { }; 
    
    double intersect(const Ray& r, int& total_ray_casted, InterStruct &interStruct) const {
     
        total_ray_casted++;
        Vector3d diff = r.C-O; 
        
        double a = 1;
        double b = 2* r.u.dot(diff);
        double c = (diff).squaredNorm() - R*R;
        double delta = b*b - 4 *a*c;

        if (delta < 0) return -1.0;
        double sqrtDelta = sqrt(delta);
        double t1 = (-b + sqrtDelta) /2; // / (2*a)
    
        if (t1 < 0) return -1.0;    
        t1 = (-b - sqrtDelta) /2; // closest point 
        
        if(interStruct.computeN){
            interStruct.P = r.computeIntersection(t1);
            interStruct.N = (interStruct.P-O) / R;
        } else 
        if(interStruct.computeP){
            interStruct.P = r.computeIntersection(t1);
        }

        return t1;
    }


    Vector3d getNormalAt(const Vector3d& P, bool normalize=true) const {
        Vector3d n = (P-O);    // normal direction
        if (normalize) n /= R; // P being on the surface, we know that norm(P-O) = R
        return n  ; 
    }
};

#endif // SPHERE_H
