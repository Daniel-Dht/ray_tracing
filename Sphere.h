#ifndef SPHERE_H
#define SPHERE_H

#include "Ray.h"
#include "Light.h"

#include "Vector.h"
#include "iostream" 
#include <limits>


using namespace std;

class Sphere {
public:
    Vector3d O;
    Vector3f col;
    double R;     
    bool isMirror=false;   
    bool isTrans=false;
    double emissivity=0;

    Sphere(const Vector3d O, double R, Vector3f col, int mode=0, double emissi=0) : O(O), R(R), col(col), emissivity(emissi) {     
        if     (mode==1)  isMirror = true;
        else if(mode==2)  isTrans = true;        
    }; 
    
    double intersect(const Ray& r, int& total_ray_casted){
        
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
        t1 = (-b - sqrtDelta) /2; // closer point 
        
        return t1;
    }
};

#endif // SPHERE_H
