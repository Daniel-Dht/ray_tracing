#ifndef SPHERE_H
#define SPHERE_H

#include "Ray.h"
#include "Light.h"

#include <algorithm> // for std::clamp (c++ 17)
#include <Eigen/Dense>

using Eigen::Vector3d;

class Sphere {
public:
    Vector3d O;
    double R;
    Light* light;
    Sphere(const Vector3d _O, double _R, Light* _light) {
        O = _O;
        R = _R;
        light = _light;
        //std::cout<< light->L <<std::endl;
    }; 
    
    //Sphere() {};

    float intersect(const Ray& r){

        double a = 1;
        double b = 2* r.u.dot(r.C-O);
        double c = (r.C-O).squaredNorm() - R*R;
        double delta = b*b - 4 *a*c;

        if (delta < 0) return -1.0;
        double t1 = (-b + sqrt(delta)) /2; // / (2*a)
        if (t1 < 0) return -1.0;
    
        t1 = (-b - sqrt(delta)) /2; // lower t
        // light:
        Vector3d P = r.C + t1*r.u;             // intersection point
        Vector3d n = (P-O);                    // normal of the object at P
        Vector3d l = (light->source -P);    // direction from P to the light source
        double d = l.squaredNorm();            // distance between P and light source
        n.normalize();
        l.normalize();

        float val = -l.dot(n) * light->power / (d); 
        // std::cout<< val<< std::endl;
        return std::clamp(val, 0.0f, 255.0f);
    }
};

#endif // SPHERE_H
