#ifndef SPHERE_H
#define SPHERE_H

#include "Ray.h"
#include "Light.h"

#include <algorithm> // for std::clamp (c++ 17)
#include <Eigen/Dense>
#include "iostream" 

using Eigen::Vector3d;
using Eigen::Vector3f;

class Sphere {
public:
    Vector3d O;
    Vector3f col;
    double R;
    Light* light;
    Sphere() { };
    Sphere(const Vector3d _O, double _R, Vector3f _col,Light* _light) {
        O = _O;
        R = _R;
        col = _col;
        light = _light;
        //std::cout<< light->L <<std::endl;
    }; 
    
    

    double intersect(const Ray& r){

        double a = 1;
        double b = 2* r.u.dot(r.C-O);
        double c = (r.C-O).squaredNorm() - R*R;
        double delta = b*b - 4 *a*c;

        if (delta < 0) return -1.0;
        double t1 = (-b + sqrt(delta)) /2; // / (2*a)
        // si le point de la boule le plus proche de la caméra est quand même
        // derrière elle, ce t1 négatif permet de savoir qu'elle l'est.
        if (t1 < 0) return -1.0;
        // mais si il est positif l'autre t1 peut être quand même négatif:
        // la caméra se trouve à l'interieur de la sphère.
        t1 = (-b - sqrt(delta)) /2; // closer point 
        if (t1 <R) return -1.0;

        // light:
        Vector3d P = r.C + t1*r.u;        // intersection point
        Vector3d n = (P-O);               // normal of the object at P
        Vector3d l = (light->source-P);  // direction from P to the light source
        double d = l.squaredNorm();       // distance between P and light source
        n.normalize();
        l.normalize();

        double val = l.dot(n) * light->power / (d); 

        //std::cout << val << std::endl;
        return std::clamp(val, 0.0, 255.0);
    }
};

#endif // SPHERE_H
