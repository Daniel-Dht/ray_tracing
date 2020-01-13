#ifndef SPHERE_H
#define SPHERE_H

#include "Ray.h"
#include "Light.h"

//  // for std::clamp (c++ 17)
#include <Eigen/Dense>
#include "iostream" 
#include <limits>

using Eigen::Vector3d;
using Eigen::Vector3f;

class Sphere {
public:
    Vector3d O;
    Vector3f col;
    double R; 
    //double n_ratio = 0.76923; // for transparency material 1/1.3
    bool isMirror=false;
    bool isTrans=false;

    Sphere() { };
    Sphere(const Vector3d _O, double _R, Vector3f _col) {
        O = _O;
        R = _R;
        col = _col;        
        if(col[0]==-1)
            isMirror = true;
        if(col[0]==-2)
            isTrans = true;        
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
        //if (t1 <R) return -1.0;

        return t1;
    }
};

#endif // SPHERE_H
