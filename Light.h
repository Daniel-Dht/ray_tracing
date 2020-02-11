#ifndef LIGHT_H
#define LIGHT_H

#include "Vector.h"

class Light {
public:
    const float size = 5; // pour les ombres douces
    const double power;    // light intensity
    const Vector3d source;
    Light(const Vector3d& src, double pow, float size) : power(pow), source(src), size(size) {};
    
};

#endif // LIGHT_H
