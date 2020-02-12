#ifndef MATERIAL_H
#define MATERIAL_H

#include "Vector.h"

class Material {
public:
    Vector3f col;
    bool isMirror=false;   
    bool isTrans=false;
    double emissivity=500000000.0;

    Material() {};
    Material(Vector3f col_, int mode=0, double emissi=0) {     
        if     (mode==1)  isMirror = true;
        else if(mode==2)  isTrans = true;        
        col = col_;
        emissivity *= emissi;
    }; 

};

#endif // MATERIAL_H
