#ifndef SCENE_H
#define SCENE_H

#include "Sphere.h"
#include "Light.h"

#include <Eigen/Dense>

using Eigen::Vector3d;

class Scene {
public:
	Light light = Light();
    Sphere* spheres[1] ; //= new Sphere[1];

    Scene() {
    	spheres[0] = new Sphere(Vector3d(0.0, 0.0, 0.0), 10.0, &light);
    }

    void intersect(Ray ray){

    	// for(i = 0; i < spheres.length; i++){

    	// }

    }
};

#endif // SCENE_H