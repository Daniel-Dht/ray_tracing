#ifndef SCENE_H
#define SCENE_H

#include "Sphere.h"
#include "Light.h"
#include "Ray.h"

#include <limits>
//#include <iostream>
#include <Eigen/Dense>

using Eigen::Vector3d;
using Eigen::Vector3f;

class Scene {
public:
	size_t n = 2;
	Light light = Light();
    std::vector<Sphere> spheres;

    Vector3f white = Vector3f(1,1,1);
    Vector3f red   = Vector3f(1,0,0);
    Vector3f green = Vector3f(0,1,0);
    Vector3f blue  = Vector3f(0,0,1);

    Scene() {
    	spheres.push_back( Sphere(Vector3d(-15.0, 0.0, 0.0), 10.0, red, &light));
    	spheres.push_back( Sphere(Vector3d( 10.0, 0.0, -5.0), 10.0, blue,  &light));
    	spheres.push_back( Sphere(Vector3d( 0.0, -20.0, -10.0), 10.0, green, &light));

    	// walls :
    	spheres.push_back( Sphere(Vector3d( 0.0, 0.0, -1000.0), 970.0, white, &light));
    	    	
    	spheres.push_back( Sphere(Vector3d( 0.0, -1000.0, 0.0), 970.0, white, &light));
    	spheres.push_back( Sphere(Vector3d( 0.0,  1000.0, 0.0), 970.0, white, &light));

    	spheres.push_back( Sphere(Vector3d( 1000.0, 0.0, 0.0), 970.0, white, &light));
    	spheres.push_back( Sphere(Vector3d(-1000.0, 0.0, 0.0), 970.0, white, &light));
    }

    Vector3f intersect(Ray r){

    	double min_t1 = std::numeric_limits<double>::infinity();;
    	int min_i = -1;

    	// intersecetion function
    	// arg : i
    	int i = -1;
    	while(i < int(spheres.size()) ) {
    		i += 1;
    		Sphere s = spheres[i];

    		double a = 1;
	        double b = 2* r.u.dot(r.C-s.O);
	        double c = (r.C-s.O).squaredNorm() - s.R*s.R;
	        double delta = b*b - 4 *a*c;

	        if (delta < 0) continue;
	        double t1 = (-b + sqrt(delta)) /2; // / (2*a)
	        // si le point de la boule le plus proche de la caméra est quand même
	        // derrière elle, ce t1 négatif permet de savoir qu'elle l'est.
	        if (t1 < 0) continue;
	        // mais si il est positif l'autre t1 peut être quand même négatif:
	        // la caméra se trouve alors à l'interieur de la sphère.
	        t1 = (-b - sqrt(delta)) /2; // closer point 
	        //if (t1 <s.R) continue; 

	        if (t1 < min_t1){
	        	min_t1 = t1;
	        	min_i = i;
	        }	        
    	}
    	if( min_i != -1) { // intersection occured:
	    	Sphere s = spheres[min_i];
	        Vector3d P = r.C + min_t1*r.u;    // intersection point
	        Vector3d n = (P-s.O);             // normal of the object at P
	        Vector3d l = (s.light->source-P); // direction from P to the light source
	        double d = l.squaredNorm();       // distance between P and light source
	        n.normalize();
	        l.normalize();

	        double scale = l.dot(n) * s.light->power / (d); 

	        scale = std::clamp(scale, 0.0, 255.0);
	        //std::cout << final_col << std::endl;
	        return s.col * scale;
	        // return Vector3f(1,1,1)*scale;
    	}
    	
    	return Vector3f(-1,-1,-1);
    }
};

#endif // SCENE_H