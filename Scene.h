#ifndef SCENE_H
#define SCENE_H

#include "Sphere.h"
#include "Light.h"
#include "Ray.h"

#include <limits>
#include <iostream>
#include <Eigen/Dense>

using Eigen::Vector3d;
using Eigen::Vector3f;

class Scene {
public:
		
    std::vector<Sphere> spheres;

    Vector3f white = Vector3f(1,1,1);
    Vector3f red   = Vector3f(1,0,0);
    Vector3f green = Vector3f(0,1,0);
    Vector3f blue  = Vector3f(0,0,1);
    Vector3f fred   = Vector3f(255,69,0)/255;
    Vector3f fgreen = Vector3f(0,250,0)/255;
    Vector3f fblue  = Vector3f(24,116,205)/255;
	Vector3f fgray  = Vector3f(176,196,	222)/255;

    Light light = Light(Vector3d(10.0, 30.0, 35.0), 55000000.0);
    Scene() {

    	spheres.push_back( Sphere(Vector3d(-12.0, 0.0, 0.0), 8.0, red, true));
		spheres.push_back( Sphere(Vector3d( 12.0, 0.0, 0.0), 8.0, red, true));

    	// walls :
    	spheres.push_back( Sphere(Vector3d( 0.0, 0.0, -1000.0), 970.0, fgreen));    	    	
    	//spheres.push_back( Sphere(Vector3d( 0.0, 0.0,  1000.0), 970.0, fgreen)); 
    	
    	spheres.push_back( Sphere(Vector3d( 0.0, -1000.0, 0.0), 990.0, fgray));
    	spheres.push_back( Sphere(Vector3d( 0.0,  1000.0, 0.0), 960.0, fgray));

    	spheres.push_back( Sphere(Vector3d( 1000.0, 0.0, 0.0), 970.0, fblue));
    	spheres.push_back( Sphere(Vector3d(-1000.0, 0.0, 0.0), 970.0, fblue));
    }	

	
    Vector3f intersect(Ray r, int count_bounce=0){		
		
		// recursive unction, n call max=5 
		if (count_bounce==10) return Vector3f(-1,-1,-1);

    	double min_t1 = std::numeric_limits<double>::infinity();
    	int min_i = -1;
    	int i = -1;

    	while(i < int(spheres.size()) ) {
    		
    		i += 1;
    		Sphere s = spheres[i];
    		double t1 = s.intersect(r);

    		if ((t1 > 0) && (t1 < min_t1)){
	        	min_t1 = t1;
	        	min_i = i;
	        }	        
    	}

    	if( min_i != -1) { // intersection occured:

			Sphere s1 = spheres[min_i];
			Vector3d P = r.C + min_t1*r.u;

			if (s1.isMirror) {
				Vector3d l = (light.source-P);  // direction from P to the light source
				Vector3d n = (P-s1.O);          // normal of the object at P	
				l.normalize();
				n.normalize();
				return intersect(Ray(P+0.001*n,  -l+2*n.dot(l)*n), count_bounce+1);
			}

    		// intersection point
    		Vector3d l = (light.source-P);  // direction from P to the light source
    		double d = l.norm();            // distance between P and light source
			l.normalize();

    		// check if object is in the shadow of another object
    		bool shadow = false;
    		for (i = 0; i < int(spheres.size()); ++i)
    		{
    			if (i==min_i) continue;

    			Sphere s2 = spheres[i];
				Vector3d n = (P-s1.O);           
		        n.normalize();	
    			double t2 = s2.intersect(Ray(P+0.001*n, l));
    			
    			if ((t2>0) && (t2<d)) {
    				shadow = true; 
    				//std::cout << d << "  "<< t2<< std::endl;
    				break;
    			}
    		}

    		// light    	 
    		if(!shadow) {
    			
		        Vector3d n = (P-s1.O);          // normal of the object at P	     
		        n.normalize();		        

		        double scale = l.dot(n) * light.power / (d*d); 
		        scale = std::clamp(std::pow(scale, 0.45), 0.0, 255.0);
		        //std::cout << final_col << std::endl;
		        return s1.col * scale;
		        // return Vector3f(1,1,1)*scale;
		    }
    	}
    
    	return Vector3f(-1,-1,-1);
    }
};

#endif // SCENE_H