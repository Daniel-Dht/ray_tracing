#ifndef SCENE_H
#define SCENE_H

#include "Sphere.h"
#include "Light.h"
#include "Ray.h"

#include <iostream>
#include <Eigen/Dense>

using Eigen::Vector3d;
using Eigen::Vector3f;

class Scene {
public:
		
    std::vector<Sphere> spheres;

    Vector3f mirror = Vector3f(-1,-1,-1);
    Vector3f trans  = Vector3f(-2,-2,-2);

    Vector3f white = Vector3f(1,1,1);
    Vector3f red   = Vector3f(1,0,0);
    Vector3f green = Vector3f(0,1,0);
    Vector3f blue  = Vector3f(0,0,1);
    Vector3f fred   = Vector3f(255,69,0)/255;
    Vector3f fgreen = Vector3f(0,250,0)/255;
    Vector3f fblue  = Vector3f(24,116,205)/255;
	Vector3f fgray  = Vector3f(176,196,	222)/255;

    Light light = Light(Vector3d(10.0, 30.0, 35.0), 150000000.0);
    Scene() {

    	spheres.push_back( Sphere(Vector3d(-12.0, 0.0, 5.0), 8.0, trans));    	
		spheres.push_back( Sphere(Vector3d( 12.0, 0.0, 0.0), 8.0, mirror));

    	// walls :
    	spheres.push_back( Sphere(Vector3d( 0.0, 0.0, -1000.0), 970.0, mirror));    	    	
    	spheres.push_back( Sphere(Vector3d( 0.0, 0.0,  1000.0), 900.0, fgray)); 
    	
    	spheres.push_back( Sphere(Vector3d( 0.0, -1000.0, 0.0), 990.0, fgray));
    	spheres.push_back( Sphere(Vector3d( 0.0,  1000.0, 0.0), 960.0, fgray));

    	spheres.push_back( Sphere(Vector3d( 1000.0, 0.0, 0.0), 970.0, fblue));
    	spheres.push_back( Sphere(Vector3d(-1000.0, 0.0, 0.0), 970.0, fblue));
    }	

    Vector3f intersect(const Ray& r, int count_bounce=0){		
		
		Vector3f pixel_col = Vector3f(0,0,0);

		// recursive unction, n call max=5 
		if (count_bounce==10) return pixel_col;//Vector3f(-1,-1,-1);


    	double min_t1 = 1E99;
    	int min_i = -1;
    	int i = -1;

    	while(i < int(spheres.size()) ) { // Scne::intersection du prof
    		
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
				Vector3d n = (P-s1.O);          // normal of the object at P				
				n.normalize();
				pixel_col = intersect(Ray(P+0.001*n,  r.u-2*n.dot(r.u)*n), count_bounce+1);
				return pixel_col;
			}
			else 
			if (s1.isTrans) {				
				Vector3d n = (P-s1.O);          // normal of the object at P					
				n.normalize();

				double dotprod = r.u.dot(n); //
				double n_ratio = 1.0/1.5;
				if(dotprod > 0.0){
					n_ratio = 1.5;
					//dotprod *=-1.0;
					//n*=-1.0;
				}
				double radical = 1.0 - n_ratio*n_ratio * (1.0-dotprod*dotprod);
				
				if (radical > 0.0) {
					Vector3d u = n_ratio * (r.u - n*dotprod) - n*sqrt(radical); // le - du prof devient plus car on ne modifie pas n.
					//return intersect(Ray(P-0.001*n,  u), count_bounce+1); // -0.001*n car on part derrière la surface de la sphère
					return 0.9*intersect(Ray(P-0.001*n,  u), count_bounce+1) +0.1*intersect(Ray(P+0.001*n,  r.u-2*dotprod*n), count_bounce+1);
				}				
				
			}

    		// intersection point
    		Vector3d l = (light.source-P);  // direction from P to the light source
    		double d_square = l.squaredNorm();            // distance between P and light source
			l.normalize();

    		// check if object is in the shadow of another object
    		bool shadow = false;
    		for (i = 0; i < int(spheres.size()); ++i) // encore Scene::intersection du prof mais en évitant la sphère courante
    		{
    			if (i==min_i) continue;

    			Sphere s2 = spheres[i];
				//Vector3d n = (P-s1.O); // ces deux lignes doivent bouger au dessus si decommentées          
		        //n.normalize();	
    			double t2 = s2.intersect(Ray(P, l)); // rajouter  +0.001*n à P si ça bug
    		
    			// si le point d'intersetion est plus proche que la lumière n'est proche
    			if ((t2>0) && (t2*t2 < d_square)) { // au carré pour éviter d'avoir a calculer une racine
    				shadow = true; 
    				break;
    			}
    		}

    		// light    	 
    		if(!shadow) {
    			
		        Vector3d n = (P-s1.O);          // normal of the object at P	     
		        n.normalize();		        

		        double scale = l.dot(n) * light.power / (d_square); 
		        scale = std::clamp(std::pow(scale, 0.45), 0.0, 255.0);
		        //pixel_col += s1.col * scale;
		        pixel_col = s1.col * scale;
		    }
    	}
    
    	return pixel_col;
    	//return Vector3f(0.1,0.1,0.8);
    }
};

#endif // SCENE_H