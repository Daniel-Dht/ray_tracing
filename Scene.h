#ifndef SCENE_H
#define SCENE_H

#include "Sphere.h"
#include "Light.h"
#include "Ray.h"

#include <iostream>
#include "Vector.h"

#include <random>
std::random_device rd;
std::default_random_engine engine(rd());
std::uniform_real_distribution<> uniform(0, 1);

#define M_PI 3.14159265359

Vector3d getRandomDir(Vector3d n){
	double r1 = uniform(engine);
	double r2 = uniform(engine);
	Vector3d local_random_dir(cos(2*M_PI*r1)*sqrt(1 - r2), sin(2*M_PI*r1)*sqrt(1 - r2), sqrt(r2));
	//Vector3d rand_vec(uniform(engine)-0.5, uniform(engine)-0.5, uniform(engine)); 
	Vector3d rand_vec= Vector3d(1,0,0);
	Vector3d tangent1 = n.cross(rand_vec).getNormalized();
	Vector3d tangent2 = tangent1.cross(n).getNormalized();
	return n*local_random_dir[2] + tangent1*local_random_dir[0] + tangent2*local_random_dir[1];
}

class Scene {
public:
	int n_rays, max_bounce;
    std::vector<Sphere> spheres;
    int n_spheres;
    int total_ray_casted = 0;

    Vector3f white = Vector3f(1,1,1);
    Vector3f red   = Vector3f(1,0,0);
    Vector3f green = Vector3f(0,1,0);
    Vector3f blue  = Vector3f(0,0,1);
    Vector3f fred   = Vector3f(255,69,0)/255;
    Vector3f fgreen = Vector3f(0,250,0)/255;
    Vector3f fblue  = Vector3f(24,116,205)/255;
	Vector3f fgray  = Vector3f(176,196,	222)/255;

    Light light = Light(Vector3d(10.0, 30.0, 35.0), 500000000.0, 5);

	// NORMAL : mode = 0
	// MIRROR : mode = 1
	// TRANS  : mode = 2
	enum modes{NORMAL, MIRROR, TRANS};
    Scene(int n_rays, int max_bounce) : n_rays(n_rays), max_bounce(max_bounce) {

		spheres.push_back( Sphere(Vector3d(-12., 0., 30.), 3., red, NORMAL , 0));
    	spheres.push_back( Sphere(Vector3d(-12., 5., 0.), 8., red, MIRROR, 0));
		spheres.push_back( Sphere(Vector3d( 12., 0., 0.), 8., red, MIRROR, 0));
		
    	// walls :
    	spheres.push_back( Sphere(Vector3d( 0., 0., -1000.), 970., fgray, NORMAL, 0));    	    	
    	spheres.push_back( Sphere(Vector3d( 0., 0.,  1000.), 900., fgray, NORMAL, 0)); 
    	
    	spheres.push_back( Sphere(Vector3d( 0., -1000., 0.), 990., fgray, NORMAL, 0));
    	spheres.push_back( Sphere(Vector3d( 0.,  1000., 0.), 960., fgray, NORMAL, 0));

    	spheres.push_back( Sphere(Vector3d( 1000., 0., 0.), 970., fblue, NORMAL, 0));
    	spheres.push_back( Sphere(Vector3d(-1000., 0., 0.), 970., fblue, NORMAL, 0));
    	n_spheres = int(spheres.size());
    }	

	// DIRECT   : mode = 0
	// INDIRECT : mode = 1
    Vector3f intersect(const Ray& r, int count_bounce=0){		
		
		Vector3f pixel_col = Vector3f(0,0,0);

		// max recursive call
		if (count_bounce==max_bounce) return pixel_col;

    	double min_t1 = 1E99;
    	int id_closer_sphere = -1;
    	int i = -1;

    	while(i < n_spheres ) { // Scene::intersection du prof
    		
    		i += 1;
    		Sphere &s = spheres[i];
    		double t1 = s.intersect(r, total_ray_casted);

    		if ((t1 > 0) && (t1 < min_t1)){
	        	min_t1 = t1;
	        	id_closer_sphere = i;
	        }
			        
    	}

    	if( id_closer_sphere == -1) {
			return pixel_col; // no intersection at all    		
			
    	}  
		
		Sphere &s1 = spheres[id_closer_sphere]; // closer sphere
		Vector3d P = r.C + r.u*min_t1;          // intersection point
		Vector3d n = (P-s1.O); n.normalize();   // normal

		if (s1.isMirror) 
		{					 				
			pixel_col = intersect(Ray(P+n*0.001,  r.u-n*2*n.dot(r.u)), count_bounce+1);
			return pixel_col;
		}
		else 
		if (s1.isTrans) 
		{				
			double dotprod = r.u.dot(n); //
			double n_ratio = 1.0/1.5;
			if(dotprod > 0.0) n_ratio = 1.5;
			
			double radical = 1.0 - n_ratio*n_ratio * (1.0-dotprod*dotprod);
			if (radical > 0.0) {
				Vector3d u = (r.u - n*dotprod)*n_ratio - n*sqrt(radical); // le - du prof devient plus car on ne modifie pas n.
				//return intersect(Ray(P-0.001*n,  u), count_bounce+1); // -0.001*n car on part derrière la surface de la sphère
				return intersect(Ray(P-n*0.001,  u), count_bounce+1)*0.9 + intersect(Ray(P+n*0.001,  r.u-n*2*dotprod), count_bounce+1)*0.1;
			}								
		}

		Vector3d l = (light.source-P);     // direction from P to the light source
		float scale_light = 1; // =cos(theta)

		if(1){ // ombres douce			
			l.normalize();
			Vector3d random_dir = getRandomDir(l);
			scale_light = random_dir.dot(l);			
			Vector3d new_source_light = light.source + random_dir * light.size;
			l = (new_source_light-P);
		} 

		double d_square = l.squaredNorm(); // distance between P and light source
		l /= sqrt(d_square);
						
	
		// check if point P is in the shadow of another object
		bool shadow = false;
		for (int i = 0; i < n_spheres; ++i) 
		{
			if (i==id_closer_sphere) continue;

			Sphere &s2 = spheres[i];				
			double t2 = s2.intersect(Ray(P+n*0.001, l), total_ray_casted); // ray partant de p, en direction de la lumière
		
			// si il existe un nouveau point d'intersection plus proche que la lumière n'est proche
			if ((t2>0) && (t2*t2 < d_square)) { // comparaisons au carrées pour éviter d'avoir a calculer une racine
				shadow = true; 
				break;
			}
		}

		// direct light    	 
		if(!shadow) {        
	        double intensity = std::max(0.0,l.dot(n)) * light.power / (d_square);        
	        pixel_col = s1.col  * intensity * scale_light / M_PI ;	        
	    }

        //ajout de la contribution indirecte
	    if(n_rays>1) { 
			
			Vector3d random_dir = getRandomDir(n);
			Ray random_ray(P + n*0.001, random_dir);				
			pixel_col += intersect(random_ray, count_bounce+1) * s1.col;	
	 	}
    	
    	return pixel_col;    	
    }
};


#endif // SCENE_H