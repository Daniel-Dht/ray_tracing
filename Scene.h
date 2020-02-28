#ifndef SCENE_H
#define SCENE_H

#include "Object.h"
#include "Sphere.h"
#include "Material.h"
#include "Ray.h"
#include "Triangle.h"
//#include "Geometry.h"
#include "Mesh.h"
#include "Vector.h"

#include <iostream>
#include <random>
#include <vector>
std::random_device rd;
std::default_random_engine engine(rd());
std::uniform_real_distribution<> uniform(0, 1);

#define M_PI 3.14159265359

typedef Vector3d Point3d;

Vector3d getRandomDir(const Vector3d& n){
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
	const int n_rays, max_bounce;
	int total_ray_casted = 0;

    std::vector<Object *> objects ;
    int n_objects=0;
    
	std::vector<Sphere *> lights;
    int n_lights=0;

	enum modes{NORMAL, MIRROR, TRANS};
    Scene(int n_rays, int max_bounce) : n_rays(n_rays), max_bounce(max_bounce) {

		Vector3f white = Vector3f(1,1,1);
		Vector3f red   = Vector3f(1,0,0);
		Vector3f green = Vector3f(0,1,0);
		Vector3f blue  = Vector3f(0,0,1);
		Vector3f fred   = Vector3f(255,69,0)/255;
		Vector3f fgreen = Vector3f(0,250,0)/255;
		Vector3f fblue  = Vector3f(24,116,205)/255;
		Vector3f fgray  = Vector3f(176,196,	222)/255;
		
		add_sphere(new Sphere(Vector3d(-12., 0., 30.), 3., Material(red, NORMAL , 0) ));
    	add_sphere(new Sphere(Vector3d(-12., 5., 0.), 8., Material(red, MIRROR, 0) ));
		add_sphere(new Sphere(Vector3d( 12., 0., 0.), 8., Material(red, MIRROR, 0) ));
		
    	// // walls :
    	add_sphere(new Sphere(Vector3d( 0., 0., -1000.), 970., Material(fgray, NORMAL, 0) ));    	    	
    	add_sphere(new Sphere(Vector3d( 0., 0.,  1000.), 900., Material(fgray, NORMAL, 0) )); 
    	
    	add_sphere(new Sphere(Vector3d( 0., -1000., 0.), 990., Material(fgray, NORMAL, 0) ));
    	add_sphere(new Sphere(Vector3d( 0.,  1000., 0.), 960., Material(fgray, NORMAL, 0) ));

    	add_sphere(new Sphere(Vector3d( 1000., 0., 0.), 970., Material(fblue, NORMAL, 0) ));
    	add_sphere(new Sphere(Vector3d(-1000., 0., 0.), 970., Material(fblue, NORMAL, 0) ));

		//add_triangle(new Triangle(Point3d(0.0099,  8.29565,  0.32215), Point3d(0.0134,  8.3601,  0.11505), Point3d(0.17135,  8.31755,  0.11635),  Material(red, NORMAL, 0) ));
		add_triangle(new Triangle(Point3d(-10., -15., 25.), Point3d(0., -10., 15.), Point3d(-10, 10, 8),  Material(red, MIRROR, 0) ));
		
		//add_triangle(new Triangle(Point3d(-5,  -5,  -5), Point3d(-5,  15,  -5), Point3d(15,  -5,  -5),  Material(red, NORMAL, 0) ));
		//add_triangle(new Triangle(Point3d(-2000,  -2000,  -5), Point3d(-2000,  5000,  -5), Point3d(5000,  -2000,  -5),  Material(red, NORMAL, 0) ));

		add_mesh(new Mesh("model.off", 20, Vector3d(0,-10,15),  Material(red, NORMAL, 0))) ;
		//add_mesh(new Mesh("cube.off", 2, Vector3d(15,10,20),  Material(red, NORMAL, 0))) ;
		//add_mesh(new Geometry("cube.obj", 5, Vector3d(0,0,0),  Material(red, NORMAL, 0))) ;

		// light
		//add_sphere(new Sphere(Vector3d(10.0, 30.0, 35.0), 10., Material(white, NORMAL, 1) ));
		add_light( new Sphere(Vector3d(10.0, 30.0, 35.0), 10., Material(white, NORMAL, 1) ));
		
    }	
	
	void add_mesh(Mesh *t){
		objects.push_back(t);
		n_objects++;
	}
	// void add_mesh(Geometry *t){
	// 	objects.push_back(t);
	// 	n_objects++;
	// }
	void add_triangle(Triangle *t){
		objects.push_back(t);
		n_objects++;
	}
	void add_sphere(Sphere *s){
		objects.push_back(s);
		n_objects++;
	}
	void add_light(Sphere *s){
		lights.push_back(s);
		n_lights++;
	}


	
    Vector3f intersect(const Ray& r, int count_bounce=0){		
		// find the closest sphere
		// if its a mirror : recast reflective ray
		// if its transparent : recast reflective and refracted rays
		// else:	
		// check if there are other sphere bewteen the light source and P: if yes -> P is in shadow
		// if no shadow: calculate the light intensity received on the intersection point P 
		//             + cast random rays for reflective light 
		
		Vector3f pixel_col = Vector3f(0,0,0);

		// max recursive call
		if (count_bounce==max_bounce) return pixel_col;

    	double min_t1 = 1E90;
    	int id_closest_object = -1;
    	int i = 0;
		InterStruct interStruct;
		InterStruct interStruct_min;
    	while(i < n_objects ) { // Scene::intersection du prof
    		    		
    		Object *s = objects[i];
    		double t1 = s->intersect(r, total_ray_casted, interStruct);

    		if ((t1>0) && (t1 < min_t1)){
	        	min_t1 = t1;
	        	id_closest_object = i;
				interStruct_min = interStruct;
	        }		
			i++;	        
    	}


    	if( id_closest_object == -1) {
			return pixel_col; // no intersection at all    					
    	}  
		
		Object *s1 = objects[id_closest_object]; // closest sphere
		if (s1->mat.emissivity > 0) return s1->mat.col * 255555.;

		Vector3d &P = interStruct_min.P;
		Vector3d &n = interStruct_min.N;
		// const Vector3d PP = r.computeIntersection(min_t1);	
		// Vector3d nn = s1->getNormalAt(P);
	
		if (n.dot(r.u) > 0) n*= -1.; // always show face even if normal is not well oriented

		if (s1->mat.isMirror) 
		{					 				
			return intersect(Ray(P+n*0.001,  r.u-n*2*n.dot(r.u)), count_bounce+1);
		}
		else 
		if (s1->mat.isTrans) 
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

		Sphere *light = dynamic_cast<Sphere*>(lights[0]);
		//Sphere *light = dynamic_cast<Sphere*>(objects[n_objects-1]);
		Vector3d l = light->O - P;  // direction from P to the light source 
		float scale_light = 1; // =cos(theta)

		if(1 && n_rays>1){ // ombres douce		
			// from the initial 	
			l*=-1;
			l.normalize();
			Vector3d random_dir = getRandomDir(l);
			scale_light = random_dir.dot(l); // = cos(tehta)		
			Vector3d new_source_light = light->O + random_dir * light->R;
			l = (new_source_light-P);
		}

		double d_square = l.squaredNorm(); // distance between P and light source
		l /= sqrt(d_square);
						
		// check if point P is in the shadow of another object
		bool shadow = false;
		InterStruct interShadow {false, false};

		for (int i = 0; i < n_objects; ++i) 
		{
			if (i==id_closest_object) continue;

			Object *s2 = objects[i];				
			double t2 = s2->intersect(Ray(P+n*0.001, l), total_ray_casted, interShadow); // ray partant de P, en direction de la source de lumière
		
			if (s2->mat.emissivity>0) continue; // if this the light..

			// si il existe un nouveau point d'intersection plus proche que la lumière n'est proche
			if ((t2>0) && (t2*t2 < d_square)) { // comparaisons au carrées pour éviter d'avoir a calculer une racine
				shadow = true; 
				break;
			}
		}

		// direct light    	 
		if (!shadow) {        
	        double intensity = std::max(0.0,l.dot(n)) * light->mat.emissivity / (d_square);        
	        pixel_col = s1->mat.col  * intensity * scale_light / M_PI ;
	    } else 
	 	if (n_rays==1 || 1) { 	
			double intensity = std::max(0.0,l.dot(n)) * light->mat.emissivity / (d_square) / 3;        
	        pixel_col = s1->mat.col  * intensity * scale_light / M_PI ;
		}

        // contribution indirecte
	    if (n_rays>1) { 		
			Vector3d random_dir = getRandomDir(n);
			Ray random_ray(P + n*0.001, random_dir);				
			pixel_col += intersect(random_ray, count_bounce+1) * s1->mat.col;	
	 	} 
    	
    	return pixel_col;    	
    }
};


#endif // SCENE_H