#ifndef OBJECT_H
#define OBJECT_H

#include "Ray.h"
#include "Material.h"

#include "Vector.h"
#include "iostream" 

class Object
{
public:
	
    Material mat; 

	Object() {};	
	Object(Material mat) : mat(mat) { }; 
	virtual double intersect(const Ray& r, int& total_ray_casted, InterStruct &interStruct) const = 0;
	//virtual Vector3d getNormalAt(const Vector3d& P, bool normalize=true) const = 0;
};

#endif //OBJECT_H
