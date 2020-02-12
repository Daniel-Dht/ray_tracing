
#include "Ray.h"
#include "Material.h"

#include "Vector.h"
#include "iostream" 

class Object
{
public:

	Vector3d O;
    double R;    
    Material mat; 

	Object() {};
	virtual double intersect(const Ray& r, int& total_ray_casted) = 0;
};

