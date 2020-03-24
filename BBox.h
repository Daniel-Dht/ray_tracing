#ifndef BBOX_H
#define BBOX_H

#include "Vector.h"
#include "Ray.h"

#include <vector>

class BBox {
public:

    Vector3d c_min, c_max; // box corners
    
    BBox() {};
    BBox(Vector3d c_min, Vector3d c_max) : c_min(c_min), c_max(c_max) { };

    bool intersect(const Ray& ray) const
    {
        double t_1, t_2, t_min, t_max;
        t_1 = (c_min[0] - ray.C[0]) / ray.u[0];
        t_2 = (c_max[0] - ray.C[0]) / ray.u[0];
        t_min = fmin(t_1, t_2);
        t_max = fmax(t_1, t_2);
        for (int k = 1; k < 3; ++k)
        {
            t_1 = (c_min[k] - ray.C[k]) / ray.u[k];
            t_2 = (c_max[k] - ray.C[k]) / ray.u[k];
            t_min = fmax(t_min, fmin(t_1, t_2));
            t_max = fmin(t_max, fmax(t_1, t_2));
        }
        return t_min > 0 && t_min < t_max;
    }
};

class BVH {
public:
    BBox bbox;
    int i0, i1;

    // ptr to the 2 childs in the tree
    BVH *bl = NULL;
    BVH *br = NULL; 
    
};






#endif // BBOX_H
