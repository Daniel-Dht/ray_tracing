
class Sphere {
 PVector center;
 float r;
 
 Sphere(PVector _center, float _r) {
   center = _center;
   r = _r;
 }
 
 boolean intersect(Ray ray){
   float a = 1;
   float b = 2* ray.dir.dot(ray.origin.sub(center));
   float c = (ray.origin.sub(center)).magSq() - r*r;
   float delta = b*b - 4*a*c;
   
   if (delta < 0) return false;
   double t1 = (-b + sqrt(delta)) / (2*a);
   if (t1 < 0) return false;
   return true;
 }
}
