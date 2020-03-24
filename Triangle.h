#ifndef TRIANGLE_H
#define TRIANGLE_H


#include <string>
#include <map>

#include "Ray.h"
#include "Material.h"

#include "Object.h"
#include "Vector.h"
#include "iostream" 
#include <limits>


using namespace std;
typedef Vector3d Point3d;

class Triangle : public Object { 
public:
    Point3d A,B,C;
	Vector3d N;  

    Triangle(Point3d A, Point3d B, Point3d C, Material mat, bool compute_normal=true, Vector3d n=Vector3d(0)) : A(A), B(B), C(C), Object(mat) { 
		if (compute_normal) N = (B - A).cross(C - A).getNormalized();
		else N=n;
	}; 

	Triangle(const Point3d &A, const Point3d &B, const Point3d &C, const Vector3d &n) : A(A), B(B), C(C), N(n) { 

	}; 

    
    double intersect(const Ray& r, int& total_ray_casted, InterStruct &interStruct) const {

		total_ray_casted++;
		
		// check si la normal du plan est collinéaire avec le rayon
		const double dot_N_ray = r.u.dot(N); 
        if(fabs(dot_N_ray) < 0.1 ) return -1.;

		// intersection plan/droite
		double t = (A - r.C).dot(N) / dot_N_ray;
		if (t < 0) return -1.;
		
		// intersection triangle/droite
		interStruct.P = r.computeIntersection(t);
		double  alpha, beta;

		const double y1 = (interStruct.P-A).dot(B-A); // y1
		const double  a = (B-A).dot(B-A); //  a = AB.AB
		const double  b = (B-A).dot(C-A); //  b = AC.AB
		const double y2 = (interStruct.P-A).dot(C-A); // y2
		const double  d = (C-A).dot(C-A); //  d = AC.AC
		
		beta   = (y1*b - y2*a) / (b*b - d*a);		
		alpha  = (y1-beta*b) / a;	
		// computeAlphaBeta(alpha, beta, interStruct.P);
					
		if (alpha > 0 && beta > 0 && (alpha + beta) < 1) {
			if(interStruct.computeN){				
				interStruct.N = N;
			} 
			return t;
		}
		return -1;
    }
	Vector3d getNormalAt(const Vector3d& P, bool normalize=true) const { 
        return N;
    }

	void computeAlphaBeta(double &alpha, double &beta, const Vector3d &P){
		// En partant de AP = alpha*AB+ beta*AC, on a:
		// y1 = AP.AB = alpha*AB.AB + beta*AC.AB
		// y2 = AP.AC = alpha*AB.AC + beta*AC.AC
		// Résolution du système : 
		// y1 = alpha * a + beta * b
		// y2 = alpha * b + beta * d
		// => beta  = (y1/a - y2/b) / (b/a - d/b)
		// => beta  = (y1.b - y2.a) / (b.b - d.a)
		// et alpha = (y1-beta*b) / a

		const double y1 = (P-A).dot(B-A); // y1
		const double  a = (B-A).dot(B-A); //  a = AB.AB
		const double  b = (B-A).dot(C-A); //  b = AC.AB
		const double y2 = (P-A).dot(C-A); // y2
		const double  d = (C-A).dot(C-A); //  d = AC.AC
		
		beta   = (y1*b - y2*a) / (b*b - d*a);		
		alpha  = (y1-beta*b) / a;	
	}

	double intersect_interpolateN(const Ray& r, int& total_ray_casted, InterStruct &interStruct, const Vector3d &n1,  const Vector3d &n2,  const Vector3d &n3, const float &min_t) const {

		total_ray_casted++;
		
		// check si la normal du plan est collinéaire avec le rayon
		const double dot_N_ray = r.u.dot(N); 
        if(fabs(dot_N_ray) < 0.01 ) return -1.;

		// intersection plan/droite
		const double t = (A - r.C).dot(N) / dot_N_ray;
		if (t < 0 || t>min_t) return -1.;
		
		// intersection triangle/droite
		interStruct.P = r.computeIntersection(t);
		double  alpha, beta, gamma;

		const double y1 = (interStruct.P-A).dot(B-A); // y1
		const double  a = (B-A).dot(B-A); //  a = AB.AB
		const double  b = (B-A).dot(C-A); //  b = AC.AB
		const double y2 = (interStruct.P-A).dot(C-A); // y2
		const double  d = (C-A).dot(C-A); //  d = AC.AC
		
		beta   = (y1*b - y2*a) / (b*b - d*a);		
		alpha  = (y1-beta*b) / a;	
		
		// computeAlphaBeta(alpha, beta, interStruct.P);

		if (alpha > 0 && beta > 0 && (alpha + beta) < 1) {
			if(interStruct.computeN){	
				gamma  = 1 - alpha - beta;			
				interStruct.N = n1*alpha + n2*beta + n3*gamma ;
				interStruct.N.normalize();
			} 
			return t;
		}
		return -1;
    }
};

#endif //TRIANGLE_H