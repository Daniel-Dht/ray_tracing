#ifndef VECTOR_H
#define VECTOR_H

#include<cmath>
#include<iostream>

using namespace std;

class Vector3d
{
public:    
    double x;
    double y;
    double z;

    Vector3d(double x_, double y_, double z_): x(x_), y(y_), z(z_) {}
    
    double dot(const Vector3d &p) const { return x*p.x + y*p.y + z*p.z;}
    double norm() const { return sqrt(x*x + y*y + z*z) ;}
    double squaredNorm() const { return x*x + y*y + z*z ;}
    void normalize() {
        double a = sqrt(x*x + y*y + z*z);
        x/=a; y/=a; z/=a;
    }
    
    void operator*=(double a)  { x*=a; y*=a; z*=a;}
    void operator/=(double a)  { x/=a; y/=a; z/=a;}
    Vector3d operator*(double a) const { return Vector3d(x*a, y*a, z*a);}    
    Vector3d operator/(double a) const { return Vector3d(x/a, y/a, z/a);}
    Vector3d operator+(const Vector3d &p) const { return Vector3d(x+p.x, y+p.y, z+p.z);}
    Vector3d operator-(const Vector3d &p) const { return Vector3d(x-p.x, y-p.y, z-p.z);}      
                                                    //_y*p._z - _z*p._y,       _z*p._x - _x*p._z,       _x*p._y - _y*p._x 

    Vector3d cross(const Vector3d &p) const {return Vector3d(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x );};
    Vector3d getNormalized() const {double a = sqrt(x*x + y*y + z*z); return Vector3d(x/a, y/a, z/a);}
    double operator[](int i) const {
        if(i==0) return x;
        else if(i==1) return y;
        else return z;
    }
    friend ostream& operator<<(ostream& os, const Vector3d& dt);
};

class Vector3f
{
public:    
    float x;
    float y;
    float z;

    Vector3f() {};
    Vector3f(float x_, float y_, float z_): x(x_), y(y_), z(z_) {}
    
    Vector3f operator*(float a) const { return Vector3f(x*a, y*a, z*a);}
    Vector3f operator/(float a) const { return Vector3f(x/a, y/a, z/a);}
    Vector3f operator+(const Vector3f &p) const { return Vector3f(x+p.x, y+p.y, z+p.z);}
    Vector3f operator-(const Vector3f &p) const { return Vector3f(x-p.x, y-p.y, z-p.z);}   
    void operator*=(float a)  { x*=a; y*=a; z*=a;}
    void operator/=(float a)  { x/=a; y/=a; z/=a;}
    void operator+=(const Vector3f &p) { x+=p.x; y+=p.y; z+=p.z;} 

    float operator[](int i) const {
        if(i==0) return x;
        else if(i==1) return y;
        else return z;
    }
    friend ostream& operator<<(ostream& os, const Vector3f& dt);
};

ostream& operator<<(ostream& os, const Vector3d& V)
{
    os << V.x << ",  " << V.y << ",  " << V.z;
    return os;
}
ostream& operator<<(ostream& os, const Vector3f& V)
{
    os << V.x << ",  " << V.y << ",  " << V.z;
    return os;
}

Vector3f operator*(const Vector3f &v1, const Vector3f &v2) {
    return Vector3f(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}


#endif // VECTOR_H
