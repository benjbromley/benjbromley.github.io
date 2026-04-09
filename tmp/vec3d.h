#include <iostream>
#include <math.h>
#include <stdlib.h>

#ifndef _vec3d_h
#define _vec3d_h 1


class vec3d {

  public:

    double x,y,z;

    // constructor; the : *,*,* {...} initializes the required args.
    vec3d(const double &X, const double &Y, const double &Z)
    : x(X), y(Y), z(Z) {}

    // zero to start stupid, but saves a few headaches....
    vec3d() { x=y=z=0.0;}


    // use either .x,.y,.z or [0],[1],[2]
    // first, we need a reference to each element. don't bother w/checks...
    double& vec3del(const vec3d& v, const int i) const {
        return (double&)(*(&v.x+i));
    }
    // now overload the [] operator...
    double& operator[](const int i) {
	return vec3del(*this,i);
    }
    double operator[](const int i) const {
	return vec3del(*this,i);
    }

    // assignment to a double array
    vec3d& operator=(const double* V) 
    { x = V[0]; y = V[1]; z = V[2]; return *this; }
    
    
    // vector arithmetic

    // vec add-assign:
    vec3d& operator+=(const vec3d &v)
    { x += v.x; y += v.y; z += v.z; return *this; }

    // mult by scalar, assign:
    vec3d& operator*=(const double m) 
    { x *= m; y *= m; z *= m; return *this; }

    // vec subract-assign:
    vec3d& operator-=(const vec3d &v) 
    {  x -= v.x; y -= v.y, z -= v.z; return *this; }

    // div by scalar, assign:
    vec3d& operator/=(const double d) 
    { double id=1.0/d; return *this *= id; }

    // norm
    double norm() const { return sqrt(x*x+y*y+z*z); }

    double norm2() const { return (x*x+y*y+z*z); }

    // cute junk
    vec3d unit() const { 
	double n = 1.0 / norm();
	return vec3d(x*n,y*n,z*n);
    }

    bool is_null() { return (!x && !y && !z); }
};


// vector add:
inline vec3d operator+(const vec3d& u, const vec3d &v)
{ vec3d s = u; return s+=v; }

//  vector subtract:
inline vec3d operator-(const vec3d& u, const vec3d &v)
{ vec3d s = u; return s-=v; }

// vec * double:
inline vec3d operator*(const vec3d &v, const double m)
{ vec3d p = v; return p*=m; }
 
// double * vec:
inline vec3d operator*(const double m, const vec3d &v)
{ vec3d p = v; return p*=m; }

// vec / double:
inline vec3d operator/(const vec3d& v, const double d)
{ vec3d q = v; return q /= d; }


// unary minus/plus, to allow vec3d v2 = -v1:
inline vec3d operator-(const vec3d& v) { return v*-1.0; }
inline vec3d operator+(const vec3d& v) { return v; }

// comparison operators, ==, !=, and ! :

inline bool operator==(const vec3d& u, const vec3d& v)
{ return (u.x==v.x && u.y==v.y && u.z==v.z); }
inline bool operator!=(const vec3d& u, const vec3d& v)
{ return !(v==u); }


// A const vec3d version of NULL:
const vec3d NULLvec3d(0,0,0);


// pure vector operations:

// cross product (new)
inline vec3d operator%(const vec3d& a, const vec3d& b)
{ return vec3d(a.y*b.z-b.y*a.z,a.z*b.x-b.z*a.x,a.x*b.y-b.x*a.y); }

// cross product (old)
inline vec3d operator^(const vec3d& a, const vec3d& b)
{ return vec3d(a.y*b.z-b.y*a.z,a.z*b.x-b.z*a.x,a.x*b.y-b.x*a.y); }

// dot product (old)
inline double operator*(const vec3d& u, const vec3d &v)
{ return (u.x*v.x+u.y*v.y+u.z*v.z); }

// for printing out a vector:

inline std::ostream& operator<< (std::ostream& s, const vec3d& a) {
    return s << "(" << a.x << "," << a.y << "," << a.z << ")";
}

// exchange values of two vectors
inline void swap(vec3d& a, vec3d& b) { 
    vec3d c;
    c = a;
    a = b;
    b = c;
}

// seed rand() from cstdlib
inline void srandvec3d(const int seed) { srand(seed); }

// generate a unit vector with random orientation
inline vec3d randunitvec3d() {
        const double twopi = 6.2831853071795864769252867665590;
	const double imax = 1.0/RAND_MAX;
        double phi = twopi * rand() * imax;
        double theta = acos(1.-2.*rand() * imax);
        return vec3d(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));
}

// generate a random vector inside a unit sphere
inline vec3d randvec3d() {
    const double imax = 1.0/RAND_MAX;
    double r = pow(rand()*imax,1.0/3.0);
    return randunitvec3d() * r;
}


class arr3x3 
{
 public:
    vec3d x,y,z;

    arr3x3() {}
    arr3x3(const vec3d& r1, const vec3d& r2, const vec3d& r3) :
	x(r1), y(r2), z(r3) {}
    // set all elements at once
    arr3x3(const double q) :
	x(vec3d(q,q,q)), y(vec3d(q,q,q)), z(vec3d(q,q,q)) {}
    // set diag
    arr3x3(const double a, const double b, const double c) :
	x(vec3d(a,0,0)), y(vec3d(0,b,0)), z(vec3d(0,0,c)) {}

    // use either .x,.y,.z or [0],[1],[2]
    // first, we need a reference to each element. don't bother w/checks...
    vec3d& arr3row(arr3x3& a, int i) {
        return (vec3d&)(*(&a.x+i));
    }
    // now overload the [] operator...
    vec3d& operator[](int i) {
        return arr3row(*this,i);
    }

    arr3x3 operator=(const arr3x3 &a)
    {
	x = a.x; y = a.y; z = a.z;
        return *this;
    }

    // basic arithemtic
    arr3x3 operator+=(const arr3x3 &a)
    {
        x += a.x;
        y += a.y;
        z += a.z;
        return *this;
    }

    arr3x3 operator-=(const arr3x3 &a)
    {
        x -= a.x;
        y -= a.y;
        z -= a.z;
        return *this;
    }

    arr3x3 operator+(const arr3x3 &a) const { 
	return arr3x3(x+a.x,y+a.y,z+a.z); }

    arr3x3 operator-(const arr3x3 &a) const { 
	return arr3x3(x-a.x,y-a.y,z-a.z); }

    vec3d operator*(const vec3d &v) const { 
	return vec3d(x*v,y*v,z*v); }

    arr3x3 operator*(const arr3x3 &a) const { 
        return arr3x3((*this)*a.x,(*this)*a.y,(*this)*a.z);
    }

    // arithmetic with scalars: * -and- / double
    arr3x3 operator*(const double q) const { 
	return arr3x3(x*q,y*q,z*q); }
    arr3x3 operator/(const double q) const { 
	return arr3x3(x/q,y/q,z/q); }
};


// transpose 
inline arr3x3 operator~(const arr3x3 &a) 
{
    return arr3x3(vec3d(a.x.x,a.y.x,a.z.x),
                  vec3d(a.x.y,a.y.y,a.z.y),
                  vec3d(a.x.z,a.y.z,a.z.z));
}


// outer vector product
inline arr3x3 operator&(const vec3d& a, const vec3d& b) {
    return arr3x3(a.x*b,a.y*b,a.z*b);
}

// mult by scalar from left
inline arr3x3 operator*(const double q, const arr3x3& v) { 
   return arr3x3(v.x*q,v.y*q,v.z*q); 
}

// mult by vec3d from left; should be done explicitly?
inline vec3d operator*(const vec3d& v, const arr3x3& a) { 
    return ~a * v;
}

// unary minus
inline arr3x3 operator-(arr3x3& v) { return v*-1.0; }

// unary plus
inline arr3x3 operator+(arr3x3& v) { return v*+1.0; }


inline std::ostream& operator<< (std::ostream& s, const arr3x3& a) {
    return s << "\n\t(" << a.x << ")\n\t(" << a.y << ")\n\t(" << a.z << ")";
}


#endif







