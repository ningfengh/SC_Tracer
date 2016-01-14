#ifndef __sphere_h__
#define __sphere_h__
#include "vec.h"
#include "parameters.h"

struct Sphere { 
	float rad;       // radius 
	Vec p;      // position, emission, color 
	Sphere( Vec p_, float rad_): 
     rad(rad_), p(p_) {} 
  	 float intersect(const Vec &origin, const Vec &direction) const { // returns distance, 0 if nohit 
     	Vec op = p-origin; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
     	float t, b=op.dot(direction), det=b*b-op.dot(op)+rad*rad; 
     	if (det<0) return 0; else det=sqrt(det); 
     	return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0); 
   	} 
}; 


#endif//__sphere_h__
