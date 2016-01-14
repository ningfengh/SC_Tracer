#ifndef __triangle_h__
#define __triangle_h__
#include <string>
#include <vector>
#include "vec.h"
#include "parameters.h"

class Triangle
{
public:	
	Vec v1;   // use one vertex and two edges to save some time and space
	Vec edge1;
	Vec edge2;
	Vec norm;
	Vec surfaceColor;
	Vec emissionColor;
	char token;		// T or D or L- transparent or diffusive surface or light source surface

	Triangle();
	Triangle(Vec a,Vec b,Vec c, char tld_token);

	void computeNormal();

	bool intersection(const Vec &origin, const Vec &dir, float & t) const
	{
		Vec pVec = dir%edge2;
		float det=edge1.dot(pVec);
		if(det>-eps && det <eps)
		{
			return false;
		}

		float invDet=1.f/det;
		Vec tVec=origin-v1;
		float u=(tVec.dot(pVec))*(invDet);
		if(u<0.f || u>1.f)
		{
			return false;
		}

		Vec qVec = tVec%edge1;
		float v = dir.dot(qVec)*(invDet);

		if (v<0.f||v+u>1.f)
		{
			return false;
		}
		t = (edge2.dot(qVec))*(invDet);
		if (t>eps){
			return true;
		}
		return false;
	}	

	void setColor(std::string type, Vec color);
	
};


inline bool intersect(const Vec &origin, const Vec &direction, const std::vector<Triangle> &triangle_list, float &t, int &id){ 
   float n=triangle_list.size(), d, inf=t=1e20; 
   for(int i=0; i<n; i++) if(triangle_list[i].intersection(origin,direction,d)&&d<t){t=d;id=i;} 
   return t<inf; 
} 

#endif//__triangle_h__