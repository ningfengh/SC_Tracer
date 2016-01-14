#ifndef __vec_h__
#define __vec_h__
#include <cmath>
#include <vector>
struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm 
   float x, y, z;                  // position, also color (r,g,b) 
   Vec(float x_=0, float y_=0, float z_=0){ x=x_; y=y_; z=z_; } 
   Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); } 
   Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); } 
   Vec operator*(float b) const { return Vec(x*b,y*b,z*b); } 
   Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); } 
   Vec& norm(){ return *this = *this * (1/std::sqrt(x*x+y*y+z*z)); } 
   float dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross: 
   Vec operator%(const Vec&b) const{return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);} 
   void matMul(float m[][3]) { 
		float v[3]; 
		v[0]=x;
		v[1]=y;
		v[2]=z;
		float ans[3];
		for(int i=0;i<3;i++)
		{
			ans[i]=0;
		}
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				//std::cout<<m[i][j]<<" "<<v[j]<<"\n";
				ans[i]+=m[i][j]*v[j];
			}
		}
		x=ans[0];
		y=ans[1];
		z=ans[2];
		}
}; 

inline float normalize(Vec v)
{
	return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

#endif//__vec_h__