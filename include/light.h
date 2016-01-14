#ifndef __light_h__
#define __light_h__

#include <vector>
#include "vec.h"
#include "triangle.h"

struct Light
{
	Vec pos;
	Vec color;
	Vec x_vec;
	Vec y_vec;
	int n_x;
	int n_y;
};



static void visualize_light(std::vector<Light>& lights,std::vector<Triangle>&surfaces) {
	for (int i = 0; i<lights.size(); i++) {
		Vec pt1(lights[i].pos-lights[i].x_vec*0.5-lights[i].y_vec*0.5);
		Vec pt2(lights[i].pos+lights[i].x_vec*0.5-lights[i].y_vec*0.5);
		Vec pt3(lights[i].pos+lights[i].x_vec*0.5+lights[i].y_vec*0.5);
		Vec pt4(lights[i].pos-lights[i].x_vec*0.5+lights[i].y_vec*0.5);
		
		Triangle t1(pt1,pt2,pt4,'L');
		Triangle t2(pt4,pt2,pt3,'L');
		t1.surfaceColor = Vec(1,1,1);
		t1.emissionColor = Vec(1,1,1);
		t2.surfaceColor = Vec(1,1,1);
		t2.emissionColor = Vec(1,1,1);
		
		surfaces.push_back(t1);
		surfaces.push_back(t2);	
	}	
}

#endif//__light_h__