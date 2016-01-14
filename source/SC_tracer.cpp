#include <iostream>
#include <string>
#include <vector>
#include <cmath>   
#include <stdlib.h> 
#include <stdio.h>  

#include "vec.h"
#include "photon_map.h"
#include "file_io.h"
#include "light.h"
#include "triangle.h"
#include "sphere.h"
#include "raytrace.h"
#include "parameters.h"

#define	AAKERNEL_SIZE	6



using namespace std;



int main(int argc, char *argv[]){

	int w=1000, h=1000; // # samples 
	float fov = 20.0;
    float aspectratio = w/h;
	float angle = tan(0.5*fov*M_PI/180.0);

	float AAFilter[AAKERNEL_SIZE][3] 	=		/* X, Y, coef */
	{
		-0.52, 0.38, 0.128,
		0.41, 0.56, 0.119,
		0.27, 0.08, 0.294,
		-0.17, -0.29, 0.249,
		0.58, -0.55, 0.104,
		-0.31, -0.71, 0.106
	};	

	float angle_x = 25;
	float Rx[3][3] =
    {
		{1,0,0},
		{0,cos(angle_x*M_PI/180),-sin(angle_x*M_PI/180)},
		{0,sin(angle_x*M_PI/180),cos(angle_x*M_PI/180)}
    };
	
	vector<Triangle> surfaces = parse(argv[1]);
	vector<Light> lights;
/*
	Light light1, light2;
	
	light1.pos = Vec(1.6,2.749,10.75);
	light1.color = Vec(0.5,0.5,0.5);
	light1.x_vec = Vec(1.2,0,0);
	light1.y_vec = Vec(0,0,1.2);
	light1.n_x = 9;
	light1.n_y = 9;
	
	lights.push_back(light1); */


	Light light1, light2;
	float total_area;
	float light_area[]={100,100};
	light1.pos = Vec(0,7.749,10.75);
	light1.color = Vec(1,1,1);
	light1.x_vec = Vec(10.0,0,0);
	light1.y_vec = Vec(0,0,10.0);
	light1.n_x = 9;
	light1.n_y = 9;
	
	light2.pos = Vec(-1,3,6);
	light2.color = Vec(1,1,1);
	light2.x_vec = Vec(1.0,0,0)*10;
	light2.y_vec = ((light2.pos-Vec(0,-2.75,10.75))%light2.x_vec).norm()*10;
	light2.n_x = 9;
	light2.n_y = 9;	
	
	lights.push_back(light1);
	lights.push_back(light2);
	
	total_area = 10*10 + 10*10;
	
	visualize_light(lights,surfaces);

	vector<Sphere> bound;
	Sphere diamond(Vec(0,-1.9233,10.75),1.107);
	bound.push_back(diamond);

	int n_photons = 1;
	int n_photons_caustic = 1000000;
	Photon_map pmap_r(n_photons),pmap_g(n_photons),pmap_b(n_photons) ;
	Photon_map pmap_caustic_r(n_photons_caustic),pmap_caustic_g(n_photons_caustic),pmap_caustic_b(n_photons_caustic) ;
		
	
	for (int i = 0; i<lights.size(); i++){
	
		while (pmap_r.stored_photons<(float)pmap_r.max_photons*light_area[i]/total_area) {
			float offset_x = myrand();
			float offset_y = myrand();
			Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
								-lights[i].y_vec*0.5+lights[i].y_vec*offset_y;
			Vec nl(0,-1,0);
			float r1=2*M_PI*myrand(), r2=myrand(), r2s=sqrt(r2); 
			Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1,0):Vec(1,0,0))%w).norm(), v=w%u; 
			Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
			Vec color(10,10,10);
			photon_tracing(pmap_r,pos,d,surfaces,color,0,2.40);
		}

		cout<<"finish generating red photon map"<<endl;

		while (pmap_caustic_r.stored_photons<(float)pmap_caustic_r.max_photons*light_area[i]/total_area) {
			float offset_x = myrand();
			float offset_y = myrand();
			Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
								-lights[i].y_vec*0.5+lights[i].y_vec*offset_y;
			Vec nl(0,-1,0);
			float r1=2*M_PI*myrand(), r2=myrand(), r2s=sqrt(r2); 
			Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1,0):Vec(1,0,0))%w).norm(), v=w%u; 
			Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
			Vec color(20,20,20);
			photon_tracing_caustic(pmap_caustic_r,pos,d,surfaces,color,0,false, 2.40, bound);
		}

		cout<<"finish generating caustic red photon map"<<endl;

		while (pmap_g.stored_photons<(float)pmap_g.max_photons*light_area[i]/total_area) {
			float offset_x = myrand();
			float offset_y = myrand();
			Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
								-lights[i].y_vec*0.5+lights[i].y_vec*offset_y;
			Vec nl(0,-1,0);
			float r1=2*M_PI*myrand(), r2=myrand(), r2s=sqrt(r2); 
			Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1,0):Vec(1,0,0))%w).norm(), v=w%u; 
			Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
			Vec color(10,10,10);
			photon_tracing(pmap_g,pos,d,surfaces,color,0,2.43);
		}

		cout<<"finish generating green photon map"<<endl;

		while (pmap_caustic_g.stored_photons<(float)pmap_caustic_g.max_photons*light_area[i]/total_area) {
			float offset_x = myrand();
			float offset_y = myrand();
			Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
									-lights[i].y_vec*0.5+lights[i].y_vec*offset_y;
			Vec nl(0,-1,0);
			float r1=2*M_PI*myrand(), r2=myrand(), r2s=sqrt(r2); 
			Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1,0):Vec(1,0,0))%w).norm(), v=w%u; 
			Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
			Vec color(20,20,20);
			photon_tracing_caustic(pmap_caustic_g,pos,d,surfaces,color,0,false, 2.43, bound);
		}


		cout<<"finish generating causitc green photon map"<<endl;

		while (pmap_b.stored_photons<(float)pmap_b.max_photons*light_area[i]/total_area) {
			float offset_x = myrand();
			float offset_y = myrand();
			Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
									-lights[i].y_vec*0.5+lights[i].y_vec*offset_y;
			Vec nl(0,-1,0);
			float r1=2*M_PI*myrand(), r2=myrand(), r2s=sqrt(r2); 
			Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1,0):Vec(1,0,0))%w).norm(), v=w%u; 
			Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
			Vec color(10,10,10);
			photon_tracing(pmap_b,pos,d,surfaces,color,0,2.46);
		}

		cout<<"finish generating blue photon map"<<endl;
		while (pmap_caustic_b.stored_photons<(float)pmap_caustic_b.max_photons*light_area[i]/total_area) {
			float offset_x = myrand();
			float offset_y = myrand();
			Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
									-lights[i].y_vec*0.5+lights[i].y_vec*offset_y;
			Vec nl(0,-1,0);
			float r1=2*M_PI*myrand(), r2=myrand(), r2s=sqrt(r2); 
			Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1,0):Vec(1,0,0))%w).norm(), v=w%u; 
			Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
			Vec color(10,10,10);
			photon_tracing_caustic(pmap_caustic_b,pos,d,surfaces,color,0,false, 2.46, bound);
		}

		cout<<"finish generating caustic blue photon map"<<endl;

		
	}
	pmap_r.balance();
	pmap_r.scale_photon_power(1.f/n_photons);	
	pmap_caustic_r.balance();
	pmap_caustic_r.scale_photon_power(1.f/n_photons_caustic);
	pmap_g.balance();
	pmap_g.scale_photon_power(1.f/n_photons);	
	pmap_caustic_g.balance();
	pmap_caustic_g.scale_photon_power(1.f/n_photons_caustic);
	pmap_b.balance();
	pmap_b.scale_photon_power(1.f/n_photons);		
	pmap_caustic_b.balance();
	pmap_caustic_b.scale_photon_power(1.f/n_photons_caustic);	

	visualize_light(lights,surfaces);
	
	

	Vec *c=new Vec[w*h];
	Vec output_color;
	
	
 #pragma omp parallel for schedule(dynamic, 1) private(output_color)      // OpenMP 
	for ( int y = 0 ; y < h ; y++ )
	{
		fprintf(stderr,"\rRendering %5.2f%%",100.*y/(h-1));
		for ( unsigned short x = 0; x < w ; x++ )
		{
			int i=(y)*w+x;
			output_color=Vec();
			for (int aa = 0; aa < AAKERNEL_SIZE; aa++)
			{
				Vec dr;
				dr.x = (2 * ((x+0.5+AAFilter[aa][0])/w) -1  )*angle*aspectratio;
				float temp = (1 - 2*((y+0.5+AAFilter[aa][1])/h)) ;
				
				dr.y=temp*angle;
				dr.z = 1;
				dr.norm();	
				dr.matMul(Rx);		
				//Vec dr_origin(0,0,0);
				Vec dr_origin = dr*9+Vec(-0.1,3,0);
				output_color = output_color+(Vec(1,0,0).mult(raytrace(pmap_r,pmap_caustic_r,dr_origin,dr,surfaces,lights,0,2.40))
										   + Vec(0,1,0).mult(raytrace(pmap_g,pmap_caustic_g,dr_origin,dr,surfaces,lights,0,2.43))
										   + Vec(0,0,1).mult(raytrace(pmap_b,pmap_caustic_b,dr_origin,dr,surfaces,lights,0,2.46)))*AAFilter[aa][2];
				
			}
			c[i] = output_color;
		}
	}

	FILE *f = fopen(argv[2], "w");         // Write image to PPM file. 
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int i=0; i<w*h; i++) {
		fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)); 
	}
	
	return 0;
}

