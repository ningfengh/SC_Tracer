#include "raytrace.h"

using namespace std;

Vec irradiance (
	Photon_map& pmap, 
	const Vec &origin, 
	const Vec &direction, 
	const vector<Triangle> &triangle_list, 
	int depth ,
	float RI)
{
	float t;                               // distance to intersection 
	int id=0;                               // id of intersected object 

	if (!intersect(origin, direction, triangle_list, t, id)) return Vec(); // if miss, return black
	
	const Triangle &obj = triangle_list[id];        // the hit object
 
	Vec origin_new = origin + direction*t;
	Vec n = obj.norm;
	Vec nl = n.dot(direction)<0?n:n*-1;
	Vec f = obj.surfaceColor;


	if (++depth>5) return Vec(); //R.R.
	
	if (obj.token == 'D'){      
		Vec col;
	    //direct visualization of the photon map
		float color[3];
		float pos[3]={origin_new.x,origin_new.y,origin_new.z};
		float normal[3] = {n.x, n.y, n.z};
		pmap.irradiance_estimate(color,pos,normal,0.1,100);
		col = Vec(color[0],color[1],color[2]);	
		return col;	
	}
	else if (obj.token == 'S')            // Ideal SPECULAR reflection 
		return f.mult(irradiance(pmap, origin_new,direction-n*2*n.dot(direction),triangle_list,depth, RI)); 
	else if (obj.token == 'T') {
		Vec reflRay(direction-n*2*n.dot(direction));     // Ideal dielectric REFRACTION 
		bool into = n.dot(nl)>0;                // Ray from outside going in? 
		float nc=1, nt=RI, nnt=into?nc/nt:nt/nc, ddn=direction.dot(nl), cos2t; 
		if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    {// Total internal reflection
			return f.mult(irradiance(pmap, origin_new,reflRay,triangle_list,depth, RI)); 
		}

		Vec tdir = (direction*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
		float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
		float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re; 
		return obj.emissionColor + f.mult(irradiance(pmap, origin_new,reflRay,triangle_list,depth, RI)*Re+irradiance(pmap, origin_new,tdir,triangle_list,depth, RI)*Tr); 
	}
	return Vec();
	 
}



Vec raytrace(
	Photon_map& pmap, 
	Photon_map& pmap_caustic, 
	const Vec &origin, 
	const Vec &direction, 
	const vector<Triangle> &triangle_list, 
	const vector<Light> &light_list, 
	int depth, 
	float RI )
{ 
 	
	float t;                               // distance to intersection 
	int id=0;                               // id of intersected object 
	
	if (!intersect(origin, direction, triangle_list, t, id)) return Vec(); // if miss, return black
	
	const Triangle &obj = triangle_list[id];        // the hit object
 
	Vec origin_new = origin + direction*t;
	Vec n = obj.norm;
	Vec nl = n.dot(direction)<0?n:n*-1;
	Vec f = obj.surfaceColor;
   
	//float p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl 
  
	if (++depth>20) return obj.emissionColor; //R.R. 
 
	if (obj.token == 'L'){
		//return Vec(0,0,0);
		return obj.emissionColor;
	}
	else if (obj.token == 'D' || obj.token == 'G'){                  // Ideal DIFFUSE reflection 
		Vec col(0,0,0);


/****************   Direct visualization of Caustic ******************/

		float color[3];
		float pos[3]={origin_new.x,origin_new.y,origin_new.z};
		float normal[3] = {n.x, n.y, n.z};
		pmap_caustic.irradiance_estimate(color,pos,normal,0.1,100);
		col = Vec(color[0],color[1],color[2]);

/*********************************************************************/

	
 /***************   Global illumination   ****************************/
		int nsamps = 200;
		for (int i = 0; i<nsamps ; i++)
		{
			double r1=2*M_PI*myrand(), r2=myrand(), r2s=sqrt(r2); 
     		Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1,0):Vec(1,0,0))%w).norm(), v=w%u; 
     		Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 		
			col = col + irradiance(pmap,origin_new,d,triangle_list, 0, RI)*(1.f/nsamps);
		}
/*********************************************************************/


/****************   Direct illumination  *****************************/
		for (int i = 0; i < light_list.size();i++){
			float factor = 1./light_list[i].n_x/light_list[i].n_y;
			for (int j = 0; j<light_list[i].n_x; j++) {
				for (int k = 0; k<light_list[i].n_y; k++) {
					Vec l_pos = light_list[i].pos - light_list[i].x_vec*0.5 + light_list[i].x_vec * (1./light_list[i].n_x*j) 
											  - light_list[i].y_vec*0.5 + light_list[i].y_vec * (1./light_list[i].n_y*k);
					Vec d = (l_pos - origin_new);
					float t_light = normalize(d);
					d = d.norm();
					int id = 0;
					if (!intersect(origin_new, d, triangle_list, t, id) || triangle_list[id].token=='L' || t>t_light) {
						col = col + f.mult(light_list[i].color)*(d.dot(obj.norm))*factor;	
	
					}
				}
			}
		}
		
		if (obj.token == 'G') col = col + f.mult(raytrace(pmap, pmap_caustic, origin_new,direction-n*2*n.dot(direction),triangle_list,light_list,depth,RI)); 
/**********************************************************************/		
		return col; 

	} else if (obj.token == 'S')            // Ideal SPECULAR reflection 
		return f.mult(raytrace(pmap, pmap_caustic, origin_new,direction-n*2*n.dot(direction),triangle_list,light_list,depth,RI)); 

   
	Vec reflRay(direction-n*2*n.dot(direction));     // Ideal dielectric REFRACTION 
	bool into = n.dot(nl)>0;                // Ray from outside going in? 
	float nc=1, nt=RI, nnt=into?nc/nt:nt/nc, ddn=direction.dot(nl), cos2t; 
	if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    {// Total internal reflection
		return obj.emissionColor + f.mult(raytrace(pmap, pmap_caustic, origin_new,reflRay,triangle_list,light_list,depth,RI)); 
	}

	Vec tdir = (direction*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
	float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
	float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re; 
	return obj.emissionColor + f.mult(raytrace(pmap, pmap_caustic, origin_new,reflRay,triangle_list,light_list,depth,RI)*Re
                                     +raytrace(pmap, pmap_caustic, origin_new,tdir,triangle_list,light_list,depth,RI)*Tr); 
} 
