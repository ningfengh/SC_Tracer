#include "photon_map.h"
#include <stdlib.h> 
#include <stdio.h>  
#include <cmath>

using namespace std;

Photon_map :: Photon_map (const int max_phot){
	stored_photons = 0;
	prev_scale = 1;
	max_photons = max_phot;
	
	photons = (Photon*) malloc(sizeof(Photon)*(max_photons+1));

	if (photons == nullptr) {
		fprintf(stderr, "Out of memory initializing photon map\n");
		exit(-1);
	}

	bbox_min[0] = bbox_min[1] = bbox_min[2] = 1e8f;
	bbox_max[0] = bbox_max[1] = bbox_max[2] = -1e8f;

	for (int i = 0; i<256; i++) {
		double angle = double(i)*(1.0/256.0)*M_PI;
		costheta[i] = cos(angle);
		sintheta[i] = sin(angle);
		cosphi[i] = cos(2.0*angle);
		sinphi[i] = sin(2.0*angle);
	}
}

Photon_map :: ~Photon_map() {
	free(photons);
}

void Photon_map :: photon_dir ( float *dir, const Photon *p) const {
	dir[0] = sintheta[p->theta]*cosphi[p->phi];
	dir[1] = sintheta[p->theta]*sinphi[p->phi];
	dir[2] = costheta[p->theta];
}

void Photon_map :: irradiance_estimate(
	float irrad[3],
	const float pos[3],
	const float normal[3],
	const float max_dist,
	const int nphotons) const 
{
	irrad[0] = irrad[1] = irrad[2] = 0.0;

	NearestPhotons np;
	np.dist2 = (float*)alloca(sizeof(float)*(nphotons+1));
	np.index = (const Photon**)alloca(sizeof(Photon*)*(nphotons+1));

	np.pos[0] = pos[0]; np.pos[1] = pos[1]; np.pos[2] = pos[2];
	np.max = nphotons;
	np.found = 0;
	np.got_heap = 0;
	np.dist2[0] = max_dist*max_dist;

	locate_photons(&np,1);
	
	if (np.found<8)
		return;
	float pdir[3];

	for (int i=1;i<=np.found; i++){
		const Photon *p = np.index[i];
		photon_dir(pdir,p);
		if ((pdir[0]*normal[0]+pdir[1]*normal[1]+pdir[2]*normal[2])<0.0f) {
			irrad[0] += p->power[0];
			irrad[1] += p->power[1];
			irrad[2] += p->power[2];
		}
	}
	
	const float tmp = (1.0f/M_PI)/(np.dist2[0]);
	
	irrad[0] *= tmp;
	irrad[1] *= tmp;
	irrad[2] *= tmp;
	
}

void Photon_map :: locate_photons(
	NearestPhotons *const np,
	const int index ) const
{
	const Photon *p = &photons[index];
	float dist1; 
	
	if (index<half_stored_photons) {
		dist1 = np->pos[p->plane] - p->pos[p->plane];
		
		if (dist1>0.0) {
			locate_photons(np,2*index+1);
			if (dist1*dist1<np->dist2[0])
				locate_photons(np,2*index);
		} else {
			locate_photons(np,2*index);
			if (dist1*dist1 <np->dist2[0])
				locate_photons(np,2*index+1);
		}
	}

	dist1 = p->pos[0] - np->pos[0];
	float dist2 = dist1*dist1;
	dist1 = p->pos[1] - np->pos[1];
	dist2 += dist1*dist1;
	dist1 = p->pos[2] - np->pos[2];
	dist2 += dist1*dist1;
	
	if (dist2<np->dist2[0]) { 
		// we found a photon, Insert it in the candidate list
		if (np->found< np->max) {
			// heap is not full; use array
			np->found++;
			np->dist2[np->found] = dist2;
			np->index[np->found] = p;
		} else {
			int j,parent;
			if (np->got_heap==0) { // need to build the heap
				float dst2;
				const Photon *phot;
				int half_found = np->found>>1;
				for (int k=half_found; k>=1;k--){
					parent = k;
					phot = np->index[k];
					dst2 = np->dist2[k];
					while (parent <= half_found) {
						j = parent + parent;
						if (j<np->found && np->dist2[j]<np->dist2[j+1])
							j++;
						if (dst2>=np->dist2[j])
							break;
						np->dist2[parent] = np->dist2[j];
						np->index[parent] = np->index[j];
						parent = j;
					}
					np->dist2[parent] = dst2;
					np->index[parent] = phot;
				}
				np->got_heap = 1;
			}
			parent = 1;
			j = 2;
			while (j<=np->found) {
				if (j<np->found && np->dist2[j]<np->dist2[j+1])
					j++;
				if (dist2 > np->dist2[j])
					break;
				np->dist2[parent] = np->dist2[j];
				np->index[parent] = np->index[j];
				parent = j;
				j+=j;
			}
			np->index[parent] = p;
			np->dist2[parent] = dist2;
			
			np->dist2[0] = np->dist2[1];
		}
	}
}


void Photon_map :: store (
	const float power[3],
	const float pos[3],
	const float dir[3] )
{
	if (stored_photons>max_photons)
		return;
	stored_photons++;
	Photon *const node = &photons[stored_photons];

	for (int i = 0; i<3 ;i++){
		node->pos[i] =  pos[i];
		
		if (node->pos[i]<bbox_min[i])
			bbox_min[i] = node->pos[i];
		if (node->pos[i]>bbox_max[i])
			bbox_max[i] = node->pos[i];

		node->power[i] = power[i];
	}

	int theta = int(acos(dir[2])*(256.0/M_PI) );
	if (theta>255)
		node->theta = 255;
	else
		node->theta = (unsigned char) theta;

	int phi = int(atan2(dir[1],dir[0])*(256.0/(2.0*M_PI)) );
	if (phi>255)
		node->phi = 255;
	else if (phi<0)
		node->phi = (unsigned char)(phi+256);
	else
		node->phi = (unsigned char)phi;
}

void Photon_map :: scale_photon_power (const float scale)
{
	for (int i = prev_scale; i<=stored_photons; i++) {
		photons[i].power[0] *= scale;
		photons[i].power[1] *= scale;
		photons[i].power[2] *= scale;
	}
	prev_scale = stored_photons;
}


void Photon_map :: balance(void)
{

	if (stored_photons>1) {
		Photon **pa1 = (Photon**) malloc(sizeof(Photon*)*(stored_photons+1));
		Photon **pa2 = (Photon**) malloc(sizeof(Photon*)*(stored_photons+1));

		for (int i = 0; i<=stored_photons; i++)
			pa2[i] = &photons[i];

		balance_segment(pa1, pa2, 1, 1, stored_photons);
		free(pa2);
		
		int d, j = 1, foo = 1;
		Photon foo_photon = photons[j];

		for (int i=1; i<=stored_photons; i++) {
			d = pa1[j]-photons;
			pa1[j] = nullptr;
			if (d != foo)
				photons[j] = photons[d];
			else {
				photons[j] = foo_photon;
				if (i<stored_photons) {
					for (;foo<=stored_photons;foo++)
						if (pa1[foo] != nullptr)
							break;
					foo_photon = photons[foo];
					j = foo;
				}
				continue;
			}
			j = d;
		}
		free(pa1);
	}
	half_stored_photons = stored_photons/2 - 1;
}

#define swap(ph,a,b) {Photon *ph2 = ph[a]; ph[a] =ph[b]; ph[b]=ph2;}


void Photon_map :: median_split(
	Photon **p,
	const int start,
	const int end,
	const int median,
	const int axis)
{
	int left = start;
	int right = end;

	while (right>left) {
		const float v = p[right]->pos[axis];
		int i=left-1;
		int j=right;
		for(;;){
			while (p[++i]->pos[axis]<v)
				;
			while (p[--j]->pos[axis]>v && j>left)
				;
			if (i>=j)
				break;
			swap(p,i,j);
		}
		swap(p,i,right);
		if (i>=median)
			right = i-1;
		if (i<=median)
			left=i+1;
	}
}

void Photon_map :: balance_segment(
	Photon **pbal,
	Photon **porg,
	const int index,
	const int start,
	const int end )
{
	int median = 1;
	while ((4*median)<=(end-start+1))
		median += median;
	if ((3*median) <= (end-start+1)) {
		median += median;
		median += start - 1;
	} else
		median = end-median +1;


	int axis = 2;
	if ((bbox_max[0]-bbox_min[0])>(bbox_max[1]-bbox_min[1]) &&
		(bbox_max[0]-bbox_min[0])>(bbox_max[2]-bbox_min[2]))
		axis = 0;
	else if ((bbox_max[1]-bbox_min[1])>(bbox_max[2]-bbox_min[2]))
		axis = 1;

	median_split (porg, start, end, median, axis);
	pbal[ index ] = porg[median];
	pbal[ index ]->plane = axis;


	if (median>start) {
		if (start<median-1) {
			const float tmp = bbox_max[axis];
			bbox_max[axis] = pbal[index]->pos[axis];
			balance_segment(pbal, porg, 2*index, start, median-1);
			bbox_max[axis] = tmp;
		} else {
			pbal[2*index] = porg[start];
		}
	}
	if (median<end) {
		if (median+1<end) {
			const float tmp = bbox_min[axis];
			bbox_min[axis] = pbal[index]->pos[axis];
			balance_segment(pbal,porg, 2*index+1, median+1, end);
			bbox_min[axis] = tmp;
		} else {
			pbal[2*index+1] = porg[end];
		}
	}
}


void photon_tracing(Photon_map &pmap, const Vec origin, Vec direction, vector<Triangle>& triangle_list, Vec color, int depth, float RI){
	float t;                               // distance to intersection 
	int id=0;                               // id of intersected object 
	if (!intersect(origin, direction, triangle_list, t, id)) return; // if miss, return black
	
	const Triangle &obj = triangle_list[id];        // the hit object
 	if (++depth>5) return; //R.R.
	Vec origin_new = origin + direction*t;
	Vec n = obj.norm;
	Vec nl = n.dot(direction)<0?n:n*-1;
	Vec f = obj.surfaceColor;
	direction.norm();
	
	
	float p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl 
	
	if (obj.token == 'D' || obj.token == 'G'){                  // Ideal DIFFUSE reflection 
		if (myrand()>p){   //absorb the photon
			float color_tmp[] = {f.x*color.x*(1.f/(1-p)),f.y*color.y*(1.f/(1-p)),f.z*color.z*(1.f/(1-p))};
			float pos_tmp[] = {origin_new.x,origin_new.y,origin_new.z};
			float dir_tmp[] = {direction.x,direction.y,direction.z};
			//cout<<color_tmp[0]<<" "<<color_tmp[1]<<" "<<color_tmp[2]<<endl;
			pmap.store(color_tmp,
			           pos_tmp,
					   dir_tmp);
			return;
						
		}
		else {  // trace another ray
			float r1=2*M_PI*myrand(), r2=myrand(), r2s=sqrt(r2); 
			Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1,0):Vec(1,0,0))%w).norm(), v=w%u; 
			Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 
			photon_tracing(pmap, origin_new,d,triangle_list, f.mult(color)*(1.f/p), depth, RI);
		}
	} else if (obj.token == 'S')    // Ideal SPECULAR reflection 
		photon_tracing(pmap,origin_new,direction-n*2*n.dot(direction),triangle_list,f.mult(color)*(1.f/p),depth, RI); 
		
	else if (obj.token == 'T')   
	{        
		Vec reflRay(direction-n*2*n.dot(direction));     // Ideal dielectric REFRACTION 
		bool into = n.dot(nl)>0;                // Ray from outside going in? 
		float nc=1, nt=RI, nnt=into?nc/nt:nt/nc, ddn=direction.dot(nl), cos2t; 
		if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    {// Total internal reflection
     //cout<<"TIR"<<endl; 
			photon_tracing(pmap,origin_new,reflRay,triangle_list, color, depth, RI); 
		}
		Vec tdir = (direction*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
		float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
		float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P); 
		myrand()<P ?   // Russian roulette 
		photon_tracing(pmap, origin_new,reflRay,triangle_list,color*RP,depth, RI):
		photon_tracing(pmap, origin_new,tdir,triangle_list,color*TP, depth, RI);
	}
}


void photon_tracing_caustic(Photon_map &pmap, const Vec origin, Vec direction, vector<Triangle>& triangle_list, Vec color, int depth, bool flag, float RI, vector<Sphere>& bound){
	float t;                               // distance to intersection 
	int id=0;                               // id of intersected object 
	
	if (!flag){
		for (int i = 0; i<bound.size(); i++){
			if (!bound[i].intersect(origin, direction)) return;
		}
	}
	
	if (!intersect(origin, direction, triangle_list, t, id)) return; // if miss, return black
	if (++depth>5) return; //R.R.
	const Triangle &obj = triangle_list[id];        // the hit object
 	
	Vec origin_new = origin + direction*t;
	Vec n = obj.norm;
	Vec nl = n.dot(direction)<0?n:n*-1;
	Vec f = obj.surfaceColor;
	direction.norm();
	float p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl 
	
	
	if (flag==false && obj.token != 'T') return;
	flag = true;
	
	if (obj.token == 'D' || obj.token == 'G'){                  // Ideal DIFFUSE reflection 
		if (depth>5 || myrand()>p){   //absorb the photon
			float color_tmp[] = {f.x*color.x*(1.f/(1-p)),f.y*color.y*(1.f/(1-p)),f.z*color.z*(1.f/(1-p))};
			float pos_tmp[] = {origin_new.x,origin_new.y,origin_new.z};
			float dir_tmp[] = {direction.x,direction.y,direction.z};
			//cout<<color_tmp[0]<<" "<<color_tmp[1]<<" "<<color_tmp[2]<<endl;
			pmap.store(color_tmp,
			           pos_tmp,
					   dir_tmp);
			return;
						
		}
		else {  // trace another ray
			float r1=2*M_PI*myrand(), r2=myrand(), r2s=sqrt(r2); 
			Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1,0):Vec(1,0,0))%w).norm(), v=w%u; 
			Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 
			photon_tracing_caustic(pmap, origin_new,d,triangle_list, f.mult(color)*(1.f/p), depth, flag,RI, bound);
		}
	} else if (obj.token == 'S')    // Ideal SPECULAR reflection 
		photon_tracing_caustic(pmap,origin_new,direction-n*2*n.dot(direction),triangle_list,f.mult(color)*(1.f/p),depth, flag,RI, bound); 
		
	else if (obj.token == 'T')   
	{        
		Vec reflRay(direction-n*2*n.dot(direction));     // Ideal dielectric REFRACTION 
		bool into = n.dot(nl)>0;                // Ray from outside going in? 
		float nc=1, nt=RI, nnt=into?nc/nt:nt/nc, ddn=direction.dot(nl), cos2t; 
		if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    {// Total internal reflection
     //cout<<"TIR"<<endl; 
			photon_tracing(pmap,origin_new,reflRay,triangle_list, color, depth, RI); 
		}
		Vec tdir = (direction*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
		float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
		float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P); 
		myrand()<P ?   // Russian roulette 
		photon_tracing_caustic(pmap, origin_new,reflRay,triangle_list,color*RP,depth, flag,RI, bound):
		photon_tracing_caustic(pmap, origin_new,tdir,triangle_list,color*TP, depth,flag,RI, bound);
	}
}
