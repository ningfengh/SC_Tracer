#ifndef __raytrace__
#define __raytrace__
#include <vector>

#include "vec.h"
#include "photon_map.h"
#include "triangle.h"
#include "light.h"
#include "parameters.h"



Vec irradiance (
	Photon_map& pmap, 
	const Vec &origin, 
	const Vec &direction, 
	const std::vector<Triangle> &triangle_list, 
	int depth ,
	float RI);
Vec raytrace(
	Photon_map& pmap, 
	Photon_map& pmap_caustic, 
	const Vec &origin, 
	const Vec &direction, 
	const std::vector<Triangle> &triangle_list, 
	const std::vector<Light> &light_list, 
	int depth, 
	float RI );


#endif//__raytrace__