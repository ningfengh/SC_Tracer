#include "triangle.h"
using namespace std;
Triangle::Triangle()
{
	v1 = Vec(0,0,0);
	token = 'T';
	surfaceColor = Vec(0,0,0);
	emissionColor =  Vec(0,0,0);
	computeNormal();
}

Triangle::Triangle(Vec a,Vec b,Vec c, char tld_token)
{
	v1 = a;
	Vec v2 = b;
	Vec v3 = c;
	edge1 = v2-v1;
	edge2 = v3-v1;
	token = tld_token;
	surfaceColor = Vec(0,0,0);
	emissionColor =  Vec(0,0,0);
	computeNormal();
}
	
void Triangle::computeNormal()
{
	norm = edge2%edge1;
	norm.norm();
}


	
	
void Triangle::setColor(string type, Vec color)
{
	if ( type == "surface" )
		surfaceColor = color;
	else if ( type == "emission" )
		emissionColor = color;
}