#include "file_io.h"

using namespace std;

vector<Triangle> parse(string f) {
	
	ifstream fin;
	fin.open(f);
	cout<<"----------------"<<endl;
	if(fin.fail())
	{
		cout<<"Could not open file"<<endl;
		exit(1);
	}
	vector<Triangle> triangles;
	string buffer;
	vector<float> v;
	while(!fin.eof())
	{
		getline(fin,buffer);
	    istringstream buf(buffer);
	    for(std::string token; getline(buf, token,' '); )
	       {
	       	if (token=="triangle")
	       		continue;
			if(token=="T" || token=="D"|| token=="S"|| token=="G")
			{
				//assert(v.size() == 15 || v.size()==0);
				if(v.size()==0)
					break;

				Vec v1,v2,v3,n;
				Vec surfaceColor, emissionColor;
				v1 = Vec (v[0],v[1],v[2]);
				v2 = Vec (v[3],v[4],v[5]);
				v3 = Vec (v[6],v[7],v[8]);
				surfaceColor = Vec (v[9],v[10],v[11]); 
				emissionColor = Vec(v[12],v[13],v[14]);

				char tld_token = ( v[12] != 0 ) ? 'L' : token[0];
				//char tld_token = token[0];
				Triangle tri = Triangle(v1,v2,v3, tld_token);
				tri.setColor("surface", surfaceColor);
				tri.setColor("emission", emissionColor);
					
					
				triangles.push_back(tri);
					
				v.clear();
			  	continue;
			}
	        v.push_back(stof(token));
   		}
  	}
	fin.close();
  	return triangles;	
}
