#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <curand.h>
#include <curand_kernel.h>

#define BLOCK_SIZE 		16
#define MAX_TRIANGLE	100
#define MAX_LIGHT 		10
#define eps 			0.0001
#define MAX_RAY_DEPTH 	20
#define	AAKERNEL_SIZE	6


using namespace std;

__device__ unsigned int WangHash(unsigned int a) {
    a = (a ^ 61) ^ (a >> 16);
    a = a + (a << 3);
    a = a ^ (a >> 4);
    a = a * 0x27d4eb2d;
    a = a ^ (a >> 15);
    return a;
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

struct Vec {        
   float x, y, z;                 
   __host__ __device__ Vec(){} 
   __host__ __device__ Vec(float x_, float y_, float z_){ x=x_; y=y_; z=z_; } 
   __host__ __device__ Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); } 
   __host__ __device__ Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); } 
   __host__ __device__ Vec operator*(float b) const { return Vec(x*b,y*b,z*b); } 
   __host__ __device__ Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); } 
   __host__ __device__ Vec& norm(){ return *this = *this * (1/sqrtf(x*x+y*y+z*z)); } 
   __host__ __device__ float dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross: 
   __host__ __device__ Vec operator%(const Vec&b) const{return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);} 
}; 

__device__ float normalize(Vec v)
{
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

struct Parameter {
	int w;
	int h;
	int samps;
	int n_triangles;
	int n_lights;
	float fov;
    float aspectratio;
	float angle;
};

struct Tracing_Stack{
	Vec o;
	Vec d;
	Vec pre_color;
	int depth;
	__device__ Tracing_Stack(){};	
	__device__ Tracing_Stack(const Vec &o,
					  const Vec &d,
					  const Vec &pre_color,
					  const int depth){
		this->o = o;
		this->d = d;
		this->pre_color = pre_color;
		this->depth = depth;
	}
};

struct Light
{
	Vec pos;
	Vec color;
	Vec x_vec;
	Vec y_vec;
	int n_x;
	int n_y;
	
	
};


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
    __host__ __device__ Triangle(){};
	__host__ Triangle(const Vec &a,
					  const Vec &b,
					  const Vec &c, 
					  const Vec &sColor, 
					  const Vec &eColor, 
					  char tld_token){
		v1 = a;
		Vec v2 = b;
		Vec v3 = c;
		edge1 = v2-v1;
		edge2 = v3-v1;
		token = tld_token;
		surfaceColor = sColor;
		emissionColor =  eColor;
		computeNormal();
	}

	__host__ void computeNormal()
	{
		norm = edge2%edge1;
		norm.norm();		
	}

	__device__ bool intersection(const Vec &origin, const Vec &dir, float & t) const
	{
		Vec pVec = dir%edge2;
		float det=edge1.dot(pVec);
		//if(det>-eps && det <eps)
		if(det==0)
		{
			return false;
		}

		float invDet=1./det;
		Vec tVec=origin-v1;
		float u=(tVec.dot(pVec))*(invDet);
		if(u<0. || u>1.)
		{
			return false;
		}

		Vec qVec = tVec%edge1;
		float v = dir.dot(qVec)*(invDet);

		if (v<0.||v+u>1.)
		{
			return false;
		}
		t = (edge2.dot(qVec))*(invDet);
		if (t>eps){
			return true;
		}
		return false;
	}	

};

__constant__ Triangle ctriangles[MAX_TRIANGLE];
__constant__ Parameter cparam[1];
__constant__ Light clights[MAX_LIGHT];
__constant__ float AAFilter[AAKERNEL_SIZE][3] 	=		/* X, Y, coef */
	{
		-0.52, 0.38, 0.128,
		0.41, 0.56, 0.119,
		0.27, 0.08, 0.294,
		-0.17, -0.29, 0.249,
		0.58, -0.55, 0.104,
		-0.31, -0.71, 0.106
	};	
__host__ void parse(string file_name, Triangle* triangles, int &n_triangles)
{
	ifstream fin;
	fin.open(file_name);
	int cnt_v  = 0;
	float v[15];

	if(fin.fail())
	{
		cout<<"Could not open file"<<endl;
		exit(1);
	}

	string buffer;
	n_triangles = 0;
	while(!fin.eof())
	{
		getline(fin,buffer);
	    istringstream buf(buffer);
	    for(string token; getline(buf, token,' '); )
        {
        	if (token=="triangle") {
        		cnt_v = 0;
			}
			else if(token=="T" || token=="D"|| token=="S")
			{
				if (n_triangles>=MAX_TRIANGLE){
					cout<<"Number of triangles should be equal or less than "<<MAX_TRIANGLE<<endl;
					exit(1);
				}
				Vec v1(v[0],v[1],v[2]);
				Vec v2(v[3],v[4],v[5]);
				Vec v3(v[6],v[7],v[8]);
				Vec sColor(v[9],v[10],v[11]);
				Vec eColor(v[12],v[13],v[14]);
				Triangle tri_tmp(v1,v2,v3,sColor,eColor,token[0]);
				triangles[n_triangles++] = tri_tmp;
				
			}
        	else v[cnt_v++] = stof(token);
		}
	}
	fin.close();
	return;
}

__host__ float clamp(float x){ return x<0 ? 0 : x>1 ? 1 : x; }
__host__ int toInt(float x){ return int(powf(clamp(x),1/2.2)*255+.5); }


__global__ void init_rand(curandState *state, unsigned int seed) {
	int idx_x = blockIdx.x * blockDim.x + threadIdx.x;
	int idx_y = blockIdx.y * blockDim.y + threadIdx.y;
	int idx = idx_y*(cparam[0].w)+idx_x;
	curand_init(seed + WangHash(idx),0 , 0, &state[idx]);
}


__device__ bool intersect(const Vec &origin, const Vec &direction, float &t, int &id){
	   int n = cparam[0].n_triangles;  
	   float d, inf=1e5;
	   t = inf;
	   for(int i=0; i<n; i++) if((ctriangles[i].intersection(origin,direction,d))&&d<t){t=d;id=i;} 
	   return (t<inf);
}

__device__ Vec raytrace(Vec &o, Vec &d) {
	
	Tracing_Stack stk[MAX_RAY_DEPTH+1];
	int stk_cnt = 1;
	stk[0] = Tracing_Stack(o,d,Vec(1,1,1),0);
	Vec color(0,0,0);
	
	do{
		float t;
		int id;
		if (!intersect(stk[stk_cnt-1].o, stk[stk_cnt-1].d, t, id)) {stk_cnt--; continue;} //no hit
		const Triangle &obj = ctriangles[id];        // the hit object
		Vec new_o = stk[stk_cnt-1].o + stk[stk_cnt-1].d*t;    // update the origin
		Vec n = obj.norm;
		Vec nl = n.dot(stk[stk_cnt-1].d)<0?n:n*-1;
		Vec f = obj.surfaceColor;
	   
		//float p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl 
		if (stk[stk_cnt-1].depth>=MAX_RAY_DEPTH) {
			 {color = color + stk[stk_cnt-1].pre_color.mult(obj.emissionColor);stk_cnt--;continue;} //R.R.
		}
		
		if (obj.token == 'D'){                  // Ideal DIFFUSE reflection 
			
			Vec col(0,0,0);
	
			for (int i = 0; i < cparam[0].n_lights;i++){
				float factor = 1./clights[i].n_x/clights[i].n_y;
				for (int j = 0; j<clights[i].n_x; j++) {
					for (int k = 0; k<clights[i].n_y; k++) {
						Vec l_pos = clights[i].pos - clights[i].x_vec*0.5 + clights[i].x_vec * (1./clights[i].n_x*j) 
											  - clights[i].y_vec*0.5 + clights[i].y_vec * (1./clights[i].n_y*k);
						Vec d = (l_pos - new_o);
						float t_light = normalize(d);
						d = d.norm();
						int id = 0;
						if (!intersect(new_o, d, t, id) || ctriangles[id].token=='L' || t>t_light) {
							col = col + f.mult(clights[i].color)*(d.dot(obj.norm))*factor;	
						}
					}
				}
			}
			color = color + stk[stk_cnt-1].pre_color.mult(col);stk_cnt--;
			
			continue;
		}
		
		
		else if (obj.token == 'S'){            // Ideal SPECULAR reflection 
   			color = color+stk[stk_cnt-1].pre_color.mult(obj.emissionColor);
			stk[stk_cnt-1].o = new_o;
			stk[stk_cnt-1].d = stk[stk_cnt-1].d - n*2*n.dot(stk[stk_cnt-1].d);
			stk[stk_cnt-1].pre_color = stk[stk_cnt-1].pre_color.mult(f);
			stk[stk_cnt-1].depth++;
			continue;
		}
		
				
   		Vec reflRay(stk[stk_cnt-1].d-n*2*n.dot(stk[stk_cnt-1].d));     // Ideal dielectric REFRACTION
		
   		bool into = n.dot(nl)>0;                // Ray from outside going in? 
   		float nc=1, nt=2.4, nnt=into?nc/nt:nt/nc, ddn=stk[stk_cnt-1].d.dot(nl), cos2t; 
   		if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    {// Total internal reflection
     //cout<<"TIR"<<endl;
			color = color + stk[stk_cnt-1].pre_color.mult(obj.emissionColor);
			stk[stk_cnt-1].o = new_o;
			stk[stk_cnt-1].d = reflRay;
			stk[stk_cnt-1].pre_color = stk[stk_cnt-1].pre_color.mult(f);
			stk[stk_cnt-1].depth++;
			continue;
   		}
	 	Vec tdir = (stk[stk_cnt-1].d*nnt - n*((into?1:-1)*(ddn*nnt+sqrtf(cos2t)))).norm(); 
  		float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
  		float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re;
		
		color = color + stk[stk_cnt-1].pre_color.mult(obj.emissionColor);
		stk[stk_cnt-1].o = new_o;
		stk[stk_cnt-1].d = reflRay;
		stk[stk_cnt-1].pre_color = stk[stk_cnt-1].pre_color.mult(f);
		stk[stk_cnt-1].depth++;
		stk_cnt++;
		stk[stk_cnt-1].o = new_o;
		stk[stk_cnt-1].d = tdir;
		stk[stk_cnt-1].pre_color = stk[stk_cnt-2].pre_color;
		stk[stk_cnt-1].depth = stk[stk_cnt-2].depth;
		
		stk[stk_cnt-2].pre_color = stk[stk_cnt-2].pre_color*Re;
		stk[stk_cnt-1].pre_color = stk[stk_cnt-1].pre_color*Tr;
		
	} while (stk_cnt);
		
	
	
	return color;
}




__global__ void path_tracing(Vec *d_c) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int idx = y*(cparam[0].w)+x;

	Vec dr;
	d_c[idx] = Vec(0,0,0);
	for (int i = 0; i<AAKERNEL_SIZE; i++){
		dr.x = (2. * ((x+0.5+AAFilter[i][0])/cparam[0].w) -1.  )*cparam[0].angle*cparam[0].aspectratio;
		float temp = (1. - 2.*((y+0.5+AAFilter[i][1])/cparam[0].h)) ;
		dr.y=temp*cparam[0].angle;
		dr.z = 1.;

		dr.norm();

		Vec dr_origin(0,0,0);

		d_c[idx] = d_c[idx] + (raytrace(dr_origin, dr)*AAFilter[i][2]);
	}
}

int main(int argc, char *argv[]){
	Parameter hparam;
	hparam.w = 1024;
	hparam.h = 1024;
	hparam.samps = argc==2 ? atoi(argv[1]) : 500; // # samples

	if (hparam.w%BLOCK_SIZE) {
		hparam.w = (hparam.w/BLOCK_SIZE+1)*BLOCK_SIZE;
		cout<<"Width has been changed to "<<hparam.w<<endl;
	}
	if (hparam.h%BLOCK_SIZE) {
		hparam.h = (hparam.h/BLOCK_SIZE+1)*BLOCK_SIZE;
		cout<<"Height has been changed to "<<hparam.h<<endl;
	}

	hparam.fov = 40.0;
    hparam.aspectratio = hparam.w/hparam.h;
	hparam.angle = tanf(0.5*hparam.fov*M_PI/180.0);
	
	
	Triangle htriangles[MAX_TRIANGLE];
	Light hlights[MAX_LIGHT];
	parse("prism_oct_no_light.asc", htriangles, hparam.n_triangles);
	
	hparam.n_lights = 1;
	
	hlights[0].pos = Vec(1.6,2.749,10.75);
	hlights[0].color = Vec(1,1,1);
	hlights[0].x_vec = Vec(1.2,0,0);
	hlights[0].y_vec = Vec(0,0,1.2);
	hlights[0].n_x = 9;
	hlights[0].n_y = 9;
	
	gpuErrchk(cudaSetDevice(0));	

	gpuErrchk(cudaMemcpyToSymbol(ctriangles,  htriangles,   sizeof(Triangle)*MAX_TRIANGLE));
	gpuErrchk(cudaMemcpyToSymbol(cparam , &hparam, sizeof(Parameter)));
	gpuErrchk(cudaMemcpyToSymbol(clights, &hlights, sizeof(Light)*MAX_LIGHT));
	
	Vec *c;
	Vec *d_c;

	c = (Vec*)malloc((hparam.w)*(hparam.h)*sizeof(Vec));

	gpuErrchk(cudaMalloc((void**) &d_c, (hparam.w)*(hparam.h)*sizeof(Vec)));

	dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE);
	dim3 dimGrid(hparam.w/BLOCK_SIZE,hparam.h/BLOCK_SIZE);

	
	curandState_t* states;
	gpuErrchk(cudaMalloc((void**) &states, (hparam.w)*(hparam.h) * sizeof(curandState_t)));
	init_rand<<<dimGrid,dimBlock>>>(states, time(0));
	gpuErrchk( cudaPeekAtLastError());
	

	path_tracing<<<dimGrid,dimBlock>>>(d_c);
	gpuErrchk( cudaPeekAtLastError());
	
	
	gpuErrchk(cudaMemcpy(c, d_c, (hparam.w)*(hparam.h)*sizeof(Vec), cudaMemcpyDeviceToHost));
	
	FILE *f = fopen("image_ray.ppm", "w");         // Write image to PPM file.
   	fprintf(f, "P3\n%d %d\n%d\n", hparam.w, hparam.h, 255);
   	for (int i=0; i<(hparam.w)*(hparam.h); i++) {
		fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	}
	fclose(f);
	
	free(c);
	cudaFree(d_c);
	cudaFree(states);
	
	
	
	return 0;
}
