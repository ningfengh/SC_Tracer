#ifndef __file_io_h__
#define __file_io_h__
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "parameters.h"
#include "triangle.h"

std::vector<Triangle> parse(std::string f);

inline float clamp(float x){ return x<0 ? 0 : x>1 ? 1 : x; } 

inline int toInt(float x){ return int(pow(clamp(x*1.5),gamma)*255+.5); } 


#endif//__file_io_h__
