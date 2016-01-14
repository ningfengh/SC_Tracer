#ifndef __parameters_h__
#define __parameters_h__
#include <random>
#include <ctime>
#include <functional> 
#define gamma (1.f/1.5)
#define eps 1e-4

static std::default_random_engine generator(time(nullptr));
static std::uniform_real_distribution<float> distribution(0,1);
static auto myrand = std::bind(distribution,generator);

#endif//__parameters_h__