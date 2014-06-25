
#pragma once

#include <cstdio>
#include <string>

inline void logerr(std::string str){
	fprintf(stdout, "%s", str.c_str());
}

