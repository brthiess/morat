
#pragma once

#include <cstdio>
#include <string>

namespace Morat {

inline void logerr(std::string str){
	fprintf(stdout, "%s", str.c_str());
}

}; // namespace Morat
