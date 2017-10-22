
#include "matrix.h"
#ifdef _WIN32
#include <algorithm>	// otherwise "'min' : is not a member of 'std'"  
#endif

using namespace std;


Matrix operator * (double lhs, Matrix const& rhs)
{
	return (rhs * lhs);
}

std::ostream& operator << (std::ostream& str, Matrix const& mat)
{
	str << "[";
	if (! mat.empty())
	{
		for (size_t r=0; r<mat.rows(); r++)
		{
			if (r != 0) str << "; ";
			str << mat(r, 0);
			for (size_t c=1; c<mat.cols(); c++) str << ", " << mat(r, c);
		}
	}
	str << "]";
	return str;
}
