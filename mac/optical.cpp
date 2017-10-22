
#include "matrix.h"
#include <cmath> 
#include <math.h>
#include <iostream>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif	

using namespace std;

// eta: 0 = Ge (4.2), 1 = Zn (2.2)
double optical_filter(Vector const& d, Vector const& eta)
{
	const unsigned int layers = d.size();
	constexpr int m = 47;
	constexpr double XAir = 1.0;
	constexpr double XSub = 4.0;
	constexpr double lower = 7.7;

	constexpr double stepsize = 1000.0 * (12.3 - lower) / (m - 1);
	constexpr double ConsI = XAir * XSub * 4.0;

	double Y11, Y12, Y21, Y22;
	double refl = 0.0;

	for (unsigned int i=1; i<=m; i++)
	{
		const double lambda = lower * 1000.0 + stepsize * (i - 1.0);

		Y11 = Y22 = 1.0;
		Y12 = Y21 = 0.0;

		for (unsigned int j=1; j<=layers; j++)
		{
			const double ref = (eta[j-1]) ? 2.2 : 4.2; // refractive index
			const double delta = 2.0 * M_PI * d[j-1] * ref / lambda;

			const double X11 = std::cos(delta);
			const double X22 = X11;
			const double AA = std::sin(delta);

			const double X21 = AA * ref;
			const double X12 = AA / ref;

			const double Z11 =  Y11 * X11 - Y12 * X21;
			const double Z12 =  Y11 * X12 + Y12 * X22;
			const double Z21 =  Y21 * X11 + Y22 * X21;
			const double Z22 = -Y21 * X12 + Y22 * X22;

			Y11 = Z11;
			Y12 = Z12;
			Y21 = Z21;
			Y22 = Z22;
		}

		const double QQ1 = XAir * Y11 + XSub * Y22;
		const double QQ2 = XAir * XSub * Y12 + Y21;
		refl = refl + pow((1.0 - ConsI / (QQ1 * QQ1 + QQ2 * QQ2)), 2);

	}
	return sqrt(refl / m);
}

std::vector<double> readRow(std::istream& str)
{
    std::vector<double> result;
    std::string                line;
    std::getline(str, line);

    std::stringstream          lineStream(line);
    std::string                cell;

    while(std::getline(lineStream, cell, ','))
    {
        result.push_back(std::stod(cell));
    }
    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty())
    {
        // If there was a trailing comma then add an empty element.
        result.push_back(0.0);
    }
    return result;
}

int main(int argc, char** argv) 
{
	std::ifstream       file("input.csv");
	Vector d = Vector(readRow(file));
	Vector eta = Vector(readRow(file));

	// for (int i = 0; i <= 9; i++)  cout << d[i] << ' ';
	// cout << endl;
	// for (int i = 0; i <= 9; i++)  cout << eta[i] << ' ';
	double res = optical_filter(d, eta);

	std::ofstream out;
	out.open("result");
	out << std::setprecision(15) << res;
	out.close();
        
        cout.precision(15);
        cout << res;
	return res;
}
