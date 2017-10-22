
#include "vector.h"
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <cmath>


// check index ranges and dimension only in debug mode
#ifdef DEBUG
#define CHECKRANGE
#endif


using namespace std;


Vector::Vector()
: m_size(0)
, m_data(nullptr)
, m_view(false)
{ }

Vector::Vector(size_t dim)
: m_size(dim)
, m_data((double*)malloc(sizeof(double) * dim))
, m_view(false)
{ }

Vector::Vector(double scalar)
: m_size(1)
, m_data((double*)malloc(sizeof(double)))
, m_view(false)
{ m_data[0] = scalar; }

Vector::Vector(size_t dim, double value)
: m_size(dim)
, m_data((double*)malloc(sizeof(double) * dim))
, m_view(false)
{ for (size_t i=0; i<dim; i++) m_data[i] = value; }

Vector::Vector(std::size_t dim, const double* data)
: m_size(dim)
, m_data((double*)malloc(sizeof(double) * dim))
, m_view(false)
{ for (size_t i=0; i<dim; i++) m_data[i] = data[i]; }

Vector::Vector(std::vector<double> const& other)
: m_size(other.size())
, m_data((double*)malloc(sizeof(double) * other.size()))
, m_view(false)
{ for (size_t i=0; i<other.size(); i++) m_data[i] = other[i]; }

Vector::Vector(double* data, std::size_t dim)
: m_size(dim)
, m_data(data)
, m_view(true)
{ }

Vector::Vector(Vector const& other)
: m_size(other.size())
, m_data(other.isView() ? const_cast<double*>(other.data()) : (double*)malloc(sizeof(double) * other.size()))
, m_view(other.isView())
{
	if (! other.isView()) memmove(m_data, other.data(), sizeof(double) * m_size);
}

Vector::Vector(std::initializer_list<double> l)
: m_size(l.size())
, m_data((double*)malloc(sizeof(double) * l.size()))
, m_view(false)
{
	std::size_t i=0;
	for (auto v : l) m_data[i++] = v;
}

Vector::~Vector()
{
	if (! m_view && m_data) free(m_data);
}


Vector& Vector::operator = (Vector const& rhs)
{
	if (m_view)
	{
#ifdef CHECKRANGE
		if (rhs.size() > m_size) throw runtime_error("dimension mismatch");
#endif
		memmove(m_data, rhs.data(), sizeof(double) * m_size);
	}
	else
	{
		free(m_data);
		m_size = rhs.size();
		m_data = (double*)malloc(sizeof(double) * m_size);
		memmove(m_data, rhs.m_data, sizeof(double) * m_size);
	}
	return *this;
}

Vector Vector::copy() const
{
	Vector ret(m_size, 0.0);
	for (size_t i=0; i<m_size; i++) ret[i] = (*this)[i];
	return ret;
}

double& Vector::operator [] (size_t index)
{
#ifdef CHECKRANGE
	if (index >= m_size) throw runtime_error("index out of bounds");
#endif
	return m_data[index];
}

double const& Vector::operator [] (size_t index) const
{
#ifdef CHECKRANGE
	if (index >= m_size) throw runtime_error("index out of bounds");
#endif
	return m_data[index];
}

double& Vector::operator () (size_t index)
{
#ifdef CHECKRANGE
	if (index >= m_size) throw runtime_error("index out of bounds");
#endif
	return m_data[index];
}

double const& Vector::operator () (size_t index) const
{
#ifdef CHECKRANGE
	if (index >= m_size) throw runtime_error("index out of bounds");
#endif
	return m_data[index];
}

Vector Vector::sub(std::size_t start, std::size_t end, bool view)
{
	if (view) return Vector(m_data + start, end - start);
	else return Vector(m_data + start, end - start).copy();
}

void Vector::operator += (Vector const& arg)
{
#ifdef CHECKRANGE
	if (arg.size() != m_size) throw runtime_error("dimension mismatch");
#endif
	for (size_t i=0; i<m_size; i++) m_data[i] += arg.m_data[i];
}

void Vector::operator -= (Vector const& arg)
{
#ifdef CHECKRANGE
	if (arg.size() != m_size) throw runtime_error("dimension mismatch");
#endif
	for (size_t i=0; i<m_size; i++) m_data[i] -= arg.m_data[i];
}

void Vector::operator *= (double arg)
{
	for (size_t i=0; i<m_size; i++) m_data[i] *= arg;
}

void Vector::operator /= (double arg)
{
	for (size_t i=0; i<m_size; i++) m_data[i] /= arg;
}

Vector Vector::operator - () const
{
	Vector ret(m_size);
	for (size_t i=0; i<m_size; i++) ret.m_data[i] = -m_data[i];
	return ret;
}

Vector Vector::operator + (Vector const& rhs) const
{
	Vector ret = copy();;
	ret += rhs;
	return ret;
}

Vector Vector::operator - (Vector const& rhs) const
{
	Vector ret = copy();;
	ret -= rhs;
	return ret;
}

Vector Vector::operator * (double rhs) const
{
	Vector ret = copy();;
	ret *= rhs;
	return ret;
}

double Vector::operator * (Vector const& rhs) const
{
#ifdef CHECKRANGE
	if (rhs.size() != m_size) throw runtime_error("dimension mismatch");
#endif
	double ret = 0.0;
	for (size_t i=0; i<m_size; i++) ret += m_data[i] * rhs.m_data[i];
	return ret;
}

Vector Vector::operator / (double rhs) const
{
	Vector ret = copy();;
	ret /= rhs;
	return ret;
}

bool Vector::operator == (Vector const& rhs) const
{
#ifdef CHECKRANGE
	if (rhs.size() != m_size) throw runtime_error("dimension mismatch");
#endif
	for (size_t i=0; i<size(); i++) if (m_data[i] != rhs[i]) return false;
	return true;
}

double Vector::twonorm() const
{ return sqrt(twonorm2()); }

double Vector::twonorm2() const
{ return (*this) * (*this); }

Vector operator * (double lhs, Vector const& rhs)
{ return rhs * lhs; }

std::ostream& operator << (std::ostream& str, Vector const& vec)
{
	str << "[";
	if (! vec.empty())
	{
		str << vec[0];
		for (size_t i=1; i<vec.size(); i++) str << ", " << vec[i];
	}
	str << "]";
	return str;
}
