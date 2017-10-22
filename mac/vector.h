
#pragma once


#include <vector>
#include <initializer_list>
#include <ostream>
#include <cstdlib>


// Numerical vector type with either memory owner or view semantics.
class Vector
{
public:
	// construction / destruction
	Vector();                                   // "empty" vector (dimension = 0)
	explicit Vector(std::size_t dim);           // vector of given dimension, uninitialized
	Vector(double scalar);                      // one-dimensional vector with given value
	Vector(std::size_t dim, double value);      // vector with all components initialized to value
	Vector(std::size_t dim, const double* data);// construct vector from linear data array
	Vector(std::vector<double> const& other);   // copy of standard vector (not a view!)
	Vector(double* data, std::size_t dim);      // view on the given data, e.g., for sub-vectors
	Vector(Vector const& other);                // copy constructor (sloppy, drops const-ness of views)
	Vector(std::initializer_list<double> l);    // initializer list constructor
	~Vector();

	// assignment
	Vector& operator = (Vector const& rhs);
	void operator = (double value)
	{ for (std::size_t i=0; i<m_size; i++) m_data[i] = value; }
	void concat(double arg);
	void concat(Vector const& arg);

	Vector copy() const;

	// iterators
	double* begin()
	{ return m_data; }
	const double* begin() const
	{ return m_data; }
	double* end()
	{ return m_data + m_size; }
	const double* end() const
	{ return m_data + m_size; }

	// size, element access, sub-array (slice) access
	bool empty() const
	{ return (m_size == 0); }
	std::size_t size() const
	{ return m_size; }
	double& operator [] (std::size_t index);
	double const& operator [] (std::size_t index) const;
	double& operator () (std::size_t index);
	double const& operator () (std::size_t index) const;
	double* data()
	{ return m_data; }
	const double* data() const
	{ return m_data; }
	Vector sub(std::size_t start, std::size_t end, bool view = true);
	bool isView() const
	{ return m_view; }

	// arithmetic operations
	void operator += (Vector const& arg);
	void operator -= (Vector const& arg);
	void operator *= (double arg);
	void operator /= (double arg);
	Vector operator - () const;
	Vector operator + (Vector const& rhs) const;
	Vector operator - (Vector const& rhs) const;
	Vector operator * (double rhs) const;
	double operator * (Vector const& rhs) const;
	Vector operator / (double rhs) const;

	// standard comparison for equality
	bool operator == (Vector const& rhs) const;
	bool operator != (Vector const& rhs) const
	{ return ! (operator == (rhs)); }

	// Euclidean norm and squared norm
	double twonorm() const;
	double twonorm2() const;

private:
	std::size_t m_size;
	double* m_data;
	bool m_view;
};

Vector operator * (double lhs, Vector const& rhs);
std::ostream& operator << (std::ostream& str, Vector const& vec);
