
#pragma once


#include "vector.h"

#include <vector>
#include <initializer_list>
#include <cassert>
#include <cmath>


// Simple numerical matrix class.
class Matrix
{
public:
	typedef std::vector<double>::iterator iterator;
	typedef std::vector<double>::const_iterator const_iterator;

	Matrix()
	: m_rows(0)
	, m_cols(0)
	, m_data(0)
	{ }

	Matrix(std::size_t rows, std::size_t cols, double value = 0.0)
	: m_rows(rows)
	, m_cols(cols)
	, m_data(rows * cols, value)
	{ }

	Matrix(std::size_t rows, std::size_t cols, const double* data)
	: m_rows(rows)
	, m_cols(cols)
	, m_data(data, data + rows * cols)
	{ }

	Matrix(std::size_t rows, std::size_t cols, std::initializer_list<double> l)
	: m_rows(rows)
	, m_cols(cols)
	, m_data(l)
	{ assert(m_data.size() == m_rows * m_cols); }

	void operator = (double value)
	{ for (std::size_t i=0; i<m_data.size(); i++) m_data[i] = value; }

	void operator = (Matrix const& other)
	{
		m_rows = other.m_rows;
		m_cols = other.m_cols;
		m_data = other.m_data;
	}

	inline std::size_t rows() const
	{ return m_rows; }
	inline std::size_t cols() const
	{ return m_cols; }
	inline std::size_t size1() const
	{ return rows(); }
	inline std::size_t size2() const
	{ return cols(); }
	inline bool empty() const
	{ return (size() == 0); }
	inline std::size_t size() const
	{ return (rows() * cols()); }
	inline double& operator () (std::size_t y, std::size_t x)
	{
		assert(y < m_rows && x < m_cols);
		return m_data[x + m_cols * y];
	}
	inline double const& operator () (std::size_t y, std::size_t x) const
	{
		assert(y < m_rows && x < m_cols);
		return m_data[x + m_cols * y];
	}

	inline double* data()
	{ return &m_data[0]; }
	inline const double* data() const
	{ return &m_data[0]; }

	inline iterator begin()
	{ return m_data.begin(); }
	inline const_iterator begin() const
	{ return m_data.begin(); }
	inline iterator end()
	{ return m_data.end(); }
	inline const_iterator end() const
	{ return m_data.end(); }

	bool finite() const
	{
		for (std::size_t i=0; i<size(); i++) if (! std::isfinite(m_data[i])) return false;
		return true;
	}

	inline void operator += (Matrix const& other)
	{
		assert(m_rows == other.rows() && m_cols == other.cols());
		std::size_t i, ic = m_rows * m_cols;
		for (i=0; i<ic; i++) m_data[i] += other.m_data[i];
	}
	inline void operator += (double scalar)
	{
		std::size_t i, ic = m_rows * m_cols;
		for (i=0; i<ic; i++) m_data[i] += scalar;
	}
	inline void operator -= (Matrix const& other)
	{
		assert(m_rows == other.rows() && m_cols == other.cols());
		std::size_t i, ic = m_rows * m_cols;
		for (i=0; i<ic; i++) m_data[i] -= other.m_data[i];
	}
	inline void operator -= (double scalar)
	{
		std::size_t i, ic = m_rows * m_cols;
		for (i=0; i<ic; i++) m_data[i] -= scalar;
	}
	inline void operator *= (double scalar)
	{
		std::size_t i, ic = m_rows * m_cols;
		for (i=0; i<ic; i++) m_data[i] *= scalar;
	}
	inline void operator /= (double scalar)
	{
		std::size_t i, ic = m_rows * m_cols;
		for (i=0; i<ic; i++) m_data[i] /= scalar;
	}
	inline Matrix operator - () const
	{
		std::size_t i, ic = size();
		Matrix ret(m_rows, m_cols);
		for (i=0; i<ic; i++) ret.m_data[i] = -m_data[i];
		return ret;
	}
	inline Matrix operator + (Matrix const& other) const
	{
		assert(m_rows == other.rows() && m_cols == other.cols());
		Matrix ret(*this);
		ret += other;
		return ret;
	}
	inline Matrix operator - (Matrix const& other) const
	{
		assert(m_rows == other.rows() && m_cols == other.cols());
		Matrix ret(*this);
		ret -= other;
		return ret;
	}
	inline Matrix operator * (double scalar) const
	{
		Matrix ret(*this);
		ret *= scalar;
		return ret;
	}
	inline Matrix operator / (double scalar) const
	{
		Matrix ret(*this);
		ret /= scalar;
		return ret;
	}
	inline Matrix operator * (Matrix const& other) const
	{
		assert(m_cols == other.rows());
		Matrix ret(m_rows, other.cols());
		std::size_t i, j, k;
		for (i=0; i<m_rows; i++)
		{
			for (j=0; j<other.cols(); j++)
			{
				double value = 0.0;
				for (k=0; k<m_cols; k++)
				{
					value += (*this)(i, k) * other(k, j);
				}
				ret(i, j) = value;
			}
		}
		return ret;
	}
	inline Vector operator * (Vector const& other) const
	{
		assert(m_cols == other.size());
		Vector ret(m_rows);
		std::size_t i, j;
		for (i=0; i<m_rows; i++)
		{
			double value = 0.0;
			for (j=0; j<m_cols; j++)
			{
				value += (*this)(i, j) * other[j];
			}
			ret[i] = value;
		}
		return ret;
	}

	inline void operator *= (Matrix const& other)
	{ *this = Matrix::operator *((Matrix const&)other); }

	inline bool operator == (Matrix const& other) const
	{
		if (other.m_cols != m_cols || other.m_rows != m_rows) return false;
		std::size_t sz = m_rows * m_cols;
		for (std::size_t i=0; i<sz; i++) if (m_data[i] != other.m_data[i]) return false;
		return true;
	}

	inline bool operator != (Matrix const& other) const
	{ return ! (operator == (other)); }

	// resize and set content to zeros
	inline void resize(std::size_t rows, std::size_t cols)
	{
		m_rows = rows;
		m_cols = cols;
		std::size_t i, ic = rows * cols;
		m_data.resize(ic);
		for (i=0; i<ic; i++) m_data[i] = 0.0;
	}

	inline void fill(double scalar)
	{
		std::fill(m_data.begin(), m_data.end(), scalar);
	}

	inline Matrix transpose() const
	{
		Matrix ret(m_cols, m_rows);
		std::size_t i, j;
		for (i=0; i<m_cols; i++)
			for (j=0; j<m_rows; j++)
				ret(i, j) = (*this)(j, i);
		return ret;
	}

	// non-const rows are views
	inline Vector row(std::size_t y)
	{ return Vector(&m_data[m_cols * y], m_cols); }

	// const rows are copies
	inline Vector row(std::size_t y) const
	{
		Vector ret(m_cols);
		for (std::size_t x=0; x<m_cols; x++) ret(x) = (*this)(y, x);
		return ret;
	}

	// columns are copies
	inline Vector col(std::size_t x) const
	{
		Vector ret(m_rows);
		for (std::size_t y=0; y<m_rows; y++) ret(y) = (*this)(y, x);
		return ret;
	}

	inline Matrix sub(std::size_t yBegin, std::size_t yNum, std::size_t xBegin, std::size_t xNum) const
	{
		Matrix ret(yNum, xNum);
		std::size_t i, j;
		for (i=0; i<yNum; i++)
			for (j=0; j<xNum; j++)
				ret(i, j) = (*this)(yBegin + i, xBegin + j);
		return ret;
	}

	// return I
	static Matrix identity(std::size_t dim)
	{
		Matrix ret(dim, dim);
		for (std::size_t i=0; i<dim; i++) ret(i, i) = 1.0;
		return ret;
	}

	// return v * w^T
	static Matrix outerProduct(Vector const& v, Vector const& w)
	{
		std::size_t i, j;
		Matrix ret(v.size(), w.size());
		for (i=0; i<v.size(); i++)
		{
			for (j=0; j<w.size(); j++)
			{
				ret(i, j) = v[i] * w[j];
			}
		}
		return ret;
	}

	// return v * v^T
	static Matrix outerProduct(Vector const& v)
	{
		std::size_t i, j, dim = v.size();
		Matrix ret(dim, dim);
		for (i=0; i<dim; i++)
		{
			for (j=0; j<i; j++)
			{
				double p = v[i] * v[j];
				ret(i, j) = p;
				ret(j, i) = p;
			}
			ret(i, i) = v[i] * v[i];
		}
		return ret;
	}

	static Matrix diag(Vector const& vec)
	{
		std::size_t i, ic = vec.size();
		Matrix ret(ic, ic);
		for (i=0; i<ic; i++) ret(i, i) = vec(i);
		return ret;
	}

	// trace
	inline double tr() const
	{
		assert(m_rows == m_cols);
		double ret = 0.0;
		for (std::size_t i=0; i<m_rows; i++) ret += (*this)(i, i);
		return ret;
	}

	inline double sum() const
	{
		double ret = 0.0;
		std::size_t i, ic = m_data.size();
		for (i=0; i<ic; i++) ret += m_data[i];
		return ret;
	}

	inline double onenorm() const
	{
		double ret = 0.0;
		const_iterator it = begin();
		for (; it != end(); it++) ret += std::abs(*it);
		return ret;
	}

	inline double twonorm2() const
	{
		double ret = 0.0;
		const_iterator it = begin();
		for (; it != end(); it++) ret += (*it) * (*it);
		return ret;
	}

	inline double twonorm() const
	{ return std::sqrt(twonorm2()); }

	inline double maxval() const
	{
		std::size_t i, ic = m_data.size();
		if (ic == 0) return -1e100;
		double ret = m_data[0];
		for (i=1; i<ic; i++)
		{
			double d = m_data[i];
			if (d > ret) ret = d;
		}
		return ret;
	}
	inline double minval() const
	{
		std::size_t i, ic = m_data.size();
		if (ic == 0) return +1e100;
		double ret = m_data[0];
		for (i=1; i<ic; i++)
		{
			double d = m_data[i];
			if (d < ret) ret = d;
		}
		return ret;
	}

protected:
	std::size_t m_rows;
	std::size_t m_cols;
	std::vector<double> m_data;
};


Matrix operator * (double lhs, Matrix const& rhs);
std::ostream& operator << (std::ostream& str, Matrix const& mat);
