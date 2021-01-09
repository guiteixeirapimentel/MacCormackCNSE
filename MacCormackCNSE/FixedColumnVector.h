#pragma once

template<class T>
class FixedColumnVector
{
public:
	FixedColumnVector() {};
	FixedColumnVector(T valInic) { for (size_t i = 0; i < cNColumns; i++) { cValues[i] = valInic; } }
	FixedColumnVector(const FixedColumnVector& fxcv) { for (size_t i = 0; i < cNColumns; i++) { cValues[i] = fxcv.cValues[i]; } }
	~FixedColumnVector() {};


	inline FixedColumnVector& operator+=(const FixedColumnVector& rhs)
	{
		for (size_t i = 0; i < cNColumns; i++)
		{
			cValues[i] += rhs.cValues[i];
		}

		return *this;
	}
	inline FixedColumnVector& operator-=(const FixedColumnVector& rhs)
	{
		for (size_t i = 0; i < cNColumns; i++)
		{
			cValues[i] -= rhs.cValues[i];
		}

		return *this;
	}
	inline FixedColumnVector& operator*=(const FixedColumnVector& rhs)
	{
		for (size_t i = 0; i < cNColumns; i++)
		{
			cValues[i] *= rhs.cValues[i];
		}

		return *this;
	}
	inline FixedColumnVector& operator/=(const FixedColumnVector& rhs)
	{
		for (size_t i = 0; i < cNColumns; i++)
		{
			cValues[i] /= rhs.cValues[i];
		}

		return *this;
	}

	inline FixedColumnVector operator+(T rhs) const 
	{
		FixedColumnVector res(*this);

		for (size_t i = 0; i < cNColumns; i++)
		{
			res.cValues[i] += rhs;
		}

		return res;
	}
	inline FixedColumnVector operator-(T rhs) const
	{
		FixedColumnVector res(*this);

		for (size_t i = 0; i < cNColumns; i++)
		{
			res.cValues[i] -= rhs;
		}

		return res;
	}
	inline FixedColumnVector operator*(T rhs) const
	{
		FixedColumnVector res(*this);

		for (size_t i = 0; i < cNColumns; i++)
		{
			res.cValues[i] *= rhs;
		}

		return res;
	}
	inline FixedColumnVector operator/(T rhs) const
	{
		FixedColumnVector res(*this);

		for (size_t i = 0; i < cNColumns; i++)
		{
			res.cValues[i] /= rhs;
		}

		return res;
	}

	inline FixedColumnVector& operator+=(T rhs)
	{
		for (size_t i = 0; i < cNColumns; i++)
		{
			cValues[i] += rhs;
		}

		return *this;
	}
	inline FixedColumnVector& operator-=(T rhs)
	{
		for (size_t i = 0; i < cNColumns; i++)
		{
			cValues[i] -= rhs;
		}

		return *this;
	}
	inline FixedColumnVector& operator*=(T rhs)
	{
		for (size_t i = 0; i < cNColumns; i++)
		{
			cValues[i] *= rhs;
		}

		return *this;
	}
	inline FixedColumnVector& operator/=(T rhs)
	{
		for (size_t i = 0; i < cNColumns; i++)
		{
			cValues[i] /= rhs;
		}

		return *this;
	}

	inline FixedColumnVector operator+(const FixedColumnVector& rhs) const
	{
		FixedColumnVector res(*this);

		for (size_t i = 0; i < cNColumns; i++)
		{
			res.cValues[i] += rhs.cValues[i];
		}

		return res;
	}
	inline FixedColumnVector operator-(const FixedColumnVector& rhs) const
	{
		FixedColumnVector res(*this);

		for (size_t i = 0; i < cNColumns; i++)
		{
			res.cValues[i] -= rhs.cValues[i];
		}

		return res;
	}
	inline FixedColumnVector operator*(const FixedColumnVector& rhs) const
	{
		FixedColumnVector res(*this);

		for (size_t i = 0; i < cNColumns; i++)
		{
			res.cValues[i] *= rhs.cValues[i];
		}

		return res;
	}
	inline FixedColumnVector operator/(const FixedColumnVector& rhs) const
	{
		FixedColumnVector res(*this);

		for (size_t i = 0; i < cNColumns; i++)
		{
			res.cValues[i] /= rhs.cValues[i];
		}

		return res;
	}

public:
	constexpr static size_t cNColumns = 4;
	T cValues[cNColumns];
};