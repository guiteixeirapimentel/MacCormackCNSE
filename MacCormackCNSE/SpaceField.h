#pragma once
#include <vector>
#include <functional>
#include <cassert>

template<class T>
class SpaceField
{
public:
	SpaceField(size_t nx, size_t ny, T dx, T dy, T valInic = static_cast<T>(0.0))
		:
		cNX(nx),
		cNY(ny),
		cdx(dx),
		cdy(dy)
	{
		cData.resize(cNX*cNY, valInic);
	}
	SpaceField(SpaceField<T>&& spaceField)
		:
		cNX(spaceField.GetNX()),
		cNY(spaceField.GetNY()),
		cdx(spaceField.Getdx()),
		cdy(spaceField.Getdy()),
		cData(std::move(spaceField.cData))
	{}
	SpaceField(const SpaceField<T>& spaceField)
		:
		cNX(spaceField.GetNX()),
		cNY(spaceField.GetNY()),
		cdx(spaceField.Getdx()),
		cdy(spaceField.Getdy()),
		cData(spaceField.cData)
	{

	}
	~SpaceField() {};

	inline size_t GetNX() const { return cNX; }
	inline size_t GetNY() const { return cNY; }
	inline size_t GetNN() const { return cData.size(); }

	inline T Getdx() const { return cdx; }
	inline T Getdy() const { return cdy; }

	inline T operator()(size_t i, size_t j) const { return cData[i + (j*cNX)]; }
	inline T operator()(size_t ii) const { return cData[ii]; }

	// i -> deltax ; j -> deltay
	inline T GetValue(size_t i, size_t j) const { return cData[i + (j*cNX)]; }
	inline T GetValue(size_t ii) const { return cData[ii]; }
	inline T SetValue(size_t i, size_t j, T value) { return (cData[i + (j*cNX)] = value); }
	inline T SetValue(size_t ii, T value) { return (cData[ii] = value); }
	
	inline T FowardXDiferenciation(size_t i, size_t j) const { return (cData[i + 1 + (j*cNX)] - cData[i + (j*cNX)])/cdx; };
	inline T BackwardXDiferenciation(size_t i, size_t j) const { return (cData[i + (j*cNX)] - cData[i - 1 + (j*cNX)]) / cdx; };

	inline T FowardYDiferenciation(size_t i, size_t j) const { return (cData[i + ((j+1)*cNX)] - cData[i + (j*cNX)]) / cdx; };
	inline T BackwardYDiferenciation(size_t i, size_t j) const { return (cData[i + (j*cNX)] - cData[i + ((j - 1)*cNX)]) / cdx; };

	inline void ForwardXDifferention(SpaceField<T>& outSpaceField) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cData.size())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif

		for (size_t j = 0; j < cNY; j++)
		{
			for (size_t i = 0; i < cNX; i++)
			{
				outSpaceField.SetValue(i, j, ForwardXDifferention(i, j));
			}
		}
	}

	inline void BackwardXDifferention(SpaceField<T>& outSpaceField) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cData.size())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif

		for (size_t j = 0; j < cNY; j++)
		{
			for (size_t i = 0; i < cNX; i++)
			{
				outSpaceField.SetValue(i, j, BackwardXDifferention(i, j));
			}
		}
	}

	inline void ForwardYDifferention(SpaceField<T>& outSpaceField) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cData.size())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif

		for (size_t j = 0; j < cNY; j++)
		{
			for (size_t i = 0; i < cNX; i++)
			{
				outSpaceField.SetValue(i, j, ForwardYDifferention(i, j));
			}
		}
	}

	inline void BackwardYDifferention(SpaceField<T>& outSpaceField) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cData.size())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif

		for (size_t j = 0; j < cNY; j++)
		{
			for (size_t i = 0; i < cNX; i++)
			{
				outSpaceField.SetValue(i, j, BackwardYDifferention(i, j));
			}
		}
	}
	
	inline SpaceField& operator*=(T rhs)
	{
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] *= rhs;
		}

		return *this;
	}

	inline SpaceField& operator/=(T rhs)
	{
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] /= rhs;
		}

		return *this;
	}

	inline SpaceField& operator+=(T rhs)
	{
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] += rhs;
		}

		return *this;
	}

	inline SpaceField& operator-=(T rhs)
	{
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] -= rhs;
		}

		return *this;
	}

	inline SpaceField operator*(T rhs) const
	{
		SpaceField res(*this);
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] *= rhs;
		}

		return res;
	}

	inline SpaceField operator/(T rhs) const
	{
		SpaceField res = *this;
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] /= rhs;
		}

		return res;
	}

	inline SpaceField& operator+(T rhs) const
	{
		SpaceField res = *this;
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] += rhs;
		}

		return res;
	}

	inline SpaceField& operator-(T rhs)const
	{
		SpaceField res = *this;
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] -= rhs;
		}

		return res;
	}
			
	inline SpaceField& operator+=(const SpaceField<T>& rhs)
	{
#ifdef _DEBUG
		if (rhs.GetNN() != cData.size())
		{
			throw "Campos de tamanhos diferentes!";
		}
#endif
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] += rhs[ii];
		}
		return *this;
	}

	inline SpaceField& operator-=(const SpaceField<T>& rhs)
	{
#ifdef _DEBUG
		if (rhs.GetNN() != cData.size())
		{
			throw "Campos de tamanhos diferentes!";
		}
#endif
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] -= rhs[ii];
		}
		return *this;
	}

	inline SpaceField& operator*=(const SpaceField<T>& rhs)
	{
#ifdef _DEBUG
		if (rhs.GetNN() != cData.size())
		{
			throw "Campos de tamanhos diferentes!";
		}
#endif
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] *= rhs[ii];
		}
		return *this;
	}

	inline SpaceField& operator/=(const SpaceField<T>& rhs)
	{
#ifdef _DEBUG
		if (rhs.GetNN() != cData.size())
		{
			throw "Campos de tamanhos diferentes!";
		}
#endif
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] /= rhs[ii];
		}
		return *this;
	}

	inline SpaceField& operator()(const std::function<T(T)>& function)
	{
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] = function(cData[ii]);
		}
		return *this;
	}

	inline SpaceField operator()(const std::function<T(T)>& function) const
	{
		SpaceField res = *this;

		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.SetValue(ii) = function(cData[ii]);
		}
		return res;
	}

	inline SpaceField& operator=(SpaceField&& spaceField)
	{
		cNX= spaceField.GetNX();

		cNY= spaceField.GetNY();

		cdx= spaceField.Getdx();

		cdy= spaceField.Getdy();

		cData(std::move(spaceField.cData));
	}
	inline SpaceField& operator=(const SpaceField& spaceField) const
	{
		cNX = spaceField.GetNX();

		cNY = spaceField.GetNY();

		cdx = spaceField.Getdx();

		cdy = spaceField.Getdy();

		cData(spaceField.cData);
	}

private:
	const size_t cNX;
	const size_t cNY;
	const T cdx;
	const T cdy;
	std::vector<T> cData;
};