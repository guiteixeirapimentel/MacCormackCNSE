#pragma once
#include <vector>
#include <functional>
#include <cassert>

template<class vectorType, class floatType>
class SpaceField
{
public:
	SpaceField(size_t nx, size_t ny, floatType dx, floatType dy, vectorType valInic = static_cast<vectorType>(0.0))
		:
		cNX(nx),
		cNY(ny),
		cdx(dx),
		cdy(dy)
	{
		cData.resize(cNX*cNY, valInic);
	}
	SpaceField(SpaceField<vectorType, floatType>&& spaceField)
		:
		cNX(spaceField.GetNX()),
		cNY(spaceField.GetNY()),
		cdx(spaceField.Getdx()),
		cdy(spaceField.Getdy()),
		cData(std::move(spaceField.cData))
	{}
	SpaceField(const SpaceField<vectorType, floatType>& spaceField)
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

	inline floatType Getdx() const { return cdx; }
	inline floatType Getdy() const { return cdy; }

	inline vectorType& operator()(size_t i, size_t j) { return cData[i + (j*cNX)]; }
	inline vectorType& operator()(size_t ii) { return cData[ii]; }

	inline vectorType operator()(size_t i, size_t j) const { return cData[i + (j*cNX)]; }
	inline vectorType operator()(size_t ii) const { return cData[ii]; }

	// i -> deltax ; j -> deltay
	inline vectorType GetValue(size_t i, size_t j) const { return cData[i + (j*cNX)]; }
	inline vectorType GetValue(size_t ii) const { return cData[ii]; }
	inline vectorType SetValue(size_t i, size_t j, const vectorType& value) { return (cData[i + (j*cNX)] = value); }
	inline vectorType SetValue(size_t ii, const vectorType& value) { return (cData[ii] = value); }
	
	inline vectorType ForwardXDifferentiation(size_t i, size_t j) const { return (cData[i + 1 + (j*cNX)] - cData[i + (j*cNX)])/cdx; };
	inline vectorType BackwardXDifferentiation(size_t i, size_t j) const { return (cData[i + (j*cNX)] - cData[i - 1 + (j*cNX)]) / cdx; };

	inline vectorType ForwardYDifferentiation(size_t i, size_t j) const { return (cData[i + ((j+1)*cNX)] - cData[i + (j*cNX)]) / cdy; };
	inline vectorType BackwardYDifferentiation(size_t i, size_t j) const { return (cData[i + (j*cNX)] - cData[i + ((j-1)*cNX)]) / cdy; };

	inline vectorType CentralXDifferentiation(size_t i, size_t j) const { return (cData[i + 1 + (j*cNX)] - cData[i - 1 + (j*cNX)]) / cdx; }
	inline vectorType CentralYDifferentiation(size_t i, size_t j) const { return (cData[i + ((j + 1 )*cNX)] - cData[i + ((j - 1 )*cNX)]) / cdy; }

	inline vectorType ForwardXDifferentiation(size_t ii) const { return (cData[ii + 1] - cData[ii]) / cdx; };
	inline vectorType BackwardXDifferentiation(size_t ii) const { return (cData[ii] - cData[ii - 1]) / cdx; };

	inline vectorType ForwardYDifferentiation(size_t ii) const { return (cData[ii + 1*cNX] - cData[ii]) / cdy; };
	inline vectorType BackwardYDifferentiation(size_t ii) const { return (cData[ii] - cData[ii - 1*cNX]) / cdy; };

	inline vectorType CentralXDifferentiation(size_t ii) const { return (cData[ii + 1] - cData[ii - 1]) / cdx; }
	inline vectorType CentralYDifferentiation(size_t ii) const { return (cData[ii + 1*cNX] - cData[ii - 1*cNX]) / cdy; }

	inline void ForwardXDifferentiation(SpaceField<vectorType, floatType>& outSpaceField) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cData.size())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif

		for (size_t j = 0; j < cNY; j++)
		{
			for (size_t i = 1; i < cNX - 1; i++)
			{
				outSpaceField.SetValue(i, j, ForwardXDifferentiation(i, j));
			}
		}
	}

	inline void BackwardXDifferentiation(SpaceField<vectorType, floatType>& outSpaceField) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cData.size())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif

		for (size_t j = 0; j < cNY; j++)
		{
			for (size_t i = 1; i < cNX - 1; i++)
			{
				outSpaceField.SetValue(i, j, BackwardXDifferentiation(i, j));
			}
		}
	}

	inline void ForwardYDifferentiation(SpaceField<vectorType, floatType>& outSpaceField) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cData.size())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif

		for (size_t j = 1; j < cNY - 1; j++)
		{
			for (size_t i = 0; i < cNX; i++)
			{
				outSpaceField.SetValue(i, j, ForwardYDifferentiation(i, j));
			}
		}
	}

	inline void BackwardYDifferentiation(SpaceField<vectorType, floatType>& outSpaceField) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cData.size())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif

		for (size_t j = 1; j < cNY-1; j++)
		{
			for (size_t i = 0; i < cNX; i++)
			{
				outSpaceField.SetValue(i, j, BackwardYDifferentiation(i, j));
			}
		}
	}

	inline void CentralXDifferentiation(SpaceField<vectorType, floatType>& outSpaceField) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cData.size())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif

		for (size_t j = 0; j < cNY; j++)
		{
			for (size_t i = 1; i < cNX - 1; i++)
			{
				outSpaceField.SetValue(i, j, CentralXDifferentiation(i, j));
			}
		}
	}

	inline void CentralYDifferentiation(SpaceField<vectorType, floatType>& outSpaceField) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cData.size())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif

		for (size_t j = 1; j < cNY - 1; j++)
		{
			for (size_t i = 0; i < cNX; i++)
			{
				outSpaceField.SetValue(i, j, CentralYDifferentiation(i, j));
			}
		}
	}
	
	inline SpaceField& operator*=(floatType rhs)
	{
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] *= rhs;
		}

		return *this;
	}

	inline SpaceField& operator/=(floatType rhs)
	{
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] /= rhs;
		}

		return *this;
	}

	inline SpaceField& operator+=(floatType rhs)
	{
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] += rhs;
		}

		return *this;
	}

	inline SpaceField& operator-=(floatType rhs)
	{
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] -= rhs;
		}

		return *this;
	}

	inline SpaceField operator*(floatType rhs) const
	{
		SpaceField res(*this);
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] *= rhs;
		}

		return res;
	}

	inline SpaceField operator/(floatType rhs) const
	{
		SpaceField res = *this;
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] /= rhs;
		}

		return res;
	}

	inline SpaceField& operator+(floatType rhs) const
	{
		SpaceField res = *this;
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] += rhs;
		}

		return res;
	}

	inline SpaceField& operator-(floatType rhs)const
	{
		SpaceField res = *this;
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] -= rhs;
		}

		return res;
	}

	inline SpaceField operator-() const
	{
		SpaceField res = *this;
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] = -res.cData[ii];
		}

		return res;
	}
			
	inline SpaceField& operator+=(const SpaceField<vectorType, floatType>& rhs)
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

	inline SpaceField& operator-=(const SpaceField<vectorType, floatType>& rhs)
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

	inline SpaceField& operator*=(const SpaceField<vectorType, floatType>& rhs)
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

	inline SpaceField& operator/=(const SpaceField<vectorType, floatType>& rhs)
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

	inline SpaceField operator+(const SpaceField<vectorType, floatType>& rhs)const
	{
		SpaceField<vectorType, floatType> res = *this;
#ifdef _DEBUG
		if (rhs.GetNN() != cData.size())
		{
			throw "Campos de tamanhos diferentes!";
		}
#endif
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] = cData[ii] + rhs(ii);
		}
		return res;
	}

	inline SpaceField operator-(const SpaceField<vectorType, floatType>& rhs)const
	{
		SpaceField<vectorType, floatType> res = *this;
#ifdef _DEBUG
		if (rhs.GetNN() != cData.size())
		{
			throw "Campos de tamanhos diferentes!";
		}
#endif
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] = cData[ii] - rhs(ii);
		}
		return res;
	}

	inline SpaceField operator*(const SpaceField<vectorType, floatType>& rhs)const
	{
		SpaceField<vectorType, floatType> res = *this;
#ifdef _DEBUG
		if (rhs.GetNN() != cData.size())
		{
			throw "Campos de tamanhos diferentes!";
		}
#endif
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] = cData[ii] * rhs(ii);
		}
		return res;
	}

	inline SpaceField operator/(const SpaceField<vectorType, floatType>& rhs)const
	{
		SpaceField<vectorType, floatType> res = *this;
#ifdef _DEBUG
		if (rhs.GetNN() != cData.size())
		{
			throw "Campos de tamanhos diferentes!";
		}
#endif
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.cData[ii] = cData[ii] / rhs(ii);
		}
		return res;
	}

	inline SpaceField& operator()(const std::function<vectorType(vectorType)>& function)
	{
		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			cData[ii] = function(cData[ii]);
		}
		return *this;
	}

	inline SpaceField operator()(const std::function<vectorType(vectorType)>& function) const
	{
		SpaceField res = *this;

		for (size_t ii = 0; ii < cData.size(); ii++)
		{
			res.SetValue(ii) = function(cData[ii]);
		}
		return res;
	}

	
	inline SpaceField& operator=(const SpaceField& spaceField)
	{
#ifdef _DEBUG

		if (cNX != spaceField.GetNX() || cNY != spaceField.GetNY() ||
			cdx != spaceField.Getdx() || cdy != spaceField.Getdy())
		{
			throw "Campos de tamanhos diferentes!";
		}
#endif

		cData = spaceField.cData;

		return *this;
	}

private:
	const size_t cNX;
	const size_t cNY;
	const floatType cdx;
	const floatType cdy;
	std::vector<vectorType> cData;
};