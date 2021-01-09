#pragma once
#include "SpaceField.h"

template <class T>
class TimeField
{
public:
	TimeField(const SpaceField<T>& spaceFieldInitial)
	{
		cFieldn = spaceFieldInitial;
		cFieldnp1 = spaceFieldInitial;
	}

	~TimeField() {};

	inline T ForwardTDifferentiation(size_t i, size_t j, T dt) const { return (cFieldnpm(i, j) - cFieldn(i, j)) / dt; }
	inline T BackwardTDifferentiation(size_t i, size_t j, T dt) const { return (cFieldnpm(i, j) - cFieldn(i, j)) / dt; }

	inline void ForwardTDifferentiation(SpaceField& outSpaceField, T dt) const 
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cFieldn.GetNN())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif
		const size_t nx = cFieldn.GetNX();
		const size_t ny = cFieldn.GetNY();

		for (size_t j = 0; j < ny; j++)
		{
			for (size_t i = 0; i < nx; i++)
			{
				outSpaceField.SetValue(i, j, ForwardTDifferentiation(i, j, dt));
			}
		}
	}
	inline void BackwardTDifferentiation(SpaceField& outSpaceField, T dt) const
	{
#ifdef _DEBUG
		if (outSpaceField.GetNN() != cFieldn.GetNN())
		{
			throw "Campo usado para armazenar diferenciacao com dimensoes erradas!";
		}
#endif
		const size_t nx = cFieldn.GetNX();
		const size_t ny = cFieldn.GetNY();

		for (size_t j = 0; j < ny; j++)
		{
			for (size_t i = 0; i < nx; i++)
			{
				outSpaceField.SetValue(i, j, BackwardTDifferentiation(i, j, dt));
			}
		}
	}

	inline SpaceField GetFieldn() const { return cFieldn; }
	inline SpaceField GetFieldnp() const { return cFieldnp; }

private:
	SpaceField cFieldn;
	SpaceField cFieldnpm;
};