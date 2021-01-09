#include <iostream>

#include "SpaceField.h"
#include "FixedColumnVector.h"

typedef float precisao;
typedef SpaceField<FixedColumnVector<precisao>, precisao> SpaceFieldF;

void ExibirSpaceField(const SpaceFieldF& spf);

int main()
{
	// valores padrao
	constexpr precisao P0 = static_cast<precisao>(101325.0);
	constexpr precisao a0 = static_cast<precisao>(340.28);
	constexpr precisao T0 = static_cast<precisao>(288.16);
	
	constexpr precisao gama = static_cast<precisao>(1.4);
	constexpr precisao Pr = static_cast<precisao>(0.71);
	constexpr precisao mu0 = static_cast<precisao>(1.7894e-5);
	constexpr precisao R = static_cast<precisao>(287.0);
	constexpr precisao Cp = static_cast<precisao>(1.0048e3);
	constexpr precisao Cv = static_cast<precisao>(Cp/gama);

	constexpr precisao beta = static_cast<precisao>(0.0);
	
	// definicao de funcoes lambda

	// sutherland law
	auto visc = [&](precisao T)
	{
		return mu0 * pow(T / T0, static_cast<precisao>(1.5)) * (T0 + static_cast<precisao>(110.0)) / (T + static_cast<precisao>(110.0));
	};

	auto condTermica = [&](precisao visc)
	{
		return visc * Cp / Pr;
	};

	auto segVisc = [&](precisao visc)
	{
		return beta - static_cast<precisao>(2.0 / 3.0)*visc;
	};



	SpaceFieldF spaceField(2, 2, 0.1f, 0.1f, 0.0f);

	spaceField += 2.0f;
	spaceField *= 2.0f;

	SpaceFieldF a = spaceField*2;
	
	ExibirSpaceField(spaceField);
	ExibirSpaceField(a);

	return 0;
}

void ExibirSpaceField(const SpaceFieldF& spf)
{
	const size_t nn = spf.GetNN();
	for (size_t ii = 0; ii < nn; ii++)
	{
		std::cout << spf(ii).cValues[0] << " ; ";
	}
}
