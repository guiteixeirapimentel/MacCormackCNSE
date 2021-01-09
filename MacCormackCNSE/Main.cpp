#include <iostream>

#include "SpaceField.h"

typedef float precisao;

void ExibirSpaceField(const SpaceField<precisao>& spf);

int main()
{
	SpaceField<precisao> spaceField(2, 2, 0.1f, 0.1f);

	spaceField += 2.0f;
	spaceField *= 2.0f;

	SpaceField<precisao> a = spaceField*2;
	
	ExibirSpaceField(spaceField);
	ExibirSpaceField(a);


	return 0;
}

void ExibirSpaceField(const SpaceField<precisao>& spf)
{
	const size_t nn = spf.GetNN();
	for (size_t ii = 0; ii < nn; ii++)
	{
		std::cout << spf(ii) << " ; ";
	}
}
