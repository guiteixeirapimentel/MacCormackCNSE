#include <iostream>
#include <fstream>
#include <string>

#include "SpaceField.h"
#include "FixedColumnVector.h"

typedef float precisao;
typedef SpaceField<FixedColumnVector<precisao>, precisao> SpaceFieldVector;
typedef SpaceField<precisao, precisao> SpaceFieldScalar;

void ExibirSpaceFieldScalar(const SpaceFieldScalar& spf);

int main()
{
	// questoes de tempo de simulacao
	constexpr precisao tmax = static_cast<precisao>(15.0);
	
	constexpr precisao KfudgeFactor = static_cast<precisao>(0.5); // pagina 457 livro anderson

	// valores da malha
	constexpr size_t NX = 70;
	constexpr size_t NY = 70;
	constexpr size_t NN = NX * NY;

	constexpr precisao dx = static_cast<precisao>(1.45e-7);//0.0000209);
	constexpr precisao dy = static_cast<precisao>(1.19e-7);//0.0000209);

	constexpr precisao L = dx * NX;

	// cond. contorno
	constexpr precisao u0 = static_cast<precisao>(340.28*4.0);//(1.0);
	constexpr precisao v0 = static_cast<precisao>(0.0);
	
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

	constexpr precisao rho0 = static_cast<precisao>(P0 / (R*T0));

	constexpr precisao beta = static_cast<precisao>(0.0);

	constexpr precisao lambda0 = static_cast<precisao>(beta - precisao(2.0 / 3.0)*beta);

	// valores globais escoamento

	constexpr precisao Re0 = rho0 * u0*L / mu0;
	
	// definicao de funcoes lambda

	// sutherland law
	auto visc = [&](precisao T)
	{
		const precisao res = mu0 * pow(T / T0, static_cast<precisao>(1.5)) * (T0 + static_cast<precisao>(110.0)) / (T + static_cast<precisao>(110.0));
		
		assert(!isnan(res));

		return res;
	};

	auto condTermica = [&](precisao visc)
	{
		return visc * Cp / Pr;
	};

	auto segVisc = [&](precisao visc)
	{
		return beta - static_cast<precisao>(2.0 / 3.0)*visc;
	};

	auto e = [&](precisao T)
	{
		return T * Cv;
	};

	auto V = [](precisao u, precisao v)
	{
		return sqrt(u*u + v*v);
	};

	auto V2 = [](precisao u, precisao v)
	{
		return u * u + v * v;
	};

	auto deltatCFL = [&](precisao u, precisao v, precisao a, precisao vlinha)
	{
		const precisao c1 = static_cast<precisao>(1.0) / (dx*dx) + static_cast<precisao>(1.0) / (dy*dy);
		return static_cast<precisao>(static_cast<precisao>(1.0) / (abs(u) / dx + abs(v) / dy + a * sqrt(c1) + static_cast<precisao>(2.0)*vlinha*c1));
	};

	auto vlinhaCFL = [&](precisao mu, precisao rho)
	{
		return static_cast<precisao>(4.0 / 3.0)*(gama*mu / Pr) / rho;
	};

	auto velSom = [&](precisao T)
	{
		return sqrt(R*gama*T);
	};

	SpaceFieldScalar un(NX, NY, dx, dy, u0);
	SpaceFieldScalar vn(NX, NY, dx, dy, v0);
	SpaceFieldScalar Pn(NX, NY, dx, dy, P0);
	SpaceFieldScalar Tn(NX, NY, dx, dy, T0);
	SpaceFieldScalar rhon(NX, NY, dx, dy, rho0);
	SpaceFieldScalar Etn(NX, NY, dx, dy, rho0*(e(T0) + V2(u0, v0) / 2));
	SpaceFieldScalar mun(NX, NY, dx, dy, mu0);
	SpaceFieldScalar kn(NX, NY, dx, dy, condTermica(mu0));
	SpaceFieldScalar lambdan(NX, NY, dx, dy, lambda0);
	SpaceFieldScalar an(NX, NY, dx, dy, a0);

	SpaceFieldScalar tauxxn(NX, NY, dx, dy);
	SpaceFieldScalar tauyyn(NX, NY, dx, dy);
	SpaceFieldScalar tauxyn(NX, NY, dx, dy);

	SpaceFieldScalar qxn(NX, NY, dx, dy);
	SpaceFieldScalar qyn(NX, NY, dx, dy);

	SpaceFieldScalar buf1(NX, NY, dx, dy);
	SpaceFieldScalar buf2(NX, NY, dx, dy);
	SpaceFieldScalar buf3(NX, NY, dx, dy);
	SpaceFieldScalar buf4(NX, NY, dx, dy);

	/*SpaceFieldScalar unp = un;
	SpaceFieldScalar vnp = vn;
	SpaceFieldScalar Pnp = Pn;
	SpaceFieldScalar Tnp = Tn;
	SpaceFieldScalar rhonp = rhon;
	SpaceFieldScalar Etnp = Etn;
	SpaceFieldScalar munp = mun;
	SpaceFieldScalar knp = kn;
	SpaceFieldScalar lambdanp = lambdan;*/

	SpaceFieldVector Un(NX, NY, dx, dy);
	SpaceFieldVector Unp(NX, NY, dx, dy);
	SpaceFieldVector E(NX, NY, dx, dy);
	SpaceFieldVector F(NX, NY, dx, dy);

	SpaceFieldVector dUdt1(NX, NY, dx, dy);
	SpaceFieldVector dUdt2(NX, NY, dx, dy);

	SpaceFieldVector bufV1(NX, NY, dx, dy);
	SpaceFieldVector bufV2(NX, NY, dx, dy);

	precisao t = static_cast<precisao>(0);

	auto garanteCondsContorno = [&](SpaceFieldScalar& u, SpaceFieldScalar& v, 
		SpaceFieldScalar& p, SpaceFieldScalar& T)
	{
		// Case 1
		u(0, (NY - 1)) = static_cast<precisao>(0);
		v(0, (NY - 1)) = static_cast<precisao>(0);
		p(0, (NY - 1)) = static_cast<precisao>(P0);
		T(0, (NY - 1)) = static_cast<precisao>(T0);

		// Case 2
		for (size_t j = 0; j < NY-1; j++)
		{
			u(0, j) = static_cast<precisao>(u0);
			v(0, j) = static_cast<precisao>(0);
			p(0, j) = static_cast<precisao>(P0);
			T(0, j) = static_cast<precisao>(T0);
		}
		for (size_t i = 0; i < NX; i++)
		{
			u(i, 0) = static_cast<precisao>(u0);
			v(i, 0) = static_cast<precisao>(0);
			p(i, 0) = static_cast<precisao>(P0);
			T(i, 0) = static_cast<precisao>(T0);
		}

		// Case 3
		for (size_t i = 1; i < NX; i++)
		{
			u(i, NY-1) = static_cast<precisao>(0);
			v(i, NY-1) = static_cast<precisao>(0);
			p(i, NY-1) = static_cast<precisao>(2) * p(i, NY-2) - p(i, NY - 3);
			T(i, NY-1) = static_cast<precisao>(T0);
		}

		// Case 4
		for (size_t j = 1; j < NY - 1; j++)
		{
			u(NX - 1, j) = static_cast<precisao>(2) * u(NX-2, j) - u(NX-3, j);
			v(NX - 1, j) = static_cast<precisao>(2) * v(NX - 2, j) - v(NX - 3, j);
			p(NX - 1, j) = static_cast<precisao>(2) * p(NX - 2, j) - p(NX - 3, j);
			T(NX - 1, j) = static_cast<precisao>(2) * T(NX - 2, j) - T(NX - 3, j);
		}
	};

	auto saveToCSV = [](const SpaceFieldScalar& field, const std::string& nome)
	{
		std::ofstream arq(nome);

		for (size_t j = 0; j < field.GetNY(); j++)
		{
			for (size_t i = 0; i < field.GetNX(); i++)
			{
				arq << field(i, j) << ';';
			}

			arq << '\n';
		}

		arq.close();
	};
	
	garanteCondsContorno(un, vn, Pn, Tn);

	while (t < tmax)
	{
		/// MACCORMACK

		// Determinar passo no tempo (método explicito)

		precisao vlinhamax = static_cast<precisao>(-10e32);

		for (size_t ii = 0; ii < NN; ii++)
		{
			vlinhamax = std::fmax(vlinhaCFL(mun(ii), rhon(ii)), vlinhamax);
		}

		precisao dtmin = static_cast<precisao>(10e32);

		for (size_t ii = 0; ii < NN; ii++)
		{
			dtmin = std::fmin(dtmin, deltatCFL(un(ii), vn(ii), an(ii), vlinhamax));
		}

		dtmin *= KfudgeFactor;
		
		std::cout << "Passo no tempo: " << dtmin << std::endl;		

		// Calcula valores do campo de vetores (U, E, F, etc) (como passo predizido, calcula as
		// derivadas lembrando da regra de trocar direçao no caso de mesma derivada espacial
		// e usar derivada central no caso de diferente
		// primeiramente usa-se derivadas para trás (backward differentiation)

		for (size_t ii = 0; ii < NN; ii++)
		{
			Un(ii).cValues[0] = rhon(ii);
			Un(ii).cValues[1] = rhon(ii)*un(ii);
			Un(ii).cValues[2] = rhon(ii)*vn(ii);
			Un(ii).cValues[3] = Etn(ii);
		}

		un.BackwardXDifferentiation(buf1);
		vn.CentralYDifferentiation(buf2);
		Tn.BackwardXDifferentiation(buf3);
			   		 
		for (size_t ii = 0; ii < NN; ii++)
		{
			tauxxn(ii) = lambdan(ii)*(buf1(ii) + buf2(ii)) + static_cast<precisao>(2) * mun(ii)*buf1(ii);
			qxn(ii) = -kn(ii) * buf3(ii);
		}

		un.CentralYDifferentiation(buf1);
		vn.BackwardXDifferentiation(buf2);

		for (size_t ii = 0; ii < NN; ii++)
		{
			tauxyn(ii) = mun(ii)*(buf1(ii) + buf2(ii));
		}

		saveToCSV(tauxyn, "tauxy.csv");

		for (size_t ii = 0; ii < NN; ii++)
		{				
			E(ii).cValues[0] = rhon(ii)*un(ii);
			E(ii).cValues[1] = rhon(ii)*un(ii)*un(ii) + Pn(ii) - tauxxn(ii);
			E(ii).cValues[2] = rhon(ii)*un(ii)*vn(ii) - tauxyn(ii);
			E(ii).cValues[3] = (Etn(ii) + Pn(ii)) * un(ii) - un(ii)*tauxxn(ii) - vn(ii)*tauxyn(ii) + qxn(ii);
		}

		un.CentralXDifferentiation(buf1);
		vn.BackwardYDifferentiation(buf2);
		Tn.BackwardYDifferentiation(buf3);

		for (size_t ii = 0; ii < NN; ii++)
		{
			tauyyn(ii) = lambdan(ii)*(buf1(ii) + buf2(ii)) + static_cast<precisao>(2) * mun(ii)*buf2(ii);
			qyn(ii) = -kn(ii) * buf3(ii);
		}

		un.BackwardYDifferentiation(buf1);
		vn.CentralXDifferentiation(buf2);

		for (size_t ii = 0; ii < NN; ii++)
		{
			tauxyn(ii) = mun(ii)*(buf1(ii) + buf2(ii));
		}

		for (size_t ii = 0; ii < NN; ii++)
		{
			F(ii).cValues[0] = rhon(ii)*vn(ii);
			F(ii).cValues[1] = rhon(ii)*un(ii)*vn(ii) - tauxyn(ii);
			F(ii).cValues[2] = rhon(ii)*vn(ii)*vn(ii) + Pn(ii) - tauyyn(ii);
			F(ii).cValues[3] = (Etn(ii) + Pn(ii)) * vn(ii) - un(ii)*tauxyn(ii) - vn(ii)*tauyyn(ii) + qyn(ii);
		}


		// Calcula dU/dt¹ predizido (usando diferenças para frente nos termos explicitos, dE/dx dF/dy, e trocado
		// quando a derivada está dentro do termo (cisalhamento) (lembrar de usar diferenças centrais quando
		// a derivada for na outra variavel).

		E.ForwardXDifferentiation(bufV1);
		F.ForwardYDifferentiation(bufV2);

		dUdt1 = -bufV1 - bufV2;
			   
		// Avançar U no tempo utilizando esses dU/dt
		Unp = Un + dUdt1 * dtmin;

		for (size_t ii = 0; ii < NN; ii++)
		{
			rhon(ii) = Unp(ii).cValues[0];
			un(ii) = Unp(ii).cValues[1] / Unp(ii).cValues[0];
			vn(ii) = Unp(ii).cValues[2] / Unp(ii).cValues[0];
			Etn(ii) = Unp(ii).cValues[3];

			Tn(ii) = (Etn(ii) / rhon(ii) - (un(ii)*un(ii) + vn(ii)*vn(ii)) / static_cast<precisao>(2.0))/Cv;
			Pn(ii) = rhon(ii)*R * Tn(ii);
			mun(ii) = visc(Tn(ii));
			kn(ii) = condTermica(mun(ii));
			lambdan(ii) = segVisc(mun(ii));
		}

		garanteCondsContorno(un, vn, Pn, Tn);

		// Calcula valores do campo de vetores (U, E, F, etc) (como passo corrigido, calcula as
		// derivadas lembrando da regra de trocar direçao no caso de mesma derivada espacial
		// e usar derivada central no caso de diferente

		un.ForwardXDifferentiation(buf1);
		vn.CentralYDifferentiation(buf2);
		Tn.ForwardXDifferentiation(buf3);

		for (size_t ii = 0; ii < NN; ii++)
		{
			tauxxn(ii) = lambdan(ii)*(buf1(ii) + buf2(ii)) + static_cast<precisao>(2) * mun(ii)*buf1(ii);
			qxn(ii) = -kn(ii) * buf3(ii);
		}

		un.CentralYDifferentiation(buf1);
		vn.ForwardXDifferentiation(buf2);

		for (size_t ii = 0; ii < NN; ii++)
		{
			tauxyn(ii) = mun(ii)*(buf1(ii) + buf2(ii));
		}

		for (size_t ii = 0; ii < NN; ii++)
		{
			E(ii).cValues[0] = rhon(ii)*un(ii);
			E(ii).cValues[1] = rhon(ii)*un(ii)*un(ii) + Pn(ii) - tauxxn(ii);
			E(ii).cValues[2] = rhon(ii)*un(ii)*vn(ii) - tauxyn(ii);
			E(ii).cValues[3] = (Etn(ii) + Pn(ii)) * un(ii) - un(ii)*tauxxn(ii) - vn(ii)*tauxyn(ii) + qxn(ii);
		}

		un.CentralXDifferentiation(buf1);
		vn.ForwardYDifferentiation(buf2);
		Tn.ForwardYDifferentiation(buf3);

		for (size_t ii = 0; ii < NN; ii++)
		{
			tauyyn(ii) = lambdan(ii)*(buf1(ii) + buf2(ii)) + static_cast<precisao>(2) * mun(ii)*buf2(ii);
			qyn(ii) = -kn(ii) * buf3(ii);
		}

		un.ForwardYDifferentiation(buf1);
		vn.CentralXDifferentiation(buf2);

		for (size_t ii = 0; ii < NN; ii++)
		{
			tauxyn(ii) = mun(ii)*(buf1(ii) + buf2(ii));
		}

		for (size_t ii = 0; ii < NN; ii++)
		{
			F(ii).cValues[0] = rhon(ii)*vn(ii);
			F(ii).cValues[1] = rhon(ii)*un(ii)*vn(ii) - tauxyn(ii);
			F(ii).cValues[2] = rhon(ii)*vn(ii)*vn(ii) + Pn(ii) - tauyyn(ii);
			F(ii).cValues[3] = (Etn(ii) + Pn(ii)) * vn(ii) - un(ii)*tauxyn(ii) - vn(ii)*tauyyn(ii) + qyn(ii);
		}

		// Calcular dU/dt² novamento utilizando os valores predizidos de U

		E.BackwardXDifferentiation(bufV1);
		F.BackwardYDifferentiation(bufV2);

		dUdt2 = -bufV1 - bufV2;		

		// Calcular o valor de U final como o avanço "médio" no tempo
		// Unp = Un + 1/2 * (dU/dt¹ + dU/dt²)

		Unp = Un + (dUdt1 + dUdt2) * dtmin*static_cast<precisao>(0.5);

		for (size_t ii = 0; ii < NN; ii++)
		{
			rhon(ii) = Unp(ii).cValues[0];
			un(ii) = Unp(ii).cValues[1] / Unp(ii).cValues[0];
			vn(ii) = Unp(ii).cValues[2] / Unp(ii).cValues[0];
			Etn(ii) = Unp(ii).cValues[3];

			Tn(ii) = (Etn(ii) / rhon(ii) - (un(ii)*un(ii) + vn(ii)*vn(ii)) / static_cast<precisao>(2.0))/Cv;
			Pn(ii) = rhon(ii)*R * Tn(ii);
			mun(ii) = visc(Tn(ii));
			kn(ii) = condTermica(mun(ii));
			lambdan(ii) = segVisc(mun(ii));

			an(ii) = velSom(Tn(ii));
		}

		Un = Unp;

		garanteCondsContorno(un, vn, Pn, Tn);
		t += dtmin;

		// salva to .csv

		saveToCSV(un, "data/un" + std::to_string(t*10e5) + ".csv");
		saveToCSV(vn, "data/vn" + std::to_string(t*10e5) + ".csv");
		saveToCSV(Tn, "data/Tn" + std::to_string(t*10e5) + ".csv");
		saveToCSV(Pn, "data/Pn" + std::to_string(t*10e5) + ".csv");

		std::cout << "Tempo: " << t << std::endl;
	}

	return 0;
}

void ExibirSpaceFieldScalar(const SpaceFieldScalar& spf)
{
	const size_t nn = spf.GetNN();
	for (size_t ii = 0; ii < nn; ii++)
	{
		//std::cout << spf(ii) << " ; ";
	}
}
