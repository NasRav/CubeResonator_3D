#pragma once

#include <iostream>
#include <cmath>

class Resonator
{
public:
	Resonator();
	Resonator(double X, double Y, double Z, double T, double l);

	double	get_X();
	double	get_Y();
	double	get_Z();
	double	get_T();
	double	get_l();
	double	get_ro0();
	double	get_mu();
	double	get_p0();
	double	get_c0();
	double	get_PI();
	double	get_omega();
	double	get_delta();
	void	set_X(double X);
	void	set_Z(double Y);
	void	set_Y(double Z);
	void	set_T(double T);
	void	set_l(double l);
	void	set_ro0(double ro0);
	void	set_mu(double mu);
	void	set_p0(double p0);
	void	set_c0(double c0);
	void	set_PI(double PI);
	void	set_omega(double omega);
	void	set_delta(double delta);

protected:
	double	_X, _Y, _Z, _T, _l, _ro0, _mu, _R_gas,
		_M_mol, _p0, _gamma, _c0, _PI, _omega, _delta;
};
