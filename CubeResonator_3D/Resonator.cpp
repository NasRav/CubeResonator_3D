#include "Resonator.h"

Resonator::Resonator(double X, double Y, double Z, double T, double l) :
	_X(X), _Y(Y), _Z(Z), _T(T), _l(l), _ro0(1.225), _mu(1.82e-5),
	_R_gas(8.31), _M_mol(0.029), _PI(acos(-1)), _gamma(1.4),
	_p0(_ro0 * _R_gas * _T / _M_mol), _c0(sqrt(_gamma * _p0 / _ro0)),
	_omega(_PI * _c0 / _X), _delta(sqrt(2 * _mu / (_omega * _ro0)))
{}

double	Resonator::get_X()
{
	return this->_X;
}

double	Resonator::get_Y()
{
	return this->_Y;
}

double	Resonator::get_Z()
{
	return this->_Z;
}

double	Resonator::get_T()
{
	return this->_T;
}

double	Resonator::get_l()
{
	return this->_l;
}

double	Resonator::get_p0()
{
	return this->_p0;
}

double	Resonator::get_ro0()
{
	return this->_ro0;
}

double	Resonator::get_c0()
{
	return this->_c0;
}

double	Resonator::get_mu()
{
	return this->_mu;
}

double	Resonator::get_PI()
{
	return this->_PI;
}

double	Resonator::get_omega()
{
	return this->_omega;
}

double	Resonator::get_delta()
{
	return this->_delta;
}
