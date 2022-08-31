#include "Resonator.h"

Resonator::Resonator() {
	_X = 0.4;
	_Y = 0.4;
	_Z = 1;
	_T = 288.15;
	_l = 0.00005;
	_ro0 = 1.225;
	_mu = 1.82e-5;
	_R_gas = 8.31;
	_M_mol = 0.029;
	_PI = acos(-1);
	_gamma = 1.4;
	_p0 = _ro0 * _R_gas * _T / _M_mol;
	_c0 = sqrt(_gamma * _p0 / _ro0);
	_omega = _PI * _c0 / _X;
	_delta = sqrt(2 * _mu / (_omega * _ro0));
}

Resonator::Resonator(double X, double Y, double Z, double T, double l) :
	_X(X), _Y(Y), _Z(Z), _T(T), _l(l), _ro0(1.225), _mu(1.82e-5),
	_R_gas(8.31), _M_mol(0.029), _PI(acos(-1)), _gamma(1.4),
	_p0(_ro0 * _R_gas * _T / _M_mol), _c0(sqrt(_gamma * _p0 / _ro0)),
	_omega(_PI * _c0 / _X), _delta(sqrt(2 * _mu / (_omega * _ro0))) {}

double	Resonator::get_X() {
	return this->_X;
}

double	Resonator::get_Y() {
	return this->_Y;
}

double	Resonator::get_Z() {
	return this->_Z;
}

double	Resonator::get_T() {
	return this->_T;
}

double	Resonator::get_l() {
	return this->_l;
}

double	Resonator::get_p0() {
	return this->_p0;
}

double	Resonator::get_ro0() {
	return this->_ro0;
}

double	Resonator::get_c0() {
	return this->_c0;
}

double	Resonator::get_mu() {
	return this->_mu;
}

double	Resonator::get_PI() {
	return this->_PI;
}

double	Resonator::get_omega() {
	return this->_omega;
}

double	Resonator::get_delta() {
	return this->_delta;
}

void	Resonator::set_X(double X) {
	this->_X = X;
}

void	Resonator::set_Y(double Y) {
	this->_Y = Y;
}

void	Resonator::set_Z(double Z) {
	this->_Z = Z;
}

void	Resonator::set_T(double T) {
	this->_T = T;
}

void	Resonator::set_l(double l) {
	this->_l = l;
}

void	Resonator::set_p0(double p0) {
	this->_p0 = p0;
}

void	Resonator::set_ro0(double ro0) {
	this->_ro0 = ro0;
}

void	Resonator::set_c0(double c0) {
	this->_c0 = c0;
}

void	Resonator::set_mu(double mu) {
	this->_mu = mu;
}

void	Resonator::set_PI(double PI) {
	this->_PI = PI;
}

void	Resonator::set_omega(double omega) {
	this->_omega = omega;
}

void	Resonator::set_delta(double delta) {
	this->_delta = delta;
}
