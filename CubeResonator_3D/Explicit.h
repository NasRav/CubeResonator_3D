#pragma once

#include "Resonator.h"
#include <vector>

class Explicit : public Resonator
{
public:
	Explicit(double X, double Y, double Z, double T, double l, int Nx, int Ny, int Nz);

	void	calculate_dt();

	double	L_xx(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_yy(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_zz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_xy(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_xz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_yz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	d_dx(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	d_dy(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	d_dz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	H(int i, int j, int k);

	void	calculate_dU();
	void	calculate_U(unsigned long it);
	double	get_Pr();
	double	get_dt();

private:
	const int		_Nx, _Ny, _Nz;
	const double	_dx, _dy, _dz, _Pr;
	double			_dt;
//	std::vector<std::vector<std::vector<double>>>	_ro;
	std::vector<std::vector<std::vector<double>>>	_u;
	std::vector<std::vector<std::vector<double>>>	_v;
	std::vector<std::vector<std::vector<double>>>	_w;
	std::vector<std::vector<std::vector<double>>>	_e;
//	std::vector<std::vector<std::vector<double>>>	d_ro;
	std::vector<std::vector<std::vector<double>>>	d_u;
	std::vector<std::vector<std::vector<double>>>	d_v;
	std::vector<std::vector<std::vector<double>>>	d_w;
	std::vector<std::vector<std::vector<double>>>	d_e;
};

