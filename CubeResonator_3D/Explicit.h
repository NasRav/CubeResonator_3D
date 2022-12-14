#pragma once

#include "Resonator.h"
#include <vector>
#include <omp.h>
#include <fstream>
#include <string>

class Explicit : public Resonator
{
public:
	Explicit();
	Explicit(double X, double Y, double Z, double T, double l, int Nx, int Ny, int Nz);

	void	calculate_dt();
	void	calculate_dU();
	void	calculate_U(unsigned long it);

	double	L_xx(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_yy(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_zz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_xy(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_xz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_yz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_x(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_y(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	L_z(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k);
	double	H(int i, int j, int k);

	double	get_Pr();
	double	get_dt();
	double	get_dx();
	double	get_dy();
	double	get_dz();
	void	write_init_file();
	void	write_in_file(std::string name);
	std::vector<std::vector<std::vector<double>>>	read_from_file(std::string name);

private:
	int		_Nx, _Ny, _Nz;
	double	_dx, _dy, _dz, _Pr, _U0, _dt;
	std::vector<std::vector<std::vector<double>>>	_ro;
	std::vector<std::vector<std::vector<double>>>	_u;
	std::vector<std::vector<std::vector<double>>>	_v;
	std::vector<std::vector<std::vector<double>>>	_w;
	std::vector<std::vector<std::vector<double>>>	_e;
	std::vector<std::vector<std::vector<double>>>	d_ro;
	std::vector<std::vector<std::vector<double>>>	d_u;
	std::vector<std::vector<std::vector<double>>>	d_v;
	std::vector<std::vector<std::vector<double>>>	d_w;
	std::vector<std::vector<std::vector<double>>>	d_e;
};
