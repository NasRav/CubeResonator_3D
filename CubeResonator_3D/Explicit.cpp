#include "Explicit.h"

Explicit::Explicit(double X, double Y, double Z, double T, double l, int Nx, int Ny, int Nz) :
	Resonator(X, Y, Z, T, l), _Nx(Nx), _Ny(Ny), _Nz(Nz),
	_dx(_X / (_Nx - 1)), _dy(_Y / (_Ny - 1)), _dz(_Z / (_Nz - 1)),
	_Pr(abs(_mu / _ro0 / (1.302e-10 * pow(_T, 2) + 5.893e-6 * _T - 7283e-6)))
{
	_ro.resize(_Nx);
	_u.resize(_Nx);
	_v.resize(_Nx);
	_w.resize(_Nx);
	_e.resize(_Nx);
	d_ro.resize(_Nx);
	d_u.resize(_Nx);
	d_v.resize(_Nx);
	d_w.resize(_Nx);
	d_e.resize(_Nx);
	for (int i = 0; i < _Nx; i++)
	{
		_ro[i].resize(_Ny);
		_u[i].resize(_Ny);
		_v[i].resize(_Ny);
		_w[i].resize(_Ny);
		_e[i].resize(_Ny);
		d_ro[i].resize(_Ny);
		d_u[i].resize(_Ny);
		d_v[i].resize(_Ny);
		d_w[i].resize(_Ny);
		d_e[i].resize(_Ny);
	}
	for (int i = 0; i < _Nx; i++)
		for (int j = 0; j < _Ny; j++)
		{
			_ro[i][j].resize(_Nz);
			_u[i][j].resize(_Nz);
			_v[i][j].resize(_Nz);
			_w[i][j].resize(_Nz);
			_e[i][j].resize(_Nz);
			d_ro[i][j].resize(_Nz);
			d_u[i][j].resize(_Nz);
			d_v[i][j].resize(_Nz);
			d_w[i][j].resize(_Nz);
			d_e[i][j].resize(_Nz);
		}
}

void Explicit::calculate_dt()
{
	double	d_min(_dx);

	if (_dy < d_min)
		d_min = _dy;
	else if (_dz < d_min)
		d_min = _dz;
	if (_Pr / (6 * _gamma) < 0.15)
		this->_dt = _Pr / (6 * _gamma) * _ro0 * pow(d_min, 2) / _mu;
	else
		this->_dt = 0.15 * _ro0 * pow(_dx, 2) / _mu;
}

void Explicit::calculate_dU()
{
	for (int i = 1; i < _Nx - 1; i++)
		for (int j = 1; j < _Ny - 1; j++)
			for (int k = 1; k < _Nz - 1; k++)
			{
				d_u[i][j][k] = _mu * _dt / (_ro0 * _dy * _dz)
					* (4 / 3 * L_xx(_u, i, j, k) + 1 / 3 * L_xy(_v, i, j, k) + 1 / 3 * L_xz(_w, i, j, k)
						+ L_yy(_u, i, j, k) + L_zz(_u, i, j, k));
				d_v[i][j][k] = _mu * _dt / (_ro0 * _dx * _dz)
					* (4 / 3 * L_yy(_v, i, j, k) + 1 / 3 * L_xy(_u, i, j, k) + 1 / 3 * L_yz(_w, i, j, k)
						+ L_xx(_v, i, j, k) + L_zz(_v, i, j, k));
				d_w[i][j][k] = _mu * _dt / (_ro0 * _dx * _dy)
					* (4 / 3 * L_zz(_w, i, j, k) + 1 / 3 * L_xz(_u, i, j, k) + 1 / 3 * L_yz(_v, i, j, k)
						+ L_xx(_w, i, j, k) + L_yy(_w, i, j, k));
//				d_e[i][j][k] = _mu * _dt / (_ro0 * _dx * _dy) * _gamma / _Pr
//					* (L_xx(_e, i, j, k) + L_yy(_e, i, j, k) + L_zz(_e, i, j, k) + H(i, j, k));
			}
}

void Explicit::calculate_U(unsigned long it)
{
	for (int i = 1; i < _Nx - 1; i++)
		for (int j = 1; j < _Ny - 1; j++)
			for (int k = 1; k < _Nz - 1; k++)
			{
				_u[i][j][k] += d_u[i][j][k];
				_v[i][j][k] += d_v[i][j][k];
				_w[i][j][k] += d_w[i][j][k];
//				_e[i][j][k] += d_e[i][j][k];
			}
	for (int i = 1; i < _Nx - 1; i++)
		for (int j = 1; j < _Ny - 1; j++)
			_u[i][j][0] = _l * _omega * sin(_omega * _dt * it);
}

double Explicit::get_Pr()
{
	return this->_Pr;
}

double Explicit::get_dt()
{
	return this->_dt;
}

double Explicit::L_xx(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k)
{
	return f[i + 1][j][k] - 2 * f[i][j][k] + f[i - 1][j][k];
}

double Explicit::L_yy(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k)
{
	return f[i][j + 1][k] - 2 * f[i][j][k] + f[i][j - 1][k];
}

double Explicit::L_zz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k)
{
	return f[i][j][k + 1] - 2 * f[i][j][k] + f[i][j][k - 1];
}

double Explicit::L_xy(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k)
{
	return 0.25 * (f[i + 1][j + 1][k] - f[i + 1][j - 1][k] - f[i - 1][j + 1][k] + f[i - 1][j - 1][k]);
}

double Explicit::L_xz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k)
{
	return 0.25 * (f[i + 1][j][k + 1] - f[i + 1][j][k - 1] - f[i - 1][j][k + 1] + f[i - 1][j][k - 1]);
}

double Explicit::L_yz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k)
{
	return 0.25 * (f[i][j + 1][k + 1] - f[i][j + 1][k - 1] - f[i][j - 1][k + 1] + f[i][j - 1][k - 1]);
}
