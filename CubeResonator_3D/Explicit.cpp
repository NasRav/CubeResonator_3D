#include "Explicit.h"

Explicit::Explicit() :
	Resonator()
{
	std::cout << "here" << std::endl;
	_Pr = 0.7;
/*	std::ifstream	in("init.txt");
	std::string		line;

	if (in.is_open())
		getline(in, line);
	in.close();
	std::cout << line << std::endl;*/
}

Explicit::Explicit(double X, double Y, double Z, double T, double l, int Nx, int Ny, int Nz) :
	Resonator(X, Y, Z, T, l), _Nx(Nx), _Ny(Ny), _Nz(Nz),
	_dx(_X / (_Nx - 1)), _dy(_Y / (_Ny - 1)), _dz(_Z / (_Nz - 1)),
	_Pr(0.704), _U0(_l * _omega)
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
#pragma omp parallel for
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
#pragma omp parallel for
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
	for (int i = 0; i < _Nx; i++)
		for (int j = 0; j < _Ny; j++)
#pragma omp parallel for
			for (int k = 0; k < _Nz; k++)
			{
				_ro[i][j][k] = _ro0;
				_e[i][j][k] = 2.5 * _R_gas * _T;
			}
}

void Explicit::calculate_dt() {
	double	d_min(_dx);

	if (_dy < d_min)
		d_min = _dy;
	else if (_dz < d_min)
		d_min = _dz;
	if (_Pr / (6 * _gamma) < 0.15)
		this->_dt = _Pr / (6 * _gamma) * _ro0 * pow(d_min, 2) / _mu;
	else
		this->_dt = 0.15 * _ro0 * pow(d_min, 2) / _mu;
}

void Explicit::calculate_dU() {
	for (int i = 1; i < _Nx - 1; i++)
		for (int j = 1; j < _Ny - 1; j++)
			for (int k = 1; k < _Nz - 1; k++)
			{
#pragma omp parallel sections
				{
#pragma omp section
					{
						d_ro[i][j][k] = -_ro[i][j][k]
							* (_u[i][j][k] * L_x(_ro, i, j, k) + _v[i][j][k] * L_y(_ro, i, j, k) + _w[i][j][k] * L_z(_ro, i, j, k))
							- (L_x(_u, i, j, k) + L_y(_v, i, j, k) + L_z(_w, i, j, k));
					}
#pragma omp section
					{
						d_u[i][j][k] = _mu * _dt / (_ro[i][j][k] * _dy * _dz)
							* (4 / 3 * L_xx(_u, i, j, k) + 1 / 3 * L_xy(_v, i, j, k) + 1 / 3 * L_xz(_w, i, j, k)
								+ L_yy(_u, i, j, k) + L_zz(_u, i, j, k));
					}
#pragma omp section
					{
						d_v[i][j][k] = _mu * _dt / (_ro[i][j][k] * _dx * _dz)
							* (4 / 3 * L_yy(_v, i, j, k) + 1 / 3 * L_xy(_u, i, j, k) + 1 / 3 * L_yz(_w, i, j, k)
								+ L_xx(_v, i, j, k) + L_zz(_v, i, j, k));
					}
#pragma omp section
					{
						d_w[i][j][k] = _mu * _dt / (_ro[i][j][k] * _dx * _dy)
							* (4 / 3 * L_zz(_w, i, j, k) + 1 / 3 * L_xz(_u, i, j, k) + 1 / 3 * L_yz(_v, i, j, k)
								+ L_xx(_w, i, j, k) + L_yy(_w, i, j, k));
					}
#pragma omp section
					{
						d_e[i][j][k] = _mu * _dt / (_ro[i][j][k] * _dx * _dy) * _gamma / _Pr
							* (L_xx(_e, i, j, k) + L_yy(_e, i, j, k) + L_zz(_e, i, j, k)) + H(i, j, k);
					}
				}
			}
}

void Explicit::calculate_U(unsigned long it) {
	double	t(_dt * it);

	for (int i = 1; i < _Nx - 1; i++)
		for (int j = 1; j < _Ny - 1; j++)
#pragma omp parallel for
			for (int k = 1; k < _Nz - 1; k++)
			{
				_ro[i][j][k] += d_ro[i][j][k];
				_u[i][j][k] += d_u[i][j][k];
				_v[i][j][k] += d_v[i][j][k];
				_w[i][j][k] += d_w[i][j][k];
				_e[i][j][k] += d_e[i][j][k];
			}
	for (int i = 1; i < _Nx - 1; i++)
#pragma omp parallel for
		for (int j = 1; j < _Ny - 1; j++)
		{
			_u[i][j][0] = _U0 * sin(_omega * t);
			_ro[i][j][0] = _ro[i][j][1];
			_e[i][j][0] = _e[i][j][1];
			_ro[i][j][_Nz - 1] = _ro[i][j][_Nz - 2];
			_e[i][j][_Nz - 1] = _e[i][j][_Nz - 2];
		}
	for (int i = 0; i < _Nx; i++)
#pragma omp parallel for
		for (int k = 0; k < _Nz; k++)
		{
			_ro[i][0][k] = _ro[i][1][k];
			_e[i][0][k] = _e[i][1][k];
			_ro[i][_Ny - 1][k] = _ro[i][_Ny - 2][k];
			_e[i][_Ny - 1][k] = _e[i][_Ny - 2][k];
		}
	for (int j = 0; j < _Ny; j++)
#pragma omp parallel for
		for (int k = 0; k < _Nz; k++)
		{
			_ro[0][j][k] = _ro[1][j][k];
			_e[0][j][k] = _e[1][j][k];
			_ro[_Nx - 1][j][k] = _ro[_Nx - 2][j][k];
			_e[_Nx - 1][j][k] = _e[_Nx - 2][j][k];
		}
}

double Explicit::L_x(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k) {
	return (f[i][j][k] - f[i - 1][j][k]) / _dx;
}

double Explicit::L_y(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k) {
	return (f[i][j][k] - f[i][j - 1][k]) / _dy;
}

double Explicit::L_z(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k) {
	return (f[i][j][k] - f[i][j][k - 1]) / _dz;
}

double Explicit::H(int i, int j, int k) {
	return 2 * ( pow(L_x(_u, i, j, k), 2) + pow(L_y(_v, i, j, k), 2) + pow(L_z(_w, i, j, k), 2)
		+ 0.5 * pow(L_y(_u, i, j, k) + L_x(_v, i, j, k), 2)
		+ 0.5 * pow(L_z(_v, i, j, k) + L_y(_w, i, j, k), 2)
		+ 0.5 * pow(L_x(_w, i, j, k) + L_z(_u, i, j, k), 2) )
		- 2 / 3 * pow(L_x(_u, i, j, k) + L_y(_v, i, j, k) + L_z(_w, i, j, k), 2);
}

double Explicit::L_xx(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k) {
	return f[i + 1][j][k] - 2 * f[i][j][k] + f[i - 1][j][k];
}

double Explicit::L_yy(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k) {
	return f[i][j + 1][k] - 2 * f[i][j][k] + f[i][j - 1][k];
}

double Explicit::L_zz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k) {
	return f[i][j][k + 1] - 2 * f[i][j][k] + f[i][j][k - 1];
}

double Explicit::L_xy(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k) {
	return 0.25 * (f[i + 1][j + 1][k] - f[i + 1][j - 1][k] - f[i - 1][j + 1][k] + f[i - 1][j - 1][k]);
}

double Explicit::L_xz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k) {
	return 0.25 * (f[i + 1][j][k + 1] - f[i + 1][j][k - 1] - f[i - 1][j][k + 1] + f[i - 1][j][k - 1]);
}

double Explicit::L_yz(std::vector<std::vector<std::vector<double>>> f, int i, int j, int k) {
	return 0.25 * (f[i][j + 1][k + 1] - f[i][j + 1][k - 1] - f[i][j - 1][k + 1] + f[i][j - 1][k - 1]);
}

double Explicit::get_Pr() {
	return this->_Pr;
}

double Explicit::get_dt() {
	return this->_dt;
}

void Explicit::write_init_file() {
	std::ofstream	f_out("init.txt");

	f_out << _X << ' ' << _Y << ' ' << _Z << ' ' << _T << ' ' << _l << ' ' << _Nx << ' ' << _Ny << ' ' << _Nz << std::endl;
}

void Explicit::write_in_file(std::string name) {
	std::vector<std::vector<std::vector<double>>>	array;
	if (name == "ro")
		array = _ro;
	else if (name == "u")
		array = _u;
	else if (name == "v")
		array = _v;
	else if (name == "w")
		array = _w;
	else if (name == "e")
		array = _e;
	else
		return;

	std::ofstream	f_out(name + "_3D_NonLinear_X=" + std::to_string(_X) + "_Y=" + std::to_string(_Y) + "_Z=" + std::to_string(_Z) + ".txt");
	for (int k = 0; k < _Nz; k++)
	{
		for (int j = 0; j < _Ny; j++)
		{
			for (int i = 0; i < _Nx; i++)
				f_out << array[i][j][k] << ' ';
			f_out << std::endl;
		}
		f_out << std::endl;
	}
	f_out.close();
}
