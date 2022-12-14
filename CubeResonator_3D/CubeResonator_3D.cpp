#include "Resonator.h"
#include "Explicit.h"

int main()
{
//	Resonator	cube(1, 0.4, 0.4, 288.15, 0.00005);
	Explicit	scheme(0.4, 0.4, 1, 288.15, 0.00005, 10, 10, 30);
//	Explicit	scheme;

//	std::cout << "c0 = " << cube.get_c0() << std::endl;
//	std::cout << "p0 = " << cube.get_p0() << std::endl;
//	std::cout << "p0 / (ro0 * U0 * c0) = "
//		<< cube.get_p0() / (cube.get_ro0() * cube.get_c0() * cube.get_l() * cube.get_omega()) << std::endl;
//	std::cout << "U0 = " << cube.get_l() * cube.get_omega() << std::endl;
//	std::cout << "c0 / U0 = " << cube.get_c0() / (cube.get_l() * cube.get_omega()) << std::endl;
//	std::cout << "Pr = " << scheme.get_Pr() << std::endl;
/*	std::cout << "X, Y, Z = " << scheme.get_X() << ' ' << scheme.get_Y() << ' ' << scheme.get_Z() << std::endl;
	std::cout << "T, l = " << scheme.get_T() << ' ' << scheme.get_l() << std::endl;
	std::cout << "dx, dy, dz = " << scheme.get_dx() << ' ' << scheme.get_dy() << ' ' << scheme.get_dz() << std::endl;*/
	scheme.calculate_dt();
	std::cout << "dt = " << scheme.get_dt() << std::endl;
	scheme.write_init_file();
	for (unsigned long it = 0; it < 100; it++)
	{
		std::cout << "it = " << it << std::endl;
		scheme.calculate_dU();
		scheme.calculate_U(it);
		scheme.write_in_file("ro");
		scheme.write_in_file("u");
		scheme.write_in_file("v");
		scheme.write_in_file("w");
		scheme.write_in_file("e");
	}
	return 0;
}
