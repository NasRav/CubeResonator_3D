#include "Resonator.h"
#include "Explicit.h"

int main()
{
	Resonator   cube(1, 0.4, 0.4, 288.15, 0.0003);
	Explicit	scheme(0.4, 0.4, 1, 288.15, 0.0003, 40, 40, 100);

	std::cout << "c0 = " << cube.get_c0() << std::endl;
	std::cout << "p0 = " << cube.get_p0() << std::endl;
	std::cout << "p0 / (ro0 * U0 * c0) = "
		<< cube.get_p0() / (cube.get_ro0() * cube.get_c0() * cube.get_l() * cube.get_omega()) << std::endl;
	std::cout << "U0 = " << cube.get_l() * cube.get_omega() << std::endl;
	std::cout << "c0 / U0 = " << cube.get_c0() / (cube.get_l() * cube.get_omega()) << std::endl;
	std::cout << "Pr = " << scheme.get_Pr() << std::endl;
	scheme.calculate_dt();
	std::cout << "dt = " << scheme.get_dt() << std::endl;
	for (unsigned long it = 0; it < 100; it++)
	{
		std::cout << "it = " << it << std::endl;
		scheme.calculate_dU();
		scheme.calculate_U(it);
	}
	return 0;
}
