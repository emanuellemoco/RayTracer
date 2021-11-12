// links importantes:
// http://blog.marcinchwedczuk.pl/ray-tracing-torus
// https://en.wikipedia.org/wiki/Quartic_equation#Solving_a_quartic_equation

#ifndef TORUS_H
#define TORUS_H

#include "hittable.h"
#include "vec3.h"
#include "quartic.h"

#include <math.h>
#include <complex>

class torus : public hittable
{
public:
	torus() {}
	torus(point3 cen, double R, double r, shared_ptr<material> m)
		: center(cen), outerRadius(R), innerRadius(r), mat_ptr(m){};

	virtual bool hit(
		const ray &r, double t_min, double t_max, hit_record &rec) const override;

public:
	point3 center;
	double outerRadius, innerRadius;
	shared_ptr<material> mat_ptr;
};

bool torus::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{
	double ox = r.origin()[0];
	double oy = r.origin()[1];
	double oz = r.origin()[2];

	double dx = r.direction()[0];
	double dy = r.direction()[1];
	double dz = r.direction()[2];

	double sum_direction_square = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
	double origin_direction = (ox * dx + oy * dy + oz * dz);

	double c4 = pow(sum_direction_square, 2);

	double c3 = 4 * sum_direction_square * origin_direction;

	double c2 = 2 * sum_direction_square *
					(pow(ox, 2) + pow(oy, 2) + pow(oz, 2) - (pow(innerRadius, 2) + pow(outerRadius, 2))) +
				4 * pow(origin_direction, 2) + 4 * pow(outerRadius, 2) * pow(dy, 2);

	double c1 = 4 * (ox * ox + oy * oy + oz * oz - (pow(innerRadius, 2) + pow(outerRadius, 2))) *
					origin_direction +
				8 * pow(outerRadius, 2) * oy * dy;

	double c0 = pow((pow(ox, 2) + pow(oy, 2) + pow(oz, 2) - (pow(outerRadius, 2) + pow(innerRadius, 2))), 2) -
				4 * pow(outerRadius, 2) * ((pow(innerRadius, 2) - pow(oy, 2)));

	double a = c3 / c4;
	double b = c2 / c4;
	double c = c1 / c4;
	double d = c0 / c4;

	// solve the algebric equation of 4th order and print the results
	std::complex<double> *solutions = solve_quartic(a, b, c, d);
	for (int i = 0; i < 4; i++)
	{
		if (solutions[i].imag() == 0)
		{

			if (solutions[i].real() > t_min && solutions[i].real() < t_max)
			{
				// std::cerr  << "Entrei no iff" << std::flush;
				rec.t = solutions[i].real();
				rec.p = r.at(rec.t);
				rec.set_face_normal(r, (rec.p - center) / outerRadius);
				rec.mat_ptr = mat_ptr;
				delete[] solutions;

				return true;
			}
		}
	}
	delete[] solutions;

	return false;
}
#endif
