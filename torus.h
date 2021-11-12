// links importantes:
// http://blog.marcinchwedczuk.pl/ray-tracing-torus
// https://en.wikipedia.org/wiki/Quartic_equation#Solving_a_quartic_equation

// Para calcular as raizes de um polinomio quadratico utilizamos a funcao do github:
// https://github.com/sasamil/Quartic


#ifndef TORUS_H
#define TORUS_H

#include "hittable.h"
#include "vec3.h"

#include <math.h>    
#include <complex>

const double PI = 3.141592653589793238463L;
const double M_2PI = 2*PI;
const double eps=1e-12;

typedef std::complex<double> DComplex;

//---------------------------------------------------------------------------
// solve cubic equation x^3 + a*x^2 + b*x + c
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
unsigned int solveP3(double *x,double a,double b,double c) {
	double a2 = a*a;
    	double q  = (a2 - 3*b)/9;
	double r  = (a*(2*a2-9*b) + 27*c)/54;
    	double r2 = r*r;
	double q3 = q*q*q;
	double A,B;
    	if(r2<q3)
    	{
    		double t=r/sqrt(q3);
    		if( t<-1) t=-1;
    		if( t> 1) t= 1;
    		t=acos(t);
    		a/=3; q=-2*sqrt(q);
    		x[0]=q*cos(t/3)-a;
    		x[1]=q*cos((t+M_2PI)/3)-a;
    		x[2]=q*cos((t-M_2PI)/3)-a;
    		return 3;
    	}
    	else
    	{
    		A =-pow(fabs(r)+sqrt(r2-q3),1./3);
    		if( r<0 ) A=-A;
    		B = (0==A ? 0 : q/A);

		a/=3;
		x[0] =(A+B)-a;
		x[1] =-0.5*(A+B)-a;
		x[2] = 0.5*sqrt(3.)*(A-B);
		if(fabs(x[2])<eps) { x[2]=x[1]; return 2; }

		return 1;
        }
}



//---------------------------------------------------------------------------
// Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
// (attention - this function returns dynamically allocated array. It has to be released afterwards)
DComplex* solve_quartic(double a, double b, double c, double d)
{
	double a3 = -b;
	double b3 =  a*c -4.*d;
	double c3 = -a*a*d - c*c + 4.*b*d;

	// cubic resolvent
	// y^3 − b*y^2 + (ac−4d)*y − a^2*d−c^2+4*b*d = 0

	double x3[3];
	unsigned int iZeroes = solveP3(x3, a3, b3, c3);

	double q1, q2, p1, p2, D, sqD, y;

	y = x3[0];
	// THE ESSENCE - choosing Y with maximal absolute value !
	if(iZeroes != 1)
	{
		if(fabs(x3[1]) > fabs(y)) y = x3[1];
		if(fabs(x3[2]) > fabs(y)) y = x3[2];
	}

	// h1+h2 = y && h1*h2 = d  <=>  h^2 -y*h + d = 0    (h === q)

	D = y*y - 4*d;
	if(fabs(D) < eps) //in other words - D==0
	{
		q1 = q2 = y * 0.5;
		// g1+g2 = a && g1+g2 = b-y   <=>   g^2 - a*g + b-y = 0    (p === g)
		D = a*a - 4*(b-y);
		if(fabs(D) < eps) //in other words - D==0
			p1 = p2 = a * 0.5;

		else
		{
			sqD = sqrt(D);
			p1 = (a + sqD) * 0.5;
			p2 = (a - sqD) * 0.5;
		}
	}
	else
	{
		sqD = sqrt(D);
		q1 = (y + sqD) * 0.5;
		q2 = (y - sqD) * 0.5;
		// g1+g2 = a && g1*h2 + g2*h1 = c       ( && g === p )  Krammer
		p1 = (a*q1-c)/(q1-q2);
		p2 = (c-a*q2)/(q1-q2);
	}

    DComplex* retval = new DComplex[4];

	// solving quadratic eq. - x^2 + p1*x + q1 = 0
	D = p1*p1 - 4*q1;
	if(D < 0.0)
	{
		retval[0].real( -p1 * 0.5 );
		retval[0].imag( sqrt(-D) * 0.5 );
		retval[1] = std::conj(retval[0]);
	}
	else
	{
		sqD = sqrt(D);
		retval[0].real( (-p1 + sqD) * 0.5 );
		retval[1].real( (-p1 - sqD) * 0.5 );
	}

	// solving quadratic eq. - x^2 + p2*x + q2 = 0
	D = p2*p2 - 4*q2;
	if(D < 0.0)
	{
		retval[2].real( -p2 * 0.5 );
		retval[2].imag( sqrt(-D) * 0.5 );
		retval[3] = std::conj(retval[2]);
	}
	else
	{
		sqD = sqrt(D);
		retval[2].real( (-p2 + sqD) * 0.5 );
		retval[3].real( (-p2 - sqD) * 0.5 );
	}

    return retval;
}


class torus : public hittable {
   public:
       torus() {}
       torus(point3 cen, double R, double r, shared_ptr<material> m)
           : center(cen), outerRadius(R), innerRadius(r), mat_ptr(m) {};

       virtual bool hit(
           const ray& r, double t_min, double t_max, hit_record& rec) const override;

   public:
       point3 center;
       double outerRadius, innerRadius;
       shared_ptr<material> mat_ptr;
};

bool torus::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {

    double c4 = pow(r.direction()[0] * r.direction()[0] + r.direction()[1] * r.direction()[1] + r.direction()[2] * r.direction()[2], 2);

    double c3 = 4 * (r.origin()[0] * r.direction()[0] + r.origin()[1] * r.direction()[1] + r.origin()[2] * r.direction()[2]) * 
                (r.origin()[0] * r.direction()[0] + r.origin()[1] * r.direction()[1] + r.origin()[2] * r.direction()[2]);

    double c2 = 2 * (r.direction()[0] * r.direction()[0] + r.direction()[1] * r.direction()[1] + r.direction()[2] * r.direction()[2]) * 
                (r.origin()[0] * r.origin()[0] + r.origin()[1] * r.origin()[1] + r.origin()[2] * r.origin()[2] - (pow(innerRadius, 2) + pow(outerRadius, 2))) 
                + 4 * pow(r.origin()[0] * r.direction()[0] + r.origin()[1] * r.direction()[1] + r.origin()[2] * r.direction()[2],2) 
                + 4 * pow(outerRadius, 2) * pow(r.direction()[1], 2);

    double c1 = 4 * (r.origin()[0] * r.origin()[0] + r.origin()[1] * r.origin()[1] + r.origin()[2] * r.origin()[2] - 
                ((innerRadius * innerRadius) + (outerRadius * outerRadius))) * 
                (r.origin()[0] * r.direction()[0] + r.origin()[1] * r.direction()[1] + r.origin()[2] * r.direction()[2]) + 
                8 * pow(outerRadius, 2) * r.origin()[1] * r.direction()[1];

    double c0 = pow((r.origin()[0] * r.origin()[0]  + r.origin()[1] * r.origin()[1]  + r.origin()[2] * r.origin()[2]  - ( pow(outerRadius, 2) + pow(innerRadius, 2))),2) - 
                4 * pow(outerRadius, 2) * ( (pow(innerRadius, 2) - r.origin()[1]*r.origin()[1]));


    double a = c3/ c4;
    double b = c2/ c4;
    double c = c1/ c4;
    double d = c0/ c4;

    // solve the algebric equation of 4th order and print the results
    std::complex<double>* solutions = solve_quartic(a, b, c, d);
    for (int i=0;i<4;i++) {
        if  (solutions[i].imag() == 0) {

            if (solutions[i].real() > t_min && solutions[i].real() < t_max) {
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
