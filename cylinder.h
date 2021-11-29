#ifndef CYLINDER_H
#define CYLINDER_H


#include "rtweekend.h"

#include "hittable.h"


class cylinder : public hittable {
    public:
        cylinder() {}

        cylinder(point3 cen, double r, double h, shared_ptr<material> m)
            : center(cen), radius(r), height(h), mat_ptr(m) {};

        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const override;

    public:
        point3 center;
        double radius;
        double height;
        shared_ptr<material> mat_ptr;

};


bool cylinder::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {

    auto o = r.origin();
    auto d = r.direction();
    auto f = center;

    auto a = d[0]*d[0] + d[2]*d[2];
    
    auto b = 2*o[0]*d[0] - 2*f[0]*d[0] + 2*o[2]*d[2] - 2*f[2]*d[2];

    auto c = o[0]*o[0] - 2*o[0]*f[0] - f[0]*f[0] + o[2]*o[2] - 2*o[2]*f[2] - f[2]*f[2] - radius*radius;

    auto discriminant = b*b - 4*a*c;
    if (discriminant < 0) return false;
    auto sqrtd = sqrt(discriminant);

    auto root = (-b - sqrtd) / (2*a);
    auto p = r.at(root);

    if (root < t_min || t_max < root || p[1] > height/2 || p[1] < -height/2) {

        root = (-b + sqrtd) / (2*a);
        p = r.at(root);

        if (root < t_min || t_max < root || p[1] > height/2 || p[1] < -height/2)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);

    rec.mat_ptr = mat_ptr;

    return true;
}


#endif
