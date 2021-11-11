#ifndef TORUS_H
#define TORUS_H

#include "hittable.h"
#include "vec3.h"

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

   vec3 oc = r.origin() - center;
   auto a = r.direction().length_squared();
   auto half_b = dot(oc, r.direction());
   auto c = oc.length_squared() - outerRadius*outerRadius;

   auto discriminant = half_b*half_b - a*c;
   if (discriminant < 0) return false;
   auto sqrtd = sqrt(discriminant);

   // Find the nearest root that lies in the acceptable range.
   auto root = (-half_b - sqrtd) / a;
   if (root < t_min || t_max < root) {
       root = (-half_b + sqrtd) / a;
       if (root < t_min || t_max < root)
           return false;
   }

   rec.t = root;
   rec.p = r.at(rec.t);
   vec3 outward_normal = (rec.p - center) / outerRadius;
   rec.set_face_normal(r, outward_normal);
   rec.mat_ptr = mat_ptr;

   return true;
}

#endif
