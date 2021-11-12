#include "rtweekend.h"

#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "torus.h"
#include "ray.h"
#include "vec3.h"
#include "camera.h"
#include "material.h"

#include <iostream>
#include <fstream> // para ler e gravar em arquivos.
#include <iomanip>

// Função auxiliar para criar um fundo de imagem colorido
// pega a direcao do raio e calcula uma interpolacao entra banco e azul
color ray_color(const ray &r, const hittable &world, int depth)
{
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

    if (world.hit(r, 0.001, infinity, rec))
    {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth - 1);
        return color(0, 0, 0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}

hittable_list random_scene()
{
    hittable_list world;

    auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    auto material_center = make_shared<lambertian>(color(0.7, 0.1, 0.7));
    auto material_left = make_shared<metal>(color(0.0, 0.8, 0.8), 0.3);
    auto material_right = make_shared<metal>(color(0.8, 0.6, 0.2), 1.0);

    point3 center(2 + 0.9 * random_double(), 1, 1 + 0.9 * random_double());
    shared_ptr<material> sphere_material;
    auto albedo = color::random() * color::random();
    sphere_material = make_shared<lambertian>(albedo);
    world.add(make_shared<torus>(center, 1.2, 0.2, material_center));



    auto material_center_pink = make_shared<lambertian>(color(1.0, 0.8, 0.9));

    world.add(make_shared<torus>(center, 0.4, 0.1, material_center_pink));




    return world;
}

int main()
{
    // Image
    const auto aspect_ratio = 3.0 / 2.0;
    const int image_width = 1200;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 10;
    const int max_depth = 50;

    // World
    auto world = random_scene();

    // Camera
    // camera cam(point3(-2,2,1), point3(0,0,-1), vec3(0,1,0), 90, aspect_ratio);
    // camera cam(point3(-2,2,1), point3(0,0,-1), vec3(0,1,0), 20, aspect_ratio); // mudando o fiel of view
    point3 lookfrom(8, -5, 2);
    point3 lookat(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;
    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);

    // Render
    std::ofstream file;     // Cria um stream para arquivos
    file.open("image.ppm"); // Abre um arquivo para saída do stream

    file << "P3\n"
         << image_width << ' ' << image_height << "\n255\n";
    for (int j = image_height - 1; j >= 0; --j)
    {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i)
        {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s)
            {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(file, pixel_color, samples_per_pixel);
        }
    }

    file.close(); // Fecha o stream para arquivo
    std::cerr << "\nDone.\n";
}
