#include "material.h"
#include "muni/camera.h"
#include "muni/common.h"
#include "muni/image.h"
#include "muni/material.h"
#include "muni/math_helpers.h"
#include "muni/obj_loader.h"
#include "muni/ray_tracer.h"
#include "muni/sampler.h"
#include "muni/scenes/box.h"
#include "muni/triangle.h"
#include "ray_tracer.h"
#include "spdlog/spdlog.h"
#include "triangle.h"
#include <cmath>
#include <iostream>

using namespace muni;

RayTracer::Octree octree{};

/** Offset the ray origin to avoid self-intersection.
    \param[in] ray_pos The original ray origin.
    \param[in] normal The normal of the surface at the hit point.
    \return The offset ray origin.
*/
Vec3f offset_ray_origin(Vec3f ray_pos, Vec3f normal) {
    return ray_pos + EPS * normal;
}

/** Check if the triangle is an emitter.
    \param[in] tri The triangle to check
    \return True if the triangle is an emitter, false otherwise.
*/
bool is_emitter(const Triangle &tri) { return tri.emission != Vec3f{0.0f}; }

/** Evaluate the radiance of the area light. We **do not** check whether the hit
 point is on the light source, so make sure
 *  the hit point is on the light source before calling this function.
    \param[in] light_dir The **outgoing** direction from the light source to the
 scene. \return The radiance of the light source.
*/
Vec3f eval_area_light(const Vec3f light_dir) {
    if (dot(light_dir, BoxScene::light_normal) > 0.0f)
        return BoxScene::light_color;
    return Vec3f{0.0f};
}

/** Sample a point on the area light with a uniform distribution.
    \param[in] samples A 2D uniform random sample.
    \return A tuple containing the sampled position, the normal of the light
 source, and the PDF value.
*/
std::tuple<Vec3f, Vec3f, float> sample_area_light(Vec2f samples) {
    float x = BoxScene::light_x + BoxScene::light_len_x * samples.x;
    float y = BoxScene::light_y + BoxScene::light_len_x * samples.y;
    float z = BoxScene::light_z;
    return std::make_tuple(Vec3f{x,y,z}, BoxScene::light_normal, BoxScene::inv_light_area);
}


Vec3f shade(Triangle tri, Vec3f p, Vec3f wo) {
    // Contribution from the light source (light sampling)
    Vec3f L_dir{0.0f};

    // Uniformly sample the light at x
    auto [x, n_x, p_x] = sample_area_light(UniformSampler::next2d());
    Vec3f wx = normalize(p - x);
    float p_wx_light = BoxScene::inv_light_area * std::pow(distance(p,x), 2.0) / dot(wx, n_x);
    float p_wx_BRDF = (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]) ? 
            std::get<Lambertian>(BoxScene::materials[tri.material_id]).pdf(wo, wx, tri.face_normal) :
            (std::holds_alternative<Microfacet>(BoxScene::materials[tri.material_id]) ? 
            std::get<Microfacet>(BoxScene::materials[tri.material_id]).pdf(wo, wx, tri.face_normal) :
            std::get<Dielectric>(BoxScene::materials[tri.material_id]).pdf(wo, wx, tri.face_normal)));
    float w_light = p_wx_light / (p_wx_light + p_wx_BRDF);

    // Shoot a ray from p to x
    const auto [_, __, emitter_tri] = 
        RayTracer::closest_hit(p, -wx, octree, BoxScene::triangles);

    // If the ray is not blocked in the middle
    if (is_emitter(emitter_tri) && dot(-wx, tri.face_normal) > 0) {
        if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]))
            L_dir = eval_area_light(wx) 
                * std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval()
                * std::max(0.0f, dot(-wx, tri.face_normal)) / p_wx_light * w_light;
        else if (std::holds_alternative<Microfacet>(BoxScene::materials[tri.material_id]))
            L_dir = eval_area_light(wx) 
                * std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, -wx, tri.face_normal)
                * std::max(0.0f,dot(-wx, tri.face_normal)) / p_wx_light * w_light;
        else if (std::holds_alternative<Dielectric>(BoxScene::materials[tri.material_id]))
            L_dir = eval_area_light(wx) 
                * std::get<Dielectric>(BoxScene::materials[tri.material_id]).eval(wo, -wx, tri.face_normal)
                * std::max(0.0f,dot(-wx, tri.face_normal)) / p_wx_light * w_light;
    }
    else if (std::holds_alternative<Dielectric>(BoxScene::materials[emitter_tri.material_id])) {
        if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]))
            L_dir = eval_area_light(wx) 
                * std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval()
                * std::max(0.0f, dot(-wx, tri.face_normal)) / p_wx_light * w_light * std::get<Dielectric>(BoxScene::materials[emitter_tri.material_id]).attenuation;
        else if (std::holds_alternative<Microfacet>(BoxScene::materials[tri.material_id]))
            L_dir = eval_area_light(wx) 
                * std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, -wx, tri.face_normal)
                * std::max(0.0f,dot(-wx, tri.face_normal)) / p_wx_light * w_light * std::get<Dielectric>(BoxScene::materials[emitter_tri.material_id]).attenuation;
        else if (std::holds_alternative<Dielectric>(BoxScene::materials[tri.material_id]))
            L_dir = eval_area_light(wx) 
                * std::get<Dielectric>(BoxScene::materials[tri.material_id]).eval(wo, -wx, tri.face_normal)
                * std::max(0.0f,dot(-wx, tri.face_normal)) / p_wx_light * w_light * std::get<Dielectric>(BoxScene::materials[emitter_tri.material_id]).attenuation;
    }

    // Contribution from other reflectors
    Vec3f L_indir{0.0f};
    
    // Test Russian Roulette with probability p_rr = 0.8f
    const float p_rr = 0.8f;
    if (UniformSampler::next1d() > p_rr) {}

    // Lambertian material
    else if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id])) {
        // Randomly choose one direction wi with BRDF sampling
        auto [wi, p_wi_BRDF] = 
            std::get<Lambertian>(BoxScene::materials[tri.material_id])
            .sample(tri.face_normal, UniformSampler::next2d());

        // Trace the new ray
        const auto [is_ray_hit, t_min, nearest_tri] = 
            RayTracer::closest_hit(p, wi, octree, BoxScene::triangles);
        const Vec3f q = p + t_min * wi;

        float p_wi_light = BoxScene::inv_light_area * std::pow(distance(p,q), 2.0) / dot(wi, n_x);
        float w_BRDF = p_wi_BRDF / (p_wi_light + p_wi_BRDF);
        
        
        // If the ray hit a non-emitting object at q
        if (is_ray_hit && !is_emitter(nearest_tri)) {
            L_indir = shade(nearest_tri, offset_ray_origin(q, nearest_tri.face_normal), -wi)
                * std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval()
                * std::max(0.0f, dot(wi, tri.face_normal)) / (p_wi_BRDF * p_rr);
        }
        
        // If the ray hit an emitting object at q
        else if (is_ray_hit && is_emitter(nearest_tri)) {
            L_dir += eval_area_light(wx) 
                * std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval()
                * std::max(0.0f, dot(-wx, tri.face_normal)) / p_wi_BRDF * w_BRDF;
        }
    }

    // Microfacet material
    else if (std::holds_alternative<Microfacet>(BoxScene::materials[tri.material_id])) {
        // Randomly choose one direction wi with BRDF sampling
        auto [wi, p_wi_BRDF] = 
            std::get<Microfacet>(BoxScene::materials[tri.material_id])
                .sample(wo, tri.face_normal, UniformSampler::next2d());
         
        // Trace the new ray
        const auto [is_ray_hit, t_min, nearest_tri] = 
            RayTracer::closest_hit(p, wi, octree, BoxScene::triangles);
        const Vec3f q = p + t_min * wi;

        // If the ray hit a non-emitting object at q
        if (is_ray_hit && !is_emitter(nearest_tri)) {
            L_indir = shade(nearest_tri, offset_ray_origin(q, nearest_tri.face_normal), -wi)
                * std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal)
                * std::max(0.0f, dot(wi, tri.face_normal)) / (p_wi_BRDF * p_rr);
        }

        // If the ray hit an emitting object at q
        else if (is_ray_hit && is_emitter(nearest_tri)) {
            float p_wi_light = BoxScene::inv_light_area * std::pow(distance(p,q), 2.0) / dot(wi, n_x);
            float w_BRDF = p_wi_BRDF / (p_wi_light + p_wi_BRDF);
            L_dir += eval_area_light(wx) 
                * std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal)
                * std::max(0.0f,dot(-wx, tri.face_normal)) / (p_wi_BRDF * p_rr) * w_BRDF;
        }
    }

    // Dielectric material
    else if (std::holds_alternative<Dielectric>(BoxScene::materials[tri.material_id])) {
        // Randomly choose one direction wi with BRDF sampling
        auto [wi, p_wi_BRDF, did_reflect] = 
            std::get<Dielectric>(BoxScene::materials[tri.material_id])
                .sample(wo, tri.face_normal, UniformSampler::next2d());
         
        // Trace the new ray
        const auto [is_ray_hit, t_min, nearest_tri] = 
            RayTracer::closest_hit(p, wi, octree, BoxScene::triangles);
        const Vec3f q = p + t_min * wi;

        // If the reflected ray hit a non-emitting object at q
        if (is_ray_hit && !is_emitter(nearest_tri) && did_reflect) {
            L_indir = shade(nearest_tri, offset_ray_origin(q, nearest_tri.face_normal), -wi)
                * std::get<Dielectric>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal) // eval_r
                * std::abs(dot(wi, tri.face_normal)) / (p_wi_BRDF * p_rr);
        }

        // If the reflected ray hit an emitting object at q
        else if (is_ray_hit && is_emitter(nearest_tri)) {
            float p_wi_light = BoxScene::inv_light_area * std::pow(distance(p,q), 2.0) / dot(wi, n_x);
            float w_BRDF = p_wi_BRDF / (p_wi_light + p_wi_BRDF);
            L_dir += eval_area_light(wx) 
                * std::get<Dielectric>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal)
                * std::max(0.0f,dot(-wx, tri.face_normal)) / (p_wi_BRDF * p_rr) * w_BRDF;
        }

        // If the transmitted ray hit an non-emitting object at q
        if (is_ray_hit && !is_emitter(nearest_tri) && !did_reflect) {
            L_indir = shade(nearest_tri, offset_ray_origin(q, nearest_tri.face_normal), -wi)
                * std::get<Dielectric>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal) // eval_t
                * std::abs(dot(wi, tri.face_normal)) / (p_wi_BRDF * p_rr);
        }

        // If the transmitted ray hit an emitting object at q
        else if (is_ray_hit && is_emitter(nearest_tri)) {
            float p_wi_light = BoxScene::inv_light_area * std::pow(distance(p,q), 2.0) / dot(wi, n_x);
            float w_BRDF = p_wi_BRDF / (p_wi_light + p_wi_BRDF);
            L_dir += eval_area_light(wx) 
                * std::get<Dielectric>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal)
                * std::max(0.0f,dot(-wx, tri.face_normal)) / (p_wi_BRDF * p_rr) * w_BRDF;
        }
    }

    return L_dir + L_indir;
}

Vec3f path_tracing(Vec3f ray_pos, Vec3f ray_dir) {
    const auto [is_ray_hit, t_min, nearest_tri] =
        RayTracer::closest_hit(ray_pos, ray_dir, octree, BoxScene::triangles);
    if (!is_ray_hit) return Vec3f{0.0f};
    const Vec3f hit_position = ray_pos + t_min * ray_dir;
    if (is_emitter(nearest_tri)) return eval_area_light(-ray_dir);

    return shade(nearest_tri, hit_position, -ray_dir);
}

void add_to_scene(std::string obj_path, int material_id, Vec3f emission = Vec3f{0.0f}) {
    std::vector<Triangle> obj_triangles = load_obj(obj_path, material_id, emission);
    BoxScene::triangles.insert(BoxScene::triangles.end(),
                               std::make_move_iterator(obj_triangles.begin()),
                               std::make_move_iterator(obj_triangles.end()));
}

int main(int argc, char **argv) {


    spdlog::info("\n"
                 "----------------------------------------------\n"
                 "Welcome to CS 190I Assignment 4: Microfacet Materials\n"
                 "----------------------------------------------");
    const unsigned int max_spp = 256;
    const unsigned int image_width = 1080; 
    const unsigned int image_height = 1080; 
    // Some prepereations
    Image image{.width = image_width,
                .height = image_height,
                .pixels = std::vector<Vec3f>(image_width * image_height)};
    Camera camera{.vertical_field_of_view = 68.6f,
                  .aspect = static_cast<float>(image_width) / image_height,
                  .focal_distance = 0.8f,
                  .position = Vec3f{0.278f, 0.5f, 0.2744f},
                  .view_direction = Vec3f{0.0f, -1.0f, 0.0f},
                  .up_direction = Vec3f{0.0f, 0.0f, 1.0f},
                  .right_direction = Vec3f{-1.0f, 0.0f, 0.0f}};
    camera.init();
    UniformSampler::init(190);


    // Change the material ID of the object
    // Diffuse
    // const int obj_material_id = 0;
    // Iron
    // const int obj_material_id = 5;
    // Gold
    // const int obj_material_id = 6;
    // Smooth Glass
    const int smooth_glass = 7;
    // Rough Glass
    const int rough_tinted_glass = 8;
    
    // Load the scene
    // If program can't find the bunny.obj file, use xmake run -w . or move the bunny.obj file to the 
    // same directory as the executable file.
    add_to_scene("./objects/baby_bunny.obj", rough_tinted_glass);
    add_to_scene("./objects/mama_bunny.obj", smooth_glass);

    octree.build_octree(BoxScene::triangles);
    
    // Path Tracing with MIS with GGX
    spdlog::info("Rendering started!");
    for (int y = 0; y < image.height; y++) {
        if (y % 50 == 0) {
            spdlog::info("Rendering row {} / {} \r", y, image.height);
        }
        #pragma omp parallel for
        for (int x = 0; x < image.width; x++) {
            image(x, y) = Vec3f{0.0f};
            #pragma omp parallel for
            for (int sample = 0; sample < max_spp; sample++) {
                const float u = (x + UniformSampler::next1d()) / image.width;
                const float v = (y + UniformSampler::next1d()) / image.height;
                Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
                image(x, y) +=
                    clamp(path_tracing(camera.position, ray_direction),
                          Vec3f(0.0f), Vec3f(50.0f));
            }
            image(x, y) /= (float)max_spp;
        }
    }
    spdlog::info("Rendering finished!");
    image.save_with_tonemapping("./output/image.png");

    return 0;
}