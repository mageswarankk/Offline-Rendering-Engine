#include "common.h"
#include "triangle.h"
#include "material.h"
#include <array>
#include <variant>

namespace muni { namespace BoxScene {

static const float light_x = 0.195f;
static const float light_y = -0.355f;
static const float light_z = 0.545f;
static const float light_len_x = 0.16f;
static const float light_len_y = 0.16f;
static const float inv_light_area = 1 / (light_len_x * light_len_y);
static const Vec3f light_color{50.0f * 0.678f, 50.0f * 0.439f, 50.0f * 1.0f};
static const Vec3f light_normal{0.0f, 0.0f, -1.0f};

// Microfacet materials
const Microfacet Gold{.roughness = 0.0005f, .n1 = Vec3f{1.0f}, .n2 =Vec3f{0.2177f, 0.42659f, 1.2425f}};
const Microfacet Iron{.roughness = 0.02f, .n1 = Vec3f{1.0f}, .n2 =Vec3f{2.8851f, 2.95f, 2.6f}};
const Dielectric Smooth_Glass{.roughness = 0.005f, .eta_i = 1.0f, .eta_t = 1.55f, .attenuation = Vec3f{0.95, 0.95f, 0.94}};
const Dielectric Frosted_Glass{.roughness = 0.04f, .eta_i = 1.0f, .eta_t = 1.55f, .attenuation = Vec3f{0.95, 0.95f, 0.94}};
const Dielectric Ice{.roughness = 0.0002f, .eta_i = 1.0f, .eta_t = 1.31f, .attenuation = Vec3f{0.95f, 0.95f, 0.95f}};
const Dielectric Water{.roughness = 0.0001f, .eta_i = 1.0f, .eta_t = 1.3333f, .attenuation = Vec3f{0.961f, 0.914f, 0.671f}};
static const std::array<std::variant<Lambertian, Microfacet, Dielectric>, 11> materials = {
    // Back
    Lambertian{.albedo = Vec3f{0.1f, 0.82f, 0.71f}},
    // Bottom
    Iron,
    // Left
    Lambertian{.albedo = Vec3f{0.49f, 1.0f, 0.412f}},
    // Right
    Lambertian{.albedo = Vec3f{1.0f, 0.467, 0.0f}},
    // Top
    Lambertian{.albedo = Vec3f{1.0f, 1.0f, 1.0f}},
    // Bunny
    Iron,
    Gold,
    Smooth_Glass,
    Frosted_Glass,
    Ice,
    Water
};

static std::vector<Triangle> triangles = {
    // Light
    Triangle{.v0 = Vec3f{light_x, light_y + light_len_y, light_z},
             .v1 = Vec3f{light_x + light_len_x, light_y, light_z},
             .v2 = Vec3f{light_x, light_y, light_z},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = light_color,
             .material_id = 0},
    Triangle{.v0 = Vec3f{light_x, light_y + light_len_y, light_z},
             .v1 = Vec3f{light_x + light_len_x, light_y + light_len_y, light_z},
             .v2 = Vec3f{light_x + light_len_x, light_y, light_z},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = light_color,
             .material_id = 0},
    // Back
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .v2 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 1.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 0},
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 1.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 0},
    // Bottom
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .v2 = Vec3f{0.555999935f, -0.000000119f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 0.0f, 1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 1},
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.555999935f, -0.000000119f, 0.000000040f},
             .v2 = Vec3f{0.000000133f, -0.000000119f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 0.0f, 1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 1},
    // Left
    Triangle{.v0 = Vec3f{0.555999935f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.000000119f, 0.000000040f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .face_normal = Vec3f{-1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 2},
    Triangle{.v0 = Vec3f{0.555999935f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{-1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 2},
    // Right
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.000000133f, -0.000000119f, 0.000000040f},
             .v2 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .face_normal = Vec3f{1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 3},
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .v2 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 3},
    // Top
    Triangle{.v0 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .v2 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 4},
    Triangle{.v0 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.000000119f, 0.548799932f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 4},
    };
}}  // namespace muni::BoxScene
