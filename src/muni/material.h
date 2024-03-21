#pragma once
#include "common.h"
#include "math_helpers.h"
#include <algorithm>
#include <cmath>

namespace muni {

struct Lambertian {
  Vec3f albedo;

  /** Evaluates the BRDF for the Lambertian material.
    \return The BRDF (fr) value.
  */
  Vec3f eval() const { return albedo * INV_PI; }

  /** Samples the BRDF for the Lambertian material.
    \param[in] normal The normal of the surface.
    \param[in] u A random number in (0,1)^2.
    \return A tuple containing the sampled direction in world space and the PDF.
  */
  std::tuple<Vec3f, float> sample(Vec3f normal, Vec2f u) const {
    float theta = std::acos(std::sqrt(u.x));
    float phi = 2 * M_PI * u.y;
    float z = std::cos(theta);
    float r = std::sqrt(1.0f - z * z);
    Vec3f wo = Vec3f{r * std::cos(phi), r * std::sin(phi), z};
    return std::make_tuple(from_local(wo, normal), z / M_PI);
  }

  /** Computes the PDF for the Lambertian material.
    \param[in] wo The outgoing direction in world space.
    \param[in] wi The light incident direction in world space.
    \return The PDF value.
  */
  float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
    auto [x, y] = coordinate_system(normal);
    Vec3f wo = to_local(wo_world, normal);
    return wo.z / M_PI;
  }
};

struct Microfacet {
  float roughness;
  // refraction indices for RGB channels
  Vec3f n1;
  Vec3f n2;

  /** Computes the Fresnel term for the microfacet material.
    \param[in] wi The light incident direction in local space.
    \return The Fresnel term.
  */
  Vec3f F(Vec3f wi) const {
    float cos_theta = wi.z;
    
    float R0_r = std::pow((n1.x - n2.x) / (n1.x + n2.x), 2);
    float R0_g = std::pow((n1.y - n2.y) / (n1.y + n2.y), 2);
    float R0_b = std::pow((n1.z - n2.z) / (n1.z + n2.z), 2);
    
    float F_r = R0_r + (1.0f - R0_r) * std::pow(1 - cos_theta, 5);
    float F_g = R0_g + (1.0f - R0_g) * std::pow(1 - cos_theta, 5);
    float F_b = R0_b + (1.0f - R0_b) * std::pow(1 - cos_theta, 5);

    return Vec3f{F_r, F_g, F_b};
  }

  /** Computes the GGX nomral distribution function for the microfacet material.
    \param[in] h The half vector in local space.
    \param[in] normal The normal to the surface.
    \return The normal distribution function.
  */
  float D(Vec3f h) const {
    float theta = std::acos(h.z);
    if (std::isinf(std::tan(theta)) || h.z <= 0) return 0.0f;
    float D = 1 / std::pow(1 + std::pow(std::tan(theta) / roughness, 2), 2);
    return D / (M_PI * std::pow(roughness, 2) * std::pow(h.z, 4));
  }

  float Lambda(Vec3f w) const {
    float theta = std::acos(w.z);
    
    float a;
    if (std::isinf(std::tan(theta))) a = 0.0f;
    else a = 1 / (roughness * std::tan(theta));
    return (-1 + std::sqrt(1 + 1 / std::pow(a, 2))) / 2;
  }

  /** Computes the shadowing-masking function for the microfacet material for the GGX distribution.
    \param[in] wo The outgoing direction in local space.
    \param[in] wi The light incident direction in local space.
    \return The shadowing-masking value.
  */
  float G(Vec3f wo, Vec3f wi) const { return (1 / (1 + Lambda(wo))) * (1 / (1 + Lambda(wi))); }

  /** Evaluates the BRDF for the microfacet material using GGX.
    \param[in] wo_world The outgoing direction in world space.
    \param[in] wi_world The light incident direction in world space.
    \param[in] normal The normal of the surface.
    \return The BRDF (fr) value.
  */
  Vec3f eval(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
    Vec3f wo = to_local(wo_world, normal);
    Vec3f wi = to_local(wi_world, normal);
    Vec3f wh = normalize(wo + wi);
    if (wo.z * wi.z < 0) return Vec3f{0.0f};
    return F(wi) * G(wo, wi) * D(wh)/ (4.0f * wo.z * wi.z);
  }
  
  /** Computes the PDF for the microfacet material.
    \param[in] wo The outgoing direction in world space.
    \param[in] wi The light incident direction in world space.
    \param[in] normal The normal of the surface.
    \return The PDF value.
  */
  float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
    Vec3f wo = to_local(wo_world, normal);
    Vec3f wi = to_local(wi_world, normal);
    Vec3f wh = normalize(wo + wi);
    float theta = std::acos(wh.z);

    if (theta == 0.0f) return 0.0f;

    float p_phi = 0.5f * INV_PI;
    float p_theta = std::exp(-1.0f * std::pow(std::tan(theta) / roughness, 2))
      * (2.0f * std::sin(theta)) / (std::pow(roughness, 2) * std::pow(wh.z, 3));

    float p_wh = D(wh) * std::cos(theta);
    float p_wi = p_wh / (4 * dot(wo, wh));
    // spdlog::info(" other pdf = {}", p_wi);
    return p_wi;
  }

  /** Samples the BRDF for the microfacet material for GGX.
    \param[in] wo_world The outgoing direction in world space.
    \param[in] normal The normal of the surface.
    \param[in] u A random number in (0,1)^2.
    \return A tuple containing the sampled direction in world space and the PDF.
  */
  std::tuple<Vec3f, float> sample(Vec3f wo_world, Vec3f normal, Vec2f u) const {
    float phi = 2.0f * M_PI * u.x;
    float p_phi = M_1_2PI;
    float theta = std::atan(roughness * std::sqrt(u.y) / std::sqrt(1 - u.y));
    float p_theta = std::exp(-1.0f * std::pow(std::tan(theta) / roughness, 2))
      * (2.0f * std::sin(theta)) / (std::pow(roughness, 2) * std::pow(std::cos(theta), 3));

    float x = std::sin(theta) * std::cos(phi);
    float y = std::sin(theta) * std::sin(phi);
    float z = std::cos(theta);

    Vec3f wh = normalize(Vec3f{x,y,z});
    Vec3f wo = to_local(wo_world, normal);
    Vec3f wi = reflect(-wo, wh);
    Vec3f wi_world = from_local(wi, normal);

    float p_wh = D(wh) * std::cos(theta);
    float p_wi = p_wh / (4 * dot(wo, wh));
    // spdlog::info("sample pdf = {}", p_wi);
    return {wi_world, p_wi};
  }
};

struct Dielectric {
  float roughness;
  float eta_i; // Index of Refraction of incoming medium
  float eta_t; // Index of Refraction of transmitting medium
  Vec3f attenuation; // Attenuation coefficient for absorption

  // Computes the Fresnel term for dielectrics with unpolarized light
  float F(Vec3f i, Vec3f m) const {
    float eta_i_local = eta_i;
    float eta_t_local = eta_t;
    if (i.z < 0) std::swap(eta_i_local, eta_t_local);

    float c = std::abs(dot(i, m));
    float g_squared = (eta_t * eta_t) / (eta_i * eta_i) - 1.0f + c * c;
    if (g_squared < 0) { return 1.0f; }
    else {
      float g = std::sqrt(g_squared);
      float F_value = 0.5f * std::pow(g - c, 2) / std::pow(g + c, 2) * 
        (1.0f + std::pow(c * (g + c) - 1.0f, 2) / std::pow(c * (g - c) + 1.0f, 2));
      return clamp(F_value, 0.0f, 1.0f);
    }
  }

  // Computes the GGX normal distribution function for the microfacet material.
  float D(Vec3f m) const {
    float theta = std::acos(m.z);
    if (std::isinf(std::tan(theta))) return 1.0f;
    if (m.z < 0) {
      return 0.0f;
    }
    return roughness * roughness * chi_plus(m.z) / (M_PI * std::pow(m.z, 4) * std::pow(roughness * roughness + std::tan(theta) * std::tan(theta), 2));
  }

  // Computes the shadowing-masking function for the microfacet material for the GGX distribution.
  float G(Vec3f v, Vec3f m) const { 
    float theta = std::acos(v.z);
    return chi_plus(dot(v, m) / v.z) * 2 / (1 + std::sqrt(1 + roughness * roughness * std::tan(theta) * std::tan(theta)));
  }

  // Computes the shadowing-masking function for the microfacet material for the GGX distribution.
  float G(Vec3f i, Vec3f o, Vec3f m) const { return G(i, m) * G(o, m); }

  /** Evaluates the BSDF for the dielectric microfacet material using GGX.
    \param[in] wo_world The outgoing direction in world space.
    \param[in] wi_world The light incident direction in world space.
    \param[in] normal The normal of the surface.
    \return The BRDF (fs) value.
  */
  Vec3f eval(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
    float eta_i_local = eta_i;
    float eta_t_local = eta_t;
    if (dot(wo_world, normal) < 0) std::swap(eta_i_local, eta_t_local);

    Vec3f wo = to_local(wo_world, normal);
    Vec3f wi = to_local(wi_world, normal);
    Vec3f wh_r = normalize(sign(wi.z) * normalize(wo + wi));
    Vec3f wh_t = normalize(-eta_i_local * wo - eta_t_local * wi);
    
    float fr = F(wi, wh_r) * G(wi, wo, wh_r) * D(wh_r) / (4.0f * std::abs(wi.z) * std::abs(wo.z));
    float ft = (std::abs(dot(wi, wh_t)) * std::abs(dot(wo, wh_t))) / (std::abs(wi.z) * std::abs(wo.z))
      * eta_t_local * eta_t_local * (1 - F(wi, wh_t)) * G(wi, wo, wh_t) * D(wh_t) / std::pow(eta_i_local * dot(wi, wh_t) + eta_t_local * dot(wo, wh_t), 2);

    return Vec3f{fr + ft} * attenuation;
  }

  Vec3f eval_r(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
    Vec3f wo = to_local(wo_world, normal);
    Vec3f wi = to_local(wi_world, normal);
    Vec3f wh = normalize(sign(dot(wi_world, normal)) * normalize(wo + wi));

    float eta_i_local = eta_i;
    float eta_t_local = eta_t;
    if (dot(wo_world, normal) < 0) std::swap(eta_i_local, eta_t_local);
    
    return Vec3f{F(wi, wh) / std::abs(wi.z)} * attenuation;
  }
  
  Vec3f eval_t(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
    float eta_i_local = eta_i;
    float eta_t_local = eta_t;
    if (dot(wo_world, normal) < 0) std::swap(eta_i_local, eta_t_local);

    Vec3f wo = to_local(wo_world, normal);
    Vec3f wi = to_local(wi_world, normal);
    Vec3f wh = normalize(-1.0f * normalize(eta_i_local * wi + eta_t_local * wo));
    
    return Vec3f{(1 - F(wi, wh)) / std::abs(wi.z) * (eta_i_local * eta_i_local) / (eta_t_local * eta_t_local)} * attenuation;
  }
  
  /** Computes the PDF for the dielectric microfacet material.
    \param[in] wo The outgoing direction in world space.
    \param[in] wi The light incident direction in world space.
    \param[in] normal The normal of the surface.
    \return The PDF value.
  */
  float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
    Vec3f wo = to_local(wo_world, normal);
    Vec3f wi = to_local(wi_world, normal);

    float eta_i_local = eta_i;
    float eta_t_local = eta_t;
    // if (dot(wo_world, normal) < 0) std::swap(eta_i_local, eta_t_local);

    Vec3f wh_r = sign(dot(wi_world, normal)) * normalize(wo + wi);
    Vec3f wh_t = -1.0f * normalize(eta_i_local * wi + eta_t_local * wo);
    float theta_r = std::acos(wh_r.z);
    float theta_t = std::acos(wh_t.z);

    float p_wi_refl = 0.0f;
    float p_wi_tran = 0.0f;
    float f_s_r = 0.0f;
    float f_s_t = 0.0f;

    if (theta_r == 0.0f) p_wi_refl = 0.0f;
    else {
      f_s_r = F(wo, wh_r);
      float p_wh_r = D(wh_r) * std::cos(theta_r);
      p_wi_refl = f_s_r * p_wh_r / (4 * dot(wo, wh_r));
    }
    
    if (theta_t == 0.0f) p_wi_tran = 0.0f;
    else {
      f_s_t = 1 - F(wo, wh_t);
      float p_wh_t = D(wh_t) * std::cos(theta_t);
      p_wi_tran = (1 - f_s_t) * p_wh_t * (eta_t_local * eta_t_local * dot(wo, wh_t)) / std::pow(eta_i_local * dot(wi, wh_t) + eta_t_local * dot(wo, wh_t), 2);
    }
    float pdf = f_s_r + f_s_t;
    return pdf;
  }

  /** Samples the BSDF for the dielectric microfacet material with GGX.
    \param[in] wo_world The outgoing direction in world space.
    \param[in] normal The normal of the surface.
    \param[in] u A random number in (0,1)^2.
    \return A tuple containing the sampled direction in world space and the PDF.
  */
  std::tuple<Vec3f, float, bool> sample(Vec3f wo_world, Vec3f normal, Vec2f u) const {
    float eta_i_local = eta_i;
    float eta_t_local = eta_t;
    if (dot(wo_world, normal) < 0) std::swap(eta_i_local, eta_t_local);

    float phi = 2.0f * M_PI * u.x;
    float theta = std::atan(roughness * std::sqrt(u.y) / std::sqrt(1 - u.y));

    float x = std::sin(theta) * std::cos(phi);
    float y = std::sin(theta) * std::sin(phi);
    float z = std::cos(theta);

    Vec3f m = normalize(Vec3f{x,y,z});
    float p_m = D(m) * std::cos(theta);

    Vec3f wo = to_local(wo_world, normal);
    float f_s = F(wo, m);

    float rand_num = (rand() + 0.5f) / (RAND_MAX + 1.0f); // U(0, 1)

    if (rand_num < f_s) {
      // Sample reflection
      Vec3f wi = reflect(-wo, m);
      Vec3f wi_world = from_local(wi, normal);
      float p_wi = p_m / (4 * std::abs(dot(wo, m)));
      return {wi_world, p_wi, true};
      // return {wi_world, f_s, true};
    } else {
      // Sample transmission
      Vec3f wi = refract(-wo, m, eta_i_local, eta_t_local);
      Vec3f wi_world = from_local(wi, normal);
      float p_wi = p_m * (eta_t_local * eta_t_local * std::abs(dot(wo, m))) / std::pow(eta_i_local * dot(wi, m) + eta_t_local * dot(wo, m), 2);
      return {wi_world, p_wi, false};
      // return {wi_world, 1 - f_s, false};
    }
  }
};
}  // namespace muni
