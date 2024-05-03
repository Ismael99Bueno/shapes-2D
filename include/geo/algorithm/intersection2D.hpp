#pragma once

#include "geo/shapes2D/circle.hpp"
#include "geo/shapes2D/polygon.hpp"
#include "geo/shapes2D/aabb2D.hpp"
#include "geo/algorithm/ray2D.hpp"
#include <glm/vec2.hpp>
#include <array>
#include <utility>

namespace geo
{
struct gjk_result2D
{
    bool intersect;
    std::array<glm::vec2, 3> simplex;
};
struct mtv_result2D
{
    bool valid;
    glm::vec2 mtv;
};
struct contact_feature2D
{
    enum class type
    {
        VERTEX = 0,
        FACE = 1
    };
    contact_feature2D() = default;
    contact_feature2D(std::size_t index1, std::size_t index2, type type1, type type2, bool flipped);

    std::uint8_t index1;
    std::uint8_t index2;
    std::uint8_t type1;
    std::uint8_t type2;
};
union contact_id2D {
    contact_feature2D feature{};
    std::uint32_t key;
};
struct contact_point2D
{
    glm::vec2 point;
    contact_id2D id{};
    float penetration;
};

gjk_result2D gjk(const shape2D &sh1, const shape2D &sh2);
mtv_result2D epa(const shape2D &sh1, const shape2D &sh2, const std::array<glm::vec2, 3> &simplex,
                 float threshold = 1.e-3f);

bool may_intersect(const shape2D &sh1, const shape2D &sh2);

bool intersects(const aabb2D &bb1, const aabb2D &bb2);
bool intersects(const aabb2D &bb, const glm::vec2 &point);
bool intersects(const circle &c1, const circle &c2);
ray2D::hit intersects(const aabb2D &bb, const ray2D &ray);
ray2D::hit intersects(const circle &circ, const ray2D &ray);
template <std::size_t Capacity> ray2D::hit intersects(const polygon<Capacity> &poly, const ray2D &ray)
{
    const glm::vec2 &origin = ray.origin();
    const auto &globals = poly.vertices.globals;
    float dot1 = glm::dot(ray.normal(), globals[0] - origin);
    for (std::size_t i = 1; i < poly.vertices.size(); i++)
    {
        const glm::vec2 offset = globals[i] - origin;
        const float dot2 = glm::dot(ray.normal(), offset);
        const glm::vec2 &normal = poly.vertices.normals[i - 1];
        if (dot1 * dot2 < 0.f && glm::dot(ray.direction(), offset) > 0.f && glm::dot(ray.direction(), normal) < 0.f)
        {
            const glm::vec2 &edge = poly.vertices.edges[i - 1];
            const float interp = -glm::dot(offset, edge) / glm::length2(edge);
            const glm::vec2 point = globals[i - 1] + interp * edge;
            const float distance = glm::distance(origin, point);
            if (ray.infinite() || ray.length() >= distance)
                return {point, normal, distance, true};
        }
        dot1 = dot2;
    }
    ray2D::hit hit;
    hit.hit = false;
    return hit;
}

mtv_result2D mtv(const circle &c1, const circle &c2);

contact_point2D mtv_support_contact_point(const shape2D &sh1, const shape2D &sh2, const glm::vec2 &mtv);
contact_point2D radius_distance_contact_point(const circle &c1, const circle &c2, const glm::vec2 &mtv);
contact_point2D radius_penetration_contact_point(const circle &circ, const glm::vec2 &mtv);

template <std::size_t Capacity>
kit::dynarray<contact_point2D, 2> clipping_contacts(const polygon<Capacity> &poly1, const polygon<Capacity> &poly2,
                                                    const glm::vec2 &mtv)
{
    float max_dot = 1.05f * glm::dot(mtv, poly1.vertices.normals[0]);
    std::size_t ref_index = 0;

    const polygon<Capacity> *ref_poly = &poly1;
    const polygon<Capacity> *inc_poly = &poly2;

    for (std::size_t i = 1; i < poly1.vertices.size(); i++)
    {
        const float dot = 1.05f * glm::dot(mtv, poly1.vertices.normals[i]);
        if (dot > max_dot)
        {
            max_dot = dot;
            ref_index = i;
        }
    }
    bool flipped = false;
    for (std::size_t i = 0; i < poly2.vertices.size(); i++)
    {
        const float dot = 0.95f * glm::dot(-mtv, poly2.vertices.normals[i]);
        if (dot > max_dot)
        {
            max_dot = dot;
            ref_index = i;
            ref_poly = &poly2;
            inc_poly = &poly1;
            flipped = true;
        }
    }

    const glm::vec2 &ref_normal = ref_poly->vertices.normals[ref_index];
    float min_dot = glm::dot(ref_normal, inc_poly->vertices.normals[0]);
    std::size_t inc_index = 0;
    for (std::size_t i = 1; i < inc_poly->vertices.size(); i++)
    {
        const float dot = glm::dot(ref_normal, inc_poly->vertices.normals[i]);
        if (dot < min_dot)
        {
            min_dot = dot;
            inc_index = i;
        }
    }

    const auto clip_contact = [ref_poly, flipped](const kit::dynarray<contact_point2D, 2> &unclipped,
                                                  const std::size_t ridx, const glm::vec2 &ref_tangent) {
        kit::dynarray<contact_point2D, 2> clipped;

        const glm::vec2 &rv = ref_poly->vertices.globals[ridx];
        const glm::vec2 &iv1 = unclipped[0].point;
        const glm::vec2 &iv2 = unclipped[1].point;

        const float side_pntr1 = glm::dot(ref_tangent, iv1 - rv);
        const float side_pntr2 = glm::dot(ref_tangent, iv2 - rv);
        if (side_pntr1 < 0.f)
            clipped.push_back(unclipped[0]);
        if (side_pntr2 < 0.f)
            clipped.push_back(unclipped[1]);

        if (side_pntr1 * side_pntr2 < 0.f)
        {
            const float lerp = side_pntr1 / (side_pntr1 - side_pntr2);
            const glm::vec2 point = iv1 + lerp * (iv2 - iv1);
            const contact_feature2D cf{ridx, unclipped[0].id.feature.index1, contact_feature2D::type::VERTEX,
                                       contact_feature2D::type::FACE, flipped};
            clipped.push_back({point, {cf}, 0.f});
        }
        return clipped;
    };

    const std::size_t iidx1 = inc_index;
    const std::size_t iidx2 = (inc_index + 1) % inc_poly->vertices.size();

    const std::size_t ridx1 = ref_index;
    const std::size_t ridx2 = (ref_index + 1) % ref_poly->vertices.size();

    const glm::vec2 ref_tangent = glm::vec2{-ref_normal.y, ref_normal.x};

    const contact_feature2D cf1{ridx1, iidx1, contact_feature2D::type::FACE, contact_feature2D::type::VERTEX, flipped};
    const contact_feature2D cf2{ridx2, iidx2, contact_feature2D::type::FACE, contact_feature2D::type::VERTEX, flipped};
    const kit::dynarray<contact_point2D, 2> unclipped{{inc_poly->vertices.globals[iidx1], {cf1}, 0.f},
                                                      {inc_poly->vertices.globals[iidx2], {cf2}, 0.f}};

    kit::dynarray<contact_point2D, 2> clipped = clip_contact(unclipped, ridx1, -ref_tangent);
    if (clipped.size() != 2)
    {
        clipped.clear();
        return clipped;
    }
    clipped = clip_contact(clipped, ridx2, ref_tangent);
    if (clipped.size() != 2)
    {
        clipped.clear();
        return clipped;
    }

    for (auto it = clipped.begin(); it != clipped.end();)
    {
        const float front_pntr = glm::dot(ref_normal, it->point - ref_poly->vertices.globals[ridx1]);
        if (front_pntr >= 0.f)
            it = clipped.erase(it);
        else
        {
            it->penetration = front_pntr;
            ++it;
        }
    }
    return clipped;
}
} // namespace geo
