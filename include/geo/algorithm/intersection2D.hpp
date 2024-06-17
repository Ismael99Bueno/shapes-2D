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
    operator bool() const;
};

struct mtv_result2D
{
    bool valid;
    glm::vec2 mtv;
    operator bool() const;
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
bool intersects(const aabb2D &bb, const ray2D &ray);

template <typename T> ray2D::hit<T> intersects(const circle &circ, const ray2D &ray, T *object = nullptr)
{
    const glm::vec2 &origin = ray.origin();
    const glm::vec2 &dir = ray.direction();
    const glm::vec2 &center = circ.gcentroid();
    const glm::vec2 OC = center - origin;
    const float R2 = circ.radius() * circ.radius();
    if (glm::dot(OC, dir) <= 0.f || glm::distance2(center, origin) <= R2)
        return {};
    if (!ray.infinite())
    {
        const glm::vec2 end = origin + dir * ray.length();
        const glm::vec2 EC = center - end;
        if (glm::dot(EC, dir) >= 0.f && glm::distance2(center, end) >= R2)
            return {};
    }

    const glm::vec2 &normal = ray.normal();
    const float pdist = -glm::dot(OC, normal);
    const glm::vec2 proj = center + pdist * normal;
    const glm::vec2 point = proj - dir * glm::sqrt(R2 - pdist * pdist);
    return {point, glm::normalize(point - center), glm::distance(origin, point), true, object};
}
template <typename T, std::size_t Capacity>
ray2D::hit<T> intersects(const polygon<Capacity> &poly, const ray2D &ray, T *object = nullptr)
{
    const glm::vec2 &origin = ray.origin();
    const glm::vec2 &normal = ray.normal();
    const glm::vec2 &dir = ray.direction();
    const auto &globals = poly.vertices.globals;
    for (std::size_t i = 0; i < poly.vertices.size(); i++)
    {
        const glm::vec2 &pnormal = poly.vertices.normals[i];
        const float perpend = glm::dot(pnormal, dir);
        if (perpend >= -1.e-6f)
            continue;

        const glm::vec2 &v = globals[i];
        const glm::vec2 &edge = poly.vertices.edges[i];
        const glm::vec2 offset = v - origin;

        const float u = -glm::dot(offset, normal) / glm::dot(edge, normal);
        if (u <= 0.f || u >= 1.f)
            continue;
        const float distance = glm::dot(offset, pnormal) / perpend;
        if (!ray.infinite() && distance > ray.length())
            continue;
        const glm::vec2 point = origin + distance * dir;
        return {point, pnormal, distance, true, object};
    }
    return {};
}

mtv_result2D mtv(const circle &c1, const circle &c2);

contact_point2D mtv_support_contact_point(const shape2D &sh1, const shape2D &sh2, const glm::vec2 &mtv);
contact_point2D radius_distance_contact_point(const circle &c1, const circle &c2, const glm::vec2 &mtv);
contact_point2D radius_penetration_contact_point(const circle &circ, const glm::vec2 &mtv);

template <std::size_t Capacity>
kit::dynarray<contact_point2D, 2> clipping_contacts(const polygon<Capacity> &poly1, const polygon<Capacity> &poly2,
                                                    const glm::vec2 &mtv)
{
    float max_dot = glm::dot(mtv, poly1.vertices.normals[0]);
    std::size_t ref_index = 0;

    const polygon<Capacity> *ref_poly = &poly1;
    const polygon<Capacity> *inc_poly = &poly2;

    for (std::size_t i = 1; i < poly1.vertices.size(); i++)
    {
        const float dot = glm::dot(mtv, poly1.vertices.normals[i]);
        if (dot > max_dot)
        {
            max_dot = dot;
            ref_index = i;
        }
    }
    bool flipped = false;
    for (std::size_t i = 0; i < poly2.vertices.size(); i++)
    {
        const float dot = glm::dot(-mtv, poly2.vertices.normals[i]);
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

    return clipping_contacts(*ref_poly, *inc_poly, ref_index, inc_index, flipped);
}

template <std::size_t Capacity>
kit::dynarray<contact_point2D, 2> clipping_contacts(const polygon<Capacity> &reference,
                                                    const polygon<Capacity> &incident, const std::size_t ref_index,
                                                    const std::size_t inc_index, const bool flipped = false)
{
    const auto clip_contact = [reference, flipped](const kit::dynarray<contact_point2D, 2> &unclipped,
                                                   const std::size_t ridx, const glm::vec2 &ref_tangent) {
        kit::dynarray<contact_point2D, 2> clipped;

        const glm::vec2 &rv = reference.vertices.globals[ridx];
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
    const std::size_t iidx2 = (inc_index + 1) % incident.vertices.size();

    const std::size_t ridx1 = ref_index;
    const std::size_t ridx2 = (ref_index + 1) % reference.vertices.size();

    const glm::vec2 &ref_normal = reference.vertices.normals[ridx1];
    const glm::vec2 ref_tangent = glm::vec2{-ref_normal.y, ref_normal.x};

    const contact_feature2D cf1{ridx1, iidx1, contact_feature2D::type::FACE, contact_feature2D::type::VERTEX, flipped};
    const contact_feature2D cf2{ridx2, iidx2, contact_feature2D::type::FACE, contact_feature2D::type::VERTEX, flipped};
    const kit::dynarray<contact_point2D, 2> unclipped{{incident.vertices.globals[iidx1], {cf1}, 0.f},
                                                      {incident.vertices.globals[iidx2], {cf2}, 0.f}};

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
        const float front_pntr = glm::dot(ref_normal, it->point - reference.vertices.globals[ridx1]);
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
