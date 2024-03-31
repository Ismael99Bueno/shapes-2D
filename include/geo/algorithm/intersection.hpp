#pragma once

#include "geo/shapes2D/circle.hpp"
#include "geo/shapes2D/polygon.hpp"
#include "geo/shapes2D/aabb2D.hpp"
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

gjk_result2D gjk(const shape2D &sh1, const shape2D &sh2);
mtv_result2D epa(const shape2D &sh1, const shape2D &sh2, const std::array<glm::vec2, 3> &simplex,
                 float threshold = 1.e-3f);

glm::vec2 mtv_support_contact_point(const shape2D &sh1, const shape2D &sh2, const glm::vec2 &mtv);
bool may_intersect(const shape2D &sh1, const shape2D &sh2);

bool intersects(const aabb2D &bb1, const aabb2D &bb2);
bool intersects(const aabb2D &bb, const glm::vec2 &point);
bool intersects(const circle &c1, const circle &c2);

mtv_result2D mtv(const circle &c1, const circle &c2);

glm::vec2 radius_distance_contact_point(const circle &c1, const circle &c2);
glm::vec2 radius_penetration_contact_point(const circle &circ, const glm::vec2 &mtv);

struct contact_feature
{
    enum class type
    {
        VERTEX = 0,
        FACE = 1
    };
    contact_feature() = default;
    contact_feature(std::size_t index1, std::size_t index2, type type1, type type2, bool flipped);

    std::uint8_t index1 = 0;
    std::uint8_t index2 = 0;
    std::uint8_t type1 = 0;
    std::uint8_t type2 = 0;
};

union contact_id {
    contact_feature feature{};
    std::uint32_t key;
};

struct contact_point2D
{
    glm::vec2 point;
    contact_id id{};
};

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
            const contact_feature cf{ridx, unclipped[0].id.feature.index1, contact_feature::type::VERTEX,
                                     contact_feature::type::FACE, flipped};
            clipped.push_back({point, {cf}});
        }
        return clipped;
    };

    const std::size_t iidx1 = inc_index;
    const std::size_t iidx2 = (inc_index + 1) % inc_poly->vertices.size();

    const std::size_t ridx1 = ref_index;
    const std::size_t ridx2 = (ref_index + 1) % ref_poly->vertices.size();

    const glm::vec2 ref_tangent = glm::vec2{-ref_normal.y, ref_normal.x};

    const contact_feature cf1{ridx1, iidx1, contact_feature::type::FACE, contact_feature::type::VERTEX, flipped};
    const contact_feature cf2{ridx2, iidx2, contact_feature::type::FACE, contact_feature::type::VERTEX, flipped};
    const kit::dynarray<contact_point2D, 2> unclipped{{inc_poly->vertices.globals[iidx1], {cf1}},
                                                      {inc_poly->vertices.globals[iidx2], {cf2}}};

    // must do something if clipped is empty
    kit::dynarray<contact_point2D, 2> clipped = clip_contact(unclipped, ridx1, -ref_tangent);
    clipped = clip_contact(clipped, ridx2, ref_tangent);

    for (auto it = clipped.begin(); it != clipped.end();)
    {
        const float front_pntr = glm::dot(ref_normal, it->point - ref_poly->vertices.globals[ridx1]);
        if (front_pntr >= 0.f)
            it = clipped.erase(it);
        else
            ++it;
    }
    return clipped;
}
} // namespace geo
