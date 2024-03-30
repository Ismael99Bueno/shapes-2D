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

template <std::size_t Capacity>
glm::vec2 radius_penetration_contact_point(const circle &circ, const polygon<Capacity> &poly, const glm::vec2 &mtv)
{
    return circ.gcentroid() + circ.radius() * glm::normalize(mtv) - mtv;
}
template <std::size_t Capacity>
glm::vec2 radius_penetration_contact_point(const polygon<Capacity> &poly, const circle &circ, const glm::vec2 &mtv)
{
    return circ.gcentroid() - circ.radius() * glm::normalize(mtv) + mtv;
}

struct contact_feature
{
    enum class type
    {
        VERTEX = 0,
        NORMAL = 1
    };
    std::uint8_t index1;
    std::uint8_t index2;
    std::uint8_t type1;
    std::uint8_t type2;
};

union contact_id {
    contact_feature feature;
    std::uint32_t key;
};

struct contact_point2D
{
    glm::vec2 point;
    contact_id id;
};

template <std::size_t MaxPoints, std::size_t Capacity>
kit::dynarray<contact_point2D, MaxPoints> clipping_contacts(const polygon<Capacity> &poly1,
                                                            const polygon<Capacity> &poly2, const glm::vec2 &mtv)
{
    float max_dot = glm::dot(mtv, poly1.vertices.normals[0]);
    std::size_t normal_index = 0;

    const polygon<Capacity> *ref_poly = &poly1;
    const polygon<Capacity> *inc_poly = &poly2;

    for (std::size_t i = 1; i < poly1.vertices.size(); i++)
    {
        const float dot = 1.05f * glm::dot(mtv, poly1.vertices.normals[i]);
        if (dot > max_dot)
        {
            max_dot = dot;
            normal_index = i;
        }
    }
    bool flipped = false;
    for (std::size_t i = 0; i < poly2.vertices.size(); i++)
    {
        const float dot = 0.95f * glm::dot(-mtv, poly2.vertices.normals[i]);
        if (dot > max_dot)
        {
            max_dot = dot;
            normal_index = i;
            ref_poly = &poly2;
            inc_poly = &poly1;
            flipped = true;
        }
    }
    const glm::vec2 &normal = ref_poly->vertices.normals[normal_index];
    const glm::vec2 &start = ref_poly->vertices.globals[normal_index];
    const auto &inc_globals = inc_poly->vertices.globals;

    kit::dynarray<contact_point2D, MaxPoints> result;
    for (std::size_t i = 0; i < inc_poly->vertices.size(); i++)
    {
        const float dot = glm::dot(inc_globals[i] - start, normal);
        if (dot < 0.f)
        {
            const contact_feature cf{flipped ? (std::uint8_t)i : (std::uint8_t)normal_index,
                                     flipped ? (std::uint8_t)normal_index : (std::uint8_t)i,
                                     (std::uint8_t)(flipped ? 0u : 1u), (std::uint8_t)(flipped ? 1u : 0u)};
            result.push_back({inc_globals[i], {cf}});
            if (result.size() == MaxPoints)
                return result;
        }
    }
    return result;
}
} // namespace geo
