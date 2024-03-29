#pragma once

#include "geo/shapes2D/circle.hpp"
#include "geo/shapes2D/polygon.hpp"
#include "geo/shapes2D/aabb2D.hpp"
#include <glm/vec2.hpp>
#include <array>
#include <utility>

namespace geo
{
struct gjk_result
{
    bool intersect;
    std::array<glm::vec2, 3> simplex;
};
struct mtv_result
{
    bool valid;
    glm::vec2 mtv;
};

gjk_result gjk(const shape2D &sh1, const shape2D &sh2);
mtv_result epa(const shape2D &sh1, const shape2D &sh2, const std::array<glm::vec2, 3> &simplex,
               float threshold = 1.e-3f);

glm::vec2 mtv_support_contact_point(const shape2D &sh1, const shape2D &sh2, const glm::vec2 &mtv);
bool may_intersect(const shape2D &sh1, const shape2D &sh2);

bool intersects(const aabb2D &bb1, const aabb2D &bb2);
bool intersects(const aabb2D &bb, const glm::vec2 &point);
bool intersects(const circle &c1, const circle &c2);

mtv_result mtv(const circle &c1, const circle &c2);
glm::vec2 radius_distance_contact_point(const circle &c1, const circle &c2);

template <std::size_t MaxPoints, std::size_t Capacity>
kit::dynarray<glm::vec2, MaxPoints> clipping_contacts(const polygon<Capacity> &poly1, const polygon<Capacity> &poly2,
                                                      const glm::vec2 &mtv, bool include_intersections = true)
{
    float max_dot = glm::dot(mtv, poly1.vertices.normals[0]);
    std::size_t normal_index = 0;

    const polygon<Capacity> *ref_poly = &poly1;
    const polygon<Capacity> *inc_poly = &poly2;

    for (std::size_t i = 1; i < poly1.vertices.size(); i++)
    {
        const float dot = glm::dot(mtv, poly1.vertices.normals[i]);
        if (dot > max_dot)
        {
            max_dot = dot;
            normal_index = i;
        }
    }
    for (std::size_t i = 0; i < poly2.vertices.size(); i++)
    {
        const float dot = glm::dot(-mtv, poly2.vertices.normals[i]);
        if (dot > max_dot)
        {
            max_dot = dot;
            normal_index = i;
            ref_poly = &poly2;
            inc_poly = &poly1;
        }
    }
    const glm::vec2 &normal = ref_poly->vertices.normals[normal_index];
    const glm::vec2 &start = ref_poly->vertices.globals[normal_index];
    const auto &inc_globals = inc_poly->vertices.globals;

    kit::dynarray<glm::vec2, MaxPoints> result;
    float current_dot = glm::dot(inc_globals[0] - start, normal);
    for (std::size_t i = 0; i < inc_poly->vertices.size(); i++)
    {
        const float next_dot = glm::dot(inc_globals[i + 1] - start, normal);
        if (current_dot <= 0.f)
        {
            result.push_back(inc_globals[i]);
            if (result.size() == MaxPoints)
                break;
        }
        if (include_intersections && current_dot * next_dot < 0.f)
        {
            const float current_abs = abs(current_dot);
            const float next_abs = abs(next_dot);
            result.push_back(inc_globals[i] +
                             (inc_globals[i + 1] - inc_globals[i]) * current_abs / (current_abs + next_abs));

            if (result.size() == MaxPoints)
                break;
        }

        current_dot = next_dot;
    }
    if (ref_poly != &poly1)
        for (std::size_t i = 0; i < result.size(); i++)
            result[i] += mtv;
    return result;
}
} // namespace geo
