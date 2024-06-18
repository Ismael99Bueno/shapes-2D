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
    bool intersects = false;
    kit::dynarray<glm::vec2, 3> simplex;
    operator bool() const;
};

struct mtv_result2D
{
    bool valid = false;
    glm::vec2 mtv{0.f};
    operator bool() const;
};

struct sat_result2D
{
    bool intersects = false;
    glm::vec2 mtv{0.f};
    bool reference_incident_flipped = false;
    std::size_t reference_index = SIZE_MAX;
    std::size_t incident_index = SIZE_MAX;
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

glm::vec2 triple_cross(const glm::vec2 &v1, const glm::vec2 &v2, const glm::vec2 &v3);
void gjk_line_case(const kit::dynarray<glm::vec2, 3> &simplex, glm::vec2 &dir);
bool gjk_triangle_case(kit::dynarray<glm::vec2, 3> &simplex, glm::vec2 &dir);

template <Shape2D S1, Shape2D S2> gjk_result2D gjk(const S1 &sh1, const S2 &sh2)
{
    KIT_PERF_FUNCTION()
    if constexpr (std::is_same_v<S1, circle> && std::is_same_v<S2, circle>)
    {
        KIT_WARN("Using gjk algorithm to check if two circles are intersecting is overkill")
    }
    if constexpr (std::is_same_v<S1, shape2D> && std::is_same_v<S2, shape2D>)
    {
        KIT_ASSERT_WARN(!dynamic_cast<const circle *>(&sh1) || !dynamic_cast<const circle *>(&sh2),
                        "Using gjk algorithm to check if two circles are intersecting is overkill")
    }

    gjk_result2D result{};

    glm::vec2 dir = sh2.gcentroid() - sh1.gcentroid();
    const glm::vec2 supp = sh1.support_point(dir) - sh2.support_point(-dir);
    result.simplex.push_back(supp);
    dir = -supp;

    for (;;)
    {
        const glm::vec2 A = sh1.support_point(dir) - sh2.support_point(-dir);
        if (glm::dot(A, dir) <= 0.f)
            return result;

        result.simplex.push_back(A);
        if (result.simplex.size() == 2)
            gjk_line_case(result.simplex, dir);
        else if (gjk_triangle_case(result.simplex, dir))
        {
            result.intersects = true;
            return result;
        }
    }
}

template <Shape2D S1, Shape2D S2>
mtv_result2D epa(const S1 &sh1, const S2 &sh2, const kit::dynarray<glm::vec2, 3> &simplex,
                 const float threshold = 1.e-3f)
{
    KIT_ASSERT_ERROR(threshold > 0.f, "EPA Threshold must be greater than 0: {0}", threshold)
    KIT_PERF_FUNCTION()

    kit::dynarray<glm::vec2, 12> hull{simplex.begin(), simplex.end()};

    float min_dist = FLT_MAX;
    mtv_result2D result{};
    for (;;)
    {
        std::size_t min_index;
        for (std::size_t i = 0; i < hull.size(); i++)
        {
            const std::size_t j = (i + 1) % hull.size();

            const glm::vec2 &p1 = hull[i];
            const glm::vec2 &p2 = hull[j];

            const glm::vec2 edge = p2 - p1;

            glm::vec2 normal = glm::normalize(glm::vec2(edge.y, -edge.x));
            float dist = glm::dot(normal, p1);
            if (dist < 0.f)
            {
                dist *= -1.f;
                normal *= -1.f;
            }
            if (dist < min_dist)
            {
                min_dist = dist;
                min_index = j;
                result.mtv = normal;
            }
        }
        if (kit::approaches_zero(glm::length2(result.mtv)))
            return result;

        const glm::vec2 support = sh1.support_point(result.mtv) - sh2.support_point(-result.mtv);
        const float sup_dist = glm::dot(result.mtv, support);
        const float diff = std::abs(sup_dist - min_dist);
        if (diff <= threshold || hull.size() == hull.capacity())
            break;
        hull.insert(hull.begin() + min_index, support);
        min_dist = FLT_MAX;
    }

    result.mtv *= min_dist;
    if (kit::approaches_zero(glm::length2(result.mtv)))
        return result;

    result.valid = true;
    return result;
}

template <std::size_t Capacity> glm::vec2 sat_project_polygon(const polygon<Capacity> &poly, const glm::vec2 &axis)
{
    float proj = glm::dot(poly.vertices.globals[0], axis);
    float mm = proj;
    float mx = proj;
    for (std::size_t i = 1; i < poly.vertices.size(); i++)
    {
        proj = glm::dot(poly.vertices.globals[i], axis);
        if (proj < mm)
            mm = proj;
        if (proj > mx)
            mx = proj;
    }
    return {mm, mx};
}

glm::vec2 sat_project_circle(const circle &circ, const glm::vec2 &axis);

template <std::size_t Capacity> sat_result2D sat(const polygon<Capacity> &poly1, const polygon<Capacity> &poly2)
{
    sat_result2D result{};

    std::size_t index1 = SIZE_MAX;
    float min_overlap = FLT_MAX;
    for (std::size_t i = 0; i < poly1.vertices.size(); i++)
    {
        const glm::vec2 &normal = poly1.vertices.normals[i];
        const glm::vec2 range1 = sat_project_polygon(poly1, normal);
        const glm::vec2 range2 = sat_project_polygon(poly2, normal);
        const float overlap = glm::min(range1.y, range2.y) - glm::max(range1.x, range2.x);
        if (overlap <= 0.f)
            return result;
        if (overlap < min_overlap)
        {
            index1 = i;
            min_overlap = overlap;
        }
    }

    std::size_t index2 = SIZE_MAX;
    bool flipped = false;
    for (std::size_t i = 0; i < poly2.vertices.size(); i++)
    {
        const glm::vec2 &normal = poly2.vertices.normals[i];
        const glm::vec2 range1 = project_polygon(poly1, normal);
        const glm::vec2 range2 = project_polygon(poly2, normal);
        const float overlap = glm::min(range1.y, range2.y) - glm::max(range1.x, range2.x);
        if (overlap <= 0.f)
            return result;
        if (overlap < min_overlap)
        {
            index2 = i;
            min_overlap = overlap;
            flipped = true;
        }
    }

    result.intersects = true;
    result.reference_index = flipped ? index2 : index1;
    result.incident_index = flipped ? index1 : index2;
    result.mtv = min_overlap * (flipped ? -poly2.vertices.normals[index2] : poly1.vertices.normals[index1]);
    result.reference_incident_flipped = flipped;
    return result;
}
template <std::size_t Capacity> sat_result2D sat(const polygon<Capacity> &poly, const circle &circ)
{
    sat_result2D result{};
    const auto closest_point_on_segment = [](const glm::vec2 &p, const glm::vec2 &a, const glm::vec2 &b,
                                             const glm::vec2 &ab) {
        const glm::vec2 ap = p - a;
        const float t = glm::dot(ap, ab) / glm::length2(ab);
        if (t <= 0.f)
            return a;
        if (t >= 1.f)
            return b;
        return a + t * ab;
    };

    glm::vec2 closest_to_center{FLT_MAX};
    float min_dist = FLT_MAX;

    std::size_t index = SIZE_MAX;
    float min_overlap = FLT_MAX;
    for (std::size_t i = 0; i < poly.vertices.size(); i++)
    {
        const glm::vec2 &normal = poly.vertices.normals[i];
        const glm::vec2 range1 = sat_project_polygon(poly, normal);
        const glm::vec2 range2 = sat_project_circle(circ, normal);
        const float overlap = glm::min(range1.y, range2.y) - glm::max(range1.x, range2.x);
        if (overlap <= 0.f)
            return result;
        if (overlap < min_overlap)
        {
            min_overlap = overlap;
            index = i;
        }

        const glm::vec2 closest = closest_point_on_segment(circ.gcentroid(), poly.vertices.globals[i],
                                                           poly.vertices.globals[i + 1], poly.vertices.edges[i]);
        const float dist = glm::distance2(circ.gcentroid(), closest);
        if (dist < min_dist)
        {
            min_dist = dist;
            closest_to_center = closest;
        }
    }
    const glm::vec2 last_axis = glm::normalize(circ.gcentroid() - closest_to_center);
    const glm::vec2 range1 = sat_project_polygon(poly, last_axis);
    const glm::vec2 range2 = sat_project_circle(circ, last_axis);
    const float overlap = glm::min(range1.y, range2.y) - glm::max(range1.x, range2.x);
    if (overlap <= 0.f)
        return result;

    result.intersects = true;
    if (overlap < min_overlap)
    {
        result.incident_index = index;
        result.mtv = overlap * last_axis;
        result.reference_incident_flipped = true;
        return result;
    }
    result.reference_index = index;
    result.mtv = min_overlap * poly.vertices.normals[index];
    return result;
}
template <std::size_t Capacity> sat_result2D sat(const circle &circ, const polygon<Capacity> &poly)
{
    sat_result2D result = sat(poly, circ);
    result.reference_incident_flipped = !result.reference_incident_flipped;
    result.mtv = -result.mtv;
    return result;
}

bool intersects(const aabb2D &bb1, const aabb2D &bb2);
bool intersects(const aabb2D &bb, const glm::vec2 &point);
bool intersects(const aabb2D &bb, const ray2D &ray);

bool intersects(const circle &c1, const circle &c2);
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
    KIT_ASSERT_ERROR(ref_index < reference.vertices.size(), "Reference index out of bounds")
    KIT_ASSERT_ERROR(inc_index < incident.vertices.size(), "Incident index out of bounds")

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

template <std::size_t Capacity>
kit::dynarray<contact_point2D, 2> clipping_contacts(const polygon<Capacity> &reference,
                                                    const polygon<Capacity> &incident, const sat_result2D &sat_result)
{
    KIT_ASSERT_ERROR(sat_result.intersects, "SAT Result must report intersection")
    KIT_ASSERT_ERROR(sat_result.reference_index != SIZE_MAX && sat_result.incident_index != SIZE_MAX,
                     "SAT Result indices must both be valid. If one of them is SIZE_MAX, this SAT result may "
                     "correspond to a circle-polygon intersection")
}

} // namespace geo
