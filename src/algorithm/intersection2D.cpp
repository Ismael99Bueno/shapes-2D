#include "geo/internal/pch.hpp"
#include "geo/algorithm/intersection2D.hpp"
#include "geo/shapes2D/polygon.hpp"

#include "kit/utility/utils.hpp"

namespace geo
{
gjk_result2D::operator bool() const
{
    return intersects;
}
mtv_result2D::operator bool() const
{
    return valid;
}
sat_result2D::operator bool() const
{
    return intersects;
}

glm::vec2 triple_cross(const glm::vec2 &v1, const glm::vec2 &v2, const glm::vec2 &v3)
{
    const float crs = kit::cross2D(v1, v2);
    return glm::vec2(-v3.y * crs, v3.x * crs);
}

void gjk_line_case(const kit::dynarray<glm::vec2, 3> &simplex, glm::vec2 &dir)
{
    const glm::vec2 AB = simplex[0] - simplex[1], AO = -simplex[1];
    dir = triple_cross(AB, AO, AB);
}

bool gjk_triangle_case(kit::dynarray<glm::vec2, 3> &simplex, glm::vec2 &dir)
{
    const glm::vec2 AB = simplex[1] - simplex[2], AC = simplex[0] - simplex[2], AO = -simplex[2];
    const glm::vec2 ABperp = triple_cross(AC, AB, AB);
    if (glm::dot(ABperp, AO) >= 0.f)
    {
        simplex.erase(simplex.begin());
        dir = ABperp;
        return false;
    }
    const glm::vec2 ACperp = triple_cross(AB, AC, AC);
    if (glm::dot(ACperp, AO) >= 0.f)
    {
        simplex.erase(simplex.begin() + 1);
        dir = ACperp;
        return false;
    }
    return true;
}

glm::vec2 sat_project_circle(const circle &circ, const glm::vec2 &axis)
{
    KIT_PERF_FUNCTION()
    const glm::vec2 &center = circ.gcentroid();
    const float proj = glm::dot(center, axis);
    return {proj - circ.radius(), proj + circ.radius()};
}

bool intersects(const aabb2D &bb1, const aabb2D &bb2)
{
    const glm::vec2 df1 = bb2.min - bb1.max, df2 = bb1.min - bb2.max;
    if (df1.x > 0.f || df1.y > 0.f)
        return false;
    if (df2.x > 0.f || df2.y > 0.f)
        return false;
    return true;
}
bool intersects(const aabb2D &bb, const glm::vec2 &point)
{
    return point.x >= bb.min.x && point.x <= bb.max.x && point.y >= bb.min.y && point.y <= bb.max.y;
}

bool intersects(const circle &c1, const circle &c2)
{
    const float R = c1.radius() + c2.radius();
    return glm::distance2(c1.gcentroid(), c2.gcentroid()) < R * R;
}
bool intersects(const aabb2D &bb, const ray2D &ray)
{
    const glm::vec2 &origin = ray.origin();
    const glm::vec2 &dir = ray.direction();

    float tmin = 0.f;
    float tmax = ray.length();
    for (int i = 0; i < 2; i++)
    {
        const float inv_dir = 1.f / dir[i];
        float t1 = (bb.min[i] - origin[i]) * inv_dir;
        float t2 = (bb.max[i] - origin[i]) * inv_dir;
        if (t1 > t2)
            std::swap(t1, t2);
        tmin = glm::max(t1, tmin);
        tmax = glm::min(t2, tmax);
        if (tmin > tmax)
            return false;
    }

    return tmin >= 0.f && tmin <= ray.length();
}

mtv_result2D mtv(const circle &c1, const circle &c2)
{
    mtv_result2D result;

    const glm::vec2 dir = c1.gcentroid() - c2.gcentroid();
    result.mtv = dir - (c1.radius() + c2.radius()) * glm::normalize(dir);
    result.valid =
        !std::isnan(result.mtv.x) && !std::isnan(result.mtv.y) && !kit::approaches_zero(glm::length2(result.mtv));
    return result;
}

contact_point2D mtv_support_contact_point(const shape2D &sh1, const shape2D &sh2, const glm::vec2 &mtv)
{
    KIT_PERF_FUNCTION()
    const glm::vec2 sup1 = sh1.support_point(mtv), sup2 = sh2.support_point(-mtv);
    const float d1 = glm::length2(sh2.closest_direction_from(sup1 - mtv)),
                d2 = glm::length2(sh1.closest_direction_from(sup2 + mtv));

    contact_point2D contact;
    contact.id.key = 0;
    contact.penetration = -glm::length(mtv);
    contact.point = d1 < d2 ? sup1 : sup2 + mtv;
    return contact;
}
contact_point2D radius_distance_contact_point(const circle &c1, const circle &c2, const glm::vec2 &mtv)
{
    const glm::vec2 dir = glm::normalize(c1.gcentroid() - c2.gcentroid());
    contact_point2D contact;
    contact.id.key = 0;
    contact.penetration = -glm::length(mtv);
    contact.point = c1.gcentroid() - dir * c1.radius();
    return contact;
}
contact_point2D radius_penetration_contact_point(const circle &circ, const glm::vec2 &mtv)
{
    contact_point2D contact;
    contact.id.key = 0;
    contact.penetration = -glm::length(mtv);
    contact.point = circ.gcentroid() + circ.radius() * glm::normalize(mtv);
    return contact;
}

contact_feature2D::contact_feature2D(std::size_t index1, std::size_t index2, type type1, type type2, bool flipped)
    : index1((std::uint8_t)index1), index2((std::uint8_t)index2), type1((std::uint8_t)type1), type2((std::uint8_t)type2)
{
    if (flipped)
    {
        std::swap(this->index1, this->index2);
        std::swap(this->type1, this->type2);
    }
}

} // namespace geo