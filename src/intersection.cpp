#include "geo/pch.hpp"
#include "geo/intersection.hpp"
#include "geo/polygon.hpp"
#include "kit/profile/perf.hpp"
#include "kit/utility/utils.hpp"

#define EPA_EPSILON 1.e-3f

namespace geo
{
static glm::vec2 triple_cross(const glm::vec2 &v1, const glm::vec2 &v2, const glm::vec2 &v3)
{
    const float crs = kit::cross2D(v1, v2);
    return glm::vec2(-v3.y * crs, v3.x * crs);
}

struct arr3
{
    std::array<glm::vec2, 3> &data;
    std::size_t size = 0;

    void push(const glm::vec2 &val)
    {
        data[size++] = val;
        KIT_ASSERT_ERROR(size <= 3, "Array size exceeds 3!")
    }
    void erase(const std::size_t index)
    {
        KIT_ASSERT_ERROR(size > 0, "Cannot erase element of empty array!")
        for (std::size_t i = index; i < size - 1; i++)
            data[i] = data[i + 1];
        --size;
    }
};

static void line_case(const arr3 &simplex, glm::vec2 &dir)
{
    const glm::vec2 AB = simplex.data[0] - simplex.data[1], AO = -simplex.data[1];
    dir = triple_cross(AB, AO, AB);
}

static bool triangle_case(arr3 &simplex, glm::vec2 &dir)
{
    const glm::vec2 AB = simplex.data[1] - simplex.data[2], AC = simplex.data[0] - simplex.data[2],
                    AO = -simplex.data[2];
    const glm::vec2 ABperp = triple_cross(AC, AB, AB);
    if (glm::dot(ABperp, AO) >= 0.f)
    {
        simplex.erase(0);
        dir = ABperp;
        return false;
    }
    const glm::vec2 ACperp = triple_cross(AB, AC, AC);
    if (glm::dot(ACperp, AO) >= 0.f)
    {
        simplex.erase(1);
        dir = ACperp;
        return false;
    }
    return true;
}

gjk_result gjk(const shape2D &sh1, const shape2D &sh2)
{
    KIT_PERF_FUNCTION()
    KIT_ASSERT_WARN(!dynamic_cast<const circle *>(&sh1) || !dynamic_cast<const circle *>(&sh2),
                    "Using gjk algorithm to check if two circles are intersecting is overkill")

    gjk_result result{false, {}};
    arr3 simplex{result.simplex};

    glm::vec2 dir = sh2.centroid() - sh1.centroid();
    const glm::vec2 supp = sh1.support_point(dir) - sh2.support_point(-dir);
    simplex.push(supp);
    dir = -supp;

    for (;;)
    {
        const glm::vec2 A = sh1.support_point(dir) - sh2.support_point(-dir);
        if (glm::dot(A, dir) <= 0.f)
            return result;

        simplex.push(A);
        if (simplex.size == 2)
            line_case(simplex, dir);
        else if (triangle_case(simplex, dir))
        {
            result.intersect = true;
            return result;
        }
    }
}

epa_result epa(const shape2D &sh1, const shape2D &sh2, const std::array<glm::vec2, 3> &simplex)
{
    KIT_PERF_FUNCTION()

    std::vector<glm::vec2> hull;
    hull.reserve(10);
    hull.insert(hull.end(), simplex.begin(), simplex.end());

    float min_dist = FLT_MAX;
    epa_result result{false, glm::vec2(0.f)};
    for (;;)
    {
        std::size_t min_index;
        for (std::size_t i = 0; i < hull.size(); i++)
        {
            const std::size_t j = (i + 1) % hull.size();

            const glm::vec2 &p1 = hull[i], &p2 = hull[j];
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
        if (diff <= EPA_EPSILON)
            break;
        hull.insert(hull.begin() + (long)min_index, support);
        min_dist = FLT_MAX;
    }

    result.mtv *= min_dist;
    if (kit::approaches_zero(glm::length2(result.mtv)))
        return result;

    result.valid = true;
    return result;
}

glm::vec2 contact_point(const shape2D &sh1, const shape2D &sh2, const glm::vec2 &mtv)
{
    KIT_PERF_FUNCTION()
    const glm::vec2 sup1 = sh1.support_point(mtv), sup2 = sh2.support_point(-mtv);
    const float d1 = glm::length2(sh2.closest_direction_from(sup1 - mtv)),
                d2 = glm::length2(sh1.closest_direction_from(sup2 + mtv));
    return d1 < d2 ? sup1 : sup2 + mtv;
}

bool may_intersect(const shape2D &sh1, const shape2D &sh2)
{
    return intersects(sh1.bounding_box(), sh2.bounding_box());
}
bool intersects(const aabb2D &bb1, const aabb2D &bb2)
{
    const glm::vec2 df1 = bb2.min() - bb1.max(), df2 = bb1.min() - bb2.max();
    if (df1.x > 0.f || df1.y > 0.f)
        return false;
    if (df2.x > 0.f || df2.y > 0.f)
        return false;
    return true;
}

bool intersects(const circle &c1, const circle &c2)
{
    const float R = c1.radius + c2.radius;
    return glm::distance2(c1.centroid(), c2.centroid()) < R * R;
}
glm::vec2 mtv(const circle &c1, const circle &c2)
{
    const glm::vec2 dir = c1.centroid() - c2.centroid();
    return dir - (c1.radius + c2.radius) * glm::normalize(dir);
}
glm::vec2 contact_point(const circle &c1, const circle &c2)
{
    const glm::vec2 dir = glm::normalize(c1.centroid() - c2.centroid());
    return c1.centroid() - dir * c1.radius;
}
} // namespace geo