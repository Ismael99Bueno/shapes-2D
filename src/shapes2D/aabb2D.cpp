#include "geo/internal/pch.hpp"
#include "geo/shapes2D/aabb2D.hpp"

namespace geo
{
aabb2D::aabb2D(const glm::vec2 &point) : aabb2D(point, point)
{
}
aabb2D::aabb2D(const glm::vec2 &min, const glm::vec2 &max) : min(min), max(max)
{
    KIT_ASSERT_ERROR(min.x <= max.x && min.y <= max.y,
                     "Vector min: {0}, {1} is greater than vector max {0}, {1} in aabb", min.x, min.y, max.x, max.y);
}

glm::vec2 aabb2D::dimension() const
{
    return max - min;
}
bool aabb2D::contains(const aabb2D &aabb) const
{
    return min.x <= aabb.min.x && min.y <= aabb.min.y && max.x >= aabb.max.x && max.y >= aabb.max.y;
}
float aabb2D::area() const
{
    glm::vec2 dim = dimension();
    return dim.x * dim.y;
}

aabb2D &aabb2D::operator+=(const aabb2D &bb)
{
    *this = *this + bb;
    return *this;
}
aabb2D &aabb2D::operator-=(const aabb2D &bb)
{
    *this = *this - bb;
    return *this;
}

aabb2D operator+(const aabb2D &bb1, const aabb2D &bb2)
{
    return aabb2D({std::min(bb1.min.x, bb2.min.x), std::min(bb1.min.y, bb2.min.y)},
                  {std::max(bb1.max.x, bb2.max.x), std::max(bb1.max.y, bb2.max.y)});
}
aabb2D operator-(const aabb2D &bb1, const aabb2D &bb2)
{
    return aabb2D({std::max(bb1.min.x, bb2.min.x), std::max(bb1.min.y, bb2.min.y)},
                  {std::min(bb1.max.x, bb2.max.x), std::min(bb1.max.y, bb2.max.y)});
}
} // namespace geo