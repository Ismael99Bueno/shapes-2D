#pragma once

#include <glm/vec2.hpp>
#include <vector>

namespace geo
{
struct aabb2D
{
    aabb2D(const glm::vec2 &point = glm::vec2(0.f));
    aabb2D(const glm::vec2 &min, const glm::vec2 &max);

    glm::vec2 dimension() const;
    bool contains(const aabb2D &aabb) const;

    glm::vec2 min;
    glm::vec2 max;

    aabb2D &operator+=(const aabb2D &bb);
    aabb2D &operator-=(const aabb2D &bb);
};

aabb2D operator+(const aabb2D &bb1, const aabb2D &bb2);
aabb2D operator-(const aabb2D &bb1, const aabb2D &bb2);

} // namespace geo
