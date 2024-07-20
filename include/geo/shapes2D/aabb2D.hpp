#pragma once

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#define GLM_FORCE_INTRINSICS
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <vector>

namespace geo
{
struct aabb2D
{
    aabb2D(const glm::vec2 &point = glm::vec2(0.f));
    aabb2D(const glm::vec2 &min, const glm::vec2 &max);

    glm::vec2 dimension() const;
    bool contains(const aabb2D &aabb) const;
    float area() const;

    void enlarge(const glm::vec2 &enlarge_vector);
    void enlarge(float buffer);

    glm::vec2 min;
    glm::vec2 max;

    aabb2D &operator+=(const aabb2D &bb);
    aabb2D &operator-=(const aabb2D &bb);
};

aabb2D operator+(const aabb2D &bb1, const aabb2D &bb2);
aabb2D operator-(const aabb2D &bb1, const aabb2D &bb2);

} // namespace geo
