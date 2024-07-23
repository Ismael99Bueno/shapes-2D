#pragma once

#include "geo/shapes2D/shape.hpp"

namespace geo
{
class circle final : public shape2D
{
  public:
    circle(float radius = 1.f);
    circle(const transform2D &ltransform, float radius = 1.f);

    using shape2D::radius;
    void radius(float radius);

    glm::vec2 support_point(const glm::vec2 &direction) const override;
    bool contains_point(const glm::vec2 &p) const override;

    glm::vec2 closest_direction_from(const glm::vec2 &p) const override;
    aabb2D create_bounding_box() const override;

#ifdef KIT_USE_YAML_CPP
    YAML::Node encode() const override;
    bool decode(const YAML::Node &node) override;
#endif

  private:
    void update_area_and_inertia();
};
} // namespace geo
