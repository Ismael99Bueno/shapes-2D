#pragma once

#include <glm/vec2.hpp>
#include <glm/mat2x2.hpp>
#include "geo/shapes2D/aabb2D.hpp"
#include "geo/algorithm/transform2D.hpp"

#include "kit/serialization/yaml/serializer.hpp"

namespace geo
{
class shape2D : public kit::yaml::serializable, public kit::yaml::deserializable
{
  public:
    shape2D(const kit::transform2D<float> &ltransform);
    shape2D() = default;
    virtual ~shape2D() = default;

    virtual glm::vec2 support_point(const glm::vec2 &direction) const = 0;
    virtual bool contains_point(const glm::vec2 &p) const = 0;
    bool contains_origin() const;

    void begin_update();
    void end_update();

    const transform2D *parent() const;
    void parent(const transform2D *parent);

    const aabb2D &bounding_box() const;
    virtual void bound() = 0;

    virtual glm::vec2 closest_direction_from(const glm::vec2 &p) const = 0;

    const transform2D &ltransform() const;
    void ltransform(const kit::transform2D<float> &ltransform);

    const glm::vec2 &lcentroid() const;
    const glm::vec2 &gcentroid() const;

    const glm::vec2 &lposition() const;
    const glm::vec2 &lscale() const;
    float lrotation() const;
    const glm::vec2 &origin() const;

    void ltranslate(const glm::vec2 &dpos);
    void gtranslate(const glm::vec2 &dpos);
    void lrotate(float drotation);

    void lcentroid(const glm::vec2 &lcentroid);
    void gcentroid(const glm::vec2 &gcentroid);

    void lposition(const glm::vec2 &lposition);
    void lscale(const glm::vec2 &lscale);
    void lrotation(float langle);
    void origin(const glm::vec2 &origin);

    float area() const;
    float inertia() const;
    bool convex() const;

    bool updating() const;
    void update();

  protected:
    transform2D m_ltransform;
    glm::vec2 m_lcentroid;
    glm::vec2 m_gcentroid;
    aabb2D m_aabb;

    float m_area = 0.f;
    float m_inertia = 0.f;
    bool m_convex = true;

    virtual void on_shape_transform_update(const glm::mat3 &ltransform, const glm::mat3 &gtransform);

  private:
    bool m_pushing_update = false;
};

} // namespace geo
