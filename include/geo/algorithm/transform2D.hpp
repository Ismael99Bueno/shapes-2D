#pragma once

#include "kit/utility/transform.hpp"

namespace geo
{
class transform2D
{
  public:
    transform2D() = default;
    transform2D(const kit::transform2D<float> &transform);

    const glm::vec2 &position() const;
    const glm::vec2 &scale() const;
    const glm::vec2 &origin() const;
    float rotation() const;
    const transform2D *parent() const;

    const glm::mat3 &ltransform() const;
    glm::mat3 gtransform() const;

    const glm::mat3 &inv_ltransform() const;
    glm::mat3 inv_gtransform() const;

    void position(const glm::vec2 &position);
    void scale(const glm::vec2 &scale);
    void origin(const glm::vec2 &origin);
    void rotation(float rotation);
    void parent(const transform2D *parent);

    void translate(const glm::vec2 &dpos);
    void rotate(float drotation);

    operator const kit::transform2D<float> &();

  private:
    kit::transform2D<float> m_transform;
    const transform2D *m_parent = nullptr;

    mutable glm::mat3 m_tmat;
    mutable glm::mat3 m_inv_tmat;

    mutable bool m_cache = false;
    mutable bool m_inv_cache = false;
};
} // namespace geo