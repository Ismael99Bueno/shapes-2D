#pragma once

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#include <glm/vec2.hpp>

namespace geo
{
class ray2D
{
  public:
    template <typename T> struct hit
    {
        glm::vec2 point;
        glm::vec2 normal;
        float distance;
        bool is_hit = false;
        T *object = nullptr;

        operator bool() const
        {
            return is_hit;
        }
    };

    ray2D() = default;
    ray2D(const glm::vec2 &origin, const glm::vec2 &direction);
    ray2D(const glm::vec2 &origin, const glm::vec2 &direction, float length);

    static ray2D from_endpoints(const glm::vec2 &start, const glm::vec2 &end);

    const glm::vec2 &origin() const;
    const glm::vec2 &direction() const;
    const glm::vec2 &normal() const;

    float length() const;
    bool infinite() const;

    void resize(float length);

  private:
    glm::vec2 m_origin;
    glm::vec2 m_dir;
    glm::vec2 m_normal;
    float m_length;
    bool m_infinte;
};
} // namespace geo