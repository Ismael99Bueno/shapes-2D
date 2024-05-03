#include "geo/internal/pch.hpp"
#include "geo/algorithm/ray2D.hpp"

namespace geo
{
ray2D::ray2D(const glm::vec2 &origin, const glm::vec2 &direction)
    : m_origin(origin), m_dir(glm::normalize(direction)), m_infinte(true)
{
}

ray2D::ray2D(const glm::vec2 &origin, const glm::vec2 &direction, float length)
    : m_origin(origin), m_dir(glm::normalize(direction)), m_length(length), m_infinte(false)
{
}

ray2D ray2D::from_endpoints(const glm::vec2 &start, const glm::vec2 &end)
{
    return ray2D(start, end - start, glm::length(end - start));
}

ray2D::hit::operator bool() const
{
    return hit;
}

const glm::vec2 &ray2D::origin() const
{
    return m_origin;
}
const glm::vec2 &ray2D::direction() const
{
    return m_dir;
}
float ray2D::length() const
{
    return m_length;
}
bool ray2D::infinite() const
{
    return m_infinte;
}
} // namespace geo