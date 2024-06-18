#include "geo/internal/pch.hpp"
#include "geo/algorithm/ray2D.hpp"

namespace geo
{
ray2D::ray2D(const glm::vec2 &origin, const glm::vec2 &direction)
    : m_origin(origin), m_dir(glm::normalize(direction)), m_normal(-m_dir.y, m_dir.x), m_infinte(true)
{
}

ray2D::ray2D(const glm::vec2 &origin, const glm::vec2 &direction, float length)
    : m_origin(origin), m_dir(glm::normalize(direction)), m_normal(-m_dir.y, m_dir.x), m_length(length),
      m_infinte(false)
{
}

ray2D ray2D::from_endpoints(const glm::vec2 &start, const glm::vec2 &end)
{
    return ray2D(start, end - start, glm::length(end - start));
}

const glm::vec2 &ray2D::origin() const
{
    return m_origin;
}
const glm::vec2 &ray2D::direction() const
{
    return m_dir;
}
const glm::vec2 &ray2D::normal() const
{
    return m_normal;
}

float ray2D::length() const
{
    return m_length;
}
bool ray2D::infinite() const
{
    return m_infinte;
}

void ray2D::resize(float length)
{
    if (length == FLT_MAX)
    {
        m_infinte = true;
        return;
    }
    m_length = length;
    m_infinte = false;
}

} // namespace geo