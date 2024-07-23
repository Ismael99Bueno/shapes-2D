#include "geo/internal/pch.hpp"
#include "geo/shapes/nsphere.hpp"
#include "geo/serialization/serialization.hpp"

namespace geo
{
circle::circle(const float radius)
{
    KIT_ASSERT_WARN(radius >= 0.f, "Creating circle with negative radius: {0}", radius);
    m_radius = radius;
    m_convex = true;
    update_area_and_inertia();
    update();
}
circle::circle(const transform2D &ltransform, const float radius) : shape2D(ltransform)
{
    KIT_ASSERT_WARN(radius >= 0.f, "Creating circle with negative radius: {0}", radius);
    m_radius = radius;
    m_convex = true;
    update_area_and_inertia();
    update();
}

void circle::radius(const float radius)
{
    KIT_ASSERT_WARN(radius >= 0.f, "Setting circle radius to negative value: {0}", radius);
    m_radius = radius;
    update_area_and_inertia();
}

glm::vec2 circle::support_point(const glm::vec2 &direction) const
{
    return m_gcentroid + glm::normalize(direction) * m_radius;
}

bool circle::contains_point(const glm::vec2 &p) const
{
    return glm::length2(p - m_gcentroid) < m_radius * m_radius;
}

void circle::update_area_and_inertia()
{
    const float r2 = m_radius * m_radius;
    m_area = glm::pi<float>() * r2;
    m_inertia = 0.5f * r2;
}

aabb2D circle::create_bounding_box() const
{
    return aabb2D(m_gcentroid - m_radius, m_gcentroid + m_radius);
}

glm::vec2 circle::closest_direction_from(const glm::vec2 &p) const
{
    const glm::vec2 dir = m_gcentroid - p;
    return dir - glm::normalize(dir) * m_radius;
}

#ifdef KIT_USE_YAML_CPP
YAML::Node circle::encode() const
{
    return kit::yaml::codec<circle>::encode(*this);
}
bool circle::decode(const YAML::Node &node)
{
    return kit::yaml::codec<circle>::decode(node, *this);
}
#endif

} // namespace geo
