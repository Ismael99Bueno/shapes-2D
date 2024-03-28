#include "geo/internal/pch.hpp"
#include "geo/algorithm/transform2D.hpp"

namespace geo
{
transform2D::transform2D(const kit::transform2D<float> &transform) : m_transform(transform)
{
}

const glm::vec2 &transform2D::position() const
{
    return m_transform.position;
}
const glm::vec2 &transform2D::scale() const
{
    return m_transform.scale;
}
const glm::vec2 &transform2D::origin() const
{
    return m_transform.origin;
}
float transform2D::rotation() const
{
    return m_transform.rotation;
}

const transform2D *transform2D::parent() const
{
    return m_parent;
}

const glm::mat3 &transform2D::ltransform() const
{
    if (!m_cache)
    {
        m_tmat = m_transform.center_scale_rotate_translate3();
        m_cache = true;
    }
    return m_tmat;
}
glm::mat3 transform2D::gtransform() const
{
    if (m_parent)
        return m_parent->gtransform() * ltransform();
    return ltransform();
}

const glm::mat3 &transform2D::inv_ltransform() const
{
    if (!m_inv_cache)
    {
        m_inv_tmat = m_transform.inverse_center_scale_rotate_translate3();
        m_inv_cache = true;
    }
    return m_inv_tmat;
}
glm::mat3 transform2D::inv_gtransform() const
{
    if (m_parent)
        return inv_ltransform() * m_parent->inv_gtransform();
    return inv_ltransform();
}

void transform2D::position(const glm::vec2 &position)
{
    m_transform.position = position;
    m_cache = false;
    m_inv_cache = false;
}
void transform2D::scale(const glm::vec2 &scale)
{
    m_transform.scale = scale;
    m_cache = false;
    m_inv_cache = false;
}
void transform2D::origin(const glm::vec2 &origin)
{
    m_transform.origin = origin;
    m_cache = false;
    m_inv_cache = false;
}
void transform2D::rotation(float rotation)
{
    m_transform.rotation = rotation;
    m_cache = false;
    m_inv_cache = false;
}
void transform2D::parent(const transform2D *parent)
{
    m_parent = parent;
}

void transform2D::translate(const glm::vec2 &dpos)
{
    m_transform.position += dpos;
    m_cache = false;
    m_inv_cache = false;
}
void transform2D::rotate(float drotation)
{
    m_transform.rotation += drotation;
    m_cache = false;
    m_inv_cache = false;
}

transform2D::operator const kit::transform2D<float> &()
{
    return m_transform;
}

} // namespace geo