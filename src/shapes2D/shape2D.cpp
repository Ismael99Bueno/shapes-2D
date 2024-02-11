#include "geo/internal/pch.hpp"
#include "geo/shapes2D/shape2D.hpp"

namespace geo
{
shape2D::shape2D(const kit::transform2D<float> &ltransform) : m_ltransform(ltransform)
{
}

const kit::transform2D<float> &shape2D::ltransform() const
{
    return m_ltransform;
}

const glm::vec2 &shape2D::lcentroid() const
{
    return m_lcentroid;
}
const glm::vec2 &shape2D::gcentroid() const
{
    return m_gcentroid;
}

const glm::vec2 &shape2D::lposition() const
{
    return m_ltransform.position;
}
const glm::vec2 &shape2D::lscale() const
{
    return m_ltransform.scale;
}
float shape2D::lrotation() const
{
    return m_ltransform.rotation;
}
const glm::vec2 &shape2D::origin() const
{
    return m_ltransform.origin;
}

float shape2D::area() const
{
    return m_area;
}
float shape2D::inertia() const
{
    return m_inertia;
}
bool shape2D::convex() const
{
    return m_convex;
}

void shape2D::lcentroid(const glm::vec2 &lcentroid)
{
    ltranslate(lcentroid - m_lcentroid);
}

const kit::transform2D<float> *shape2D::parent() const
{
    return m_ltransform.parent;
}
void shape2D::parent(const kit::transform2D<float> *parent)
{
    if (m_ltransform.parent == parent)
        return;
    m_ltransform.parent = parent;
    update();
}

void shape2D::update()
{
    if (m_pushing_update)
        return;
    const glm::mat3 ltransform = m_ltransform.center_scale_rotate_translate3(true);
    if (m_ltransform.parent)
    {
        const glm::mat3 gtransform = m_ltransform.parent->center_scale_rotate_translate3() * ltransform;
        on_shape_transform_update(ltransform, gtransform);
    }
    else
        on_shape_transform_update(ltransform, ltransform);

    bound();
}

void shape2D::on_shape_transform_update(const glm::mat3 &ltransform, const glm::mat3 &gtransform)
{
    m_lcentroid = ltransform[2];
    m_gcentroid = gtransform[2];
}

void shape2D::ltranslate(const glm::vec2 &dpos)
{
    m_ltransform.position += dpos;
    update();
}
void shape2D::lrotate(const float drotation)
{
    m_ltransform.rotation += drotation;
    update();
}

void shape2D::lposition(const glm::vec2 &lposition)
{
    m_ltransform.position = lposition;
    update();
}
void shape2D::lscale(const glm::vec2 &lscale)
{
    m_ltransform.scale = lscale;
    update();
}

void shape2D::lrotation(const float lrotation)
{
    m_ltransform.rotation = lrotation;
    update();
}

void shape2D::origin(const glm::vec2 &origin)
{
    m_ltransform.origin = origin;
    update();
}

bool shape2D::contains_origin() const
{
    return contains_point(glm::vec2(0.f));
}
const aabb2D &shape2D::bounding_box() const
{
    return m_aabb;
}

bool shape2D::updating() const
{
    return m_pushing_update;
}

void shape2D::begin_update()
{
    KIT_ASSERT_ERROR(!m_pushing_update, "Cannot begin update while already updating");
    m_pushing_update = true;
}
void shape2D::end_update()
{
    m_pushing_update = false;
    update();
}
} // namespace geo