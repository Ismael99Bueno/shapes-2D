#pragma once

#include "geo/shapes2D/shape.hpp"
#include "geo/serialization/serialization.hpp"
#include "geo/shapes2D/vertices.hpp"
#include "kit/utility/utils.hpp"
#include <vector>
#include <array>
#include <utility>

namespace geo
{
template <std::size_t Capacity> class polygon final : public shape2D
{
  public:
    struct vertex_container
    {
        vertices2D<Capacity> locals;
        vertices2D<Capacity> globals;
        vertices2D<Capacity> edges;
        vertices2D<Capacity> normals;
        vertices2D<Capacity> model;
        std::size_t size() const
        {
            return locals.size();
        }
    };

    template <std::input_iterator It>
    polygon(It it1, It it2)
        : vertices{.locals{it1, it2}, .globals{it1, it2}, .edges{it1, it2}, .normals{it1, it2}, .model{it1, it2}}
    {
        m_ltransform.position = initialize_properties_and_vertices();
        update();
    }

    template <std::size_t Size = 4>
        requires(Size >= 3 && Size <= Capacity)
    polygon(const kit::dynarray<glm::vec2, Size> &verts = square(1.f))
        : vertices{.locals{verts},
                   .globals{verts.size()},
                   .edges{verts.size()},
                   .normals{verts.size()},
                   .model{verts.size()}}
    {
        m_ltransform.position = initialize_properties_and_vertices();
        update();
    }
    polygon(std::initializer_list<glm::vec2> verts)
        : vertices{.locals{verts},
                   .globals{verts.size()},
                   .edges{verts.size()},
                   .normals{verts.size()},
                   .model{verts.size()}}
    {
        m_ltransform.position = initialize_properties_and_vertices();
        update();
    }

    template <std::input_iterator It>
    polygon(const transform2D &ltransform, It it1, It it2)
        : shape2D(ltransform),
          vertices{.locals{it1, it2}, .globals{it1, it2}, .edges{it1, it2}, .normals{it1, it2}, .model{it1, it2}}
    {
        initialize_properties_and_vertices();
        update();
    }

    template <std::size_t Size>
        requires(Size >= 3 && Size <= Capacity)
    polygon(const transform2D &ltransform, const kit::dynarray<glm::vec2, Size> &verts = square(1.f))
        : shape2D(ltransform), vertices{.locals{verts},
                                        .globals{verts.size()},
                                        .edges{verts.size()},
                                        .normals{verts.size()},
                                        .model{verts.size()}}
    {
        initialize_properties_and_vertices();
        update();
    }
    polygon(const transform2D &ltransform, std::initializer_list<glm::vec2> verts)
        : shape2D(ltransform), vertices{.locals{verts},
                                        .globals{verts.size()},
                                        .edges{verts.size()},
                                        .normals{verts.size()},
                                        .model{verts.size()}}
    {
        initialize_properties_and_vertices();
        update();
    }

    vertex_container vertices;

    glm::vec2 support_point(const glm::vec2 &direction) const override
    {
        std::size_t support = 0;
        float max_dot = glm::dot(direction, vertices.globals[support] - m_gcentroid);
        for (std::size_t i = 1; i < vertices.size(); i++)
        {
            const float dot = glm::dot(direction, vertices.globals[i] - m_gcentroid);
            if (dot > max_dot)
            {
                max_dot = dot;
                support = i;
            }
        }
        return vertices.globals[support];
    }

    bool contains_point(const glm::vec2 &p) const override
    {
        KIT_ASSERT_WARN(m_convex,
                        "Checking if a point is contained in a non convex polygon yields undefined behaviour.")
        for (std::size_t i = 0; i < vertices.size(); i++)
        {
            const glm::vec2 &normal = vertices.normals[i];
            const glm::vec2 side = p - vertices.globals[i];
            if (glm::dot(normal, side) > 0.f)
                return false;
        }
        return true;
    }

    glm::vec2 closest_direction_from(const glm::vec2 &p) const override
    {
        float min_dist = FLT_MAX;
        glm::vec2 closest(0.f);
        for (std::size_t i = 0; i < vertices.size(); i++)
        {
            const float interp = glm::clamp(
                glm::dot(p - vertices.globals[i], vertices.edges[i]) / glm::length2(vertices.edges[i]), 0.f, 1.f);
            const glm::vec2 towards = vertices.globals[i] + interp * vertices.edges[i] - p;
            const float dist = glm::length2(towards);
            if (min_dist > dist)
            {
                min_dist = dist;
                closest = towards;
            }
        }
        return closest;
    }

    aabb2D create_bounding_box() const override
    {
        glm::vec2 min = vertices.globals[0];
        glm::vec2 max = vertices.globals[0];
        for (std::size_t i = 1; i < vertices.size(); i++)
        {
            const glm::vec2 &v = vertices.globals[i];
            if (min.x > v.x)
                min.x = v.x;
            if (min.y > v.y)
                min.y = v.y;
            if (max.x < v.x)
                max.x = v.x;
            if (max.y < v.y)
                max.y = v.y;
        }
        return aabb2D(min, max);
    }

    static kit::dynarray<glm::vec2, 4> square(const float size)
    {
        const float hsize = 0.5f * size;
        return kit::dynarray<glm::vec2, 4>{
            {glm::vec2(-hsize, -hsize), glm::vec2(hsize, -hsize), glm::vec2(hsize, hsize), glm::vec2(-hsize, hsize)}};
    }
    static kit::dynarray<glm::vec2, 4> rect(const float width, const float height)
    {
        const float hw = 0.5f * width;
        const float hh = 0.5f * height;
        return kit::dynarray<glm::vec2, 4>{
            {glm::vec2(-hw, -hh), glm::vec2(hw, -hh), glm::vec2(hw, hh), glm::vec2(-hw, hh)}};
    }

    template <std::size_t MaxEdges = Capacity>
    static kit::dynarray<glm::vec2, MaxEdges> ngon(const float radius, const std::uint32_t edges)
    {
        KIT_ASSERT_ERROR(edges >= 3, "Cannot make polygon with less than 3 edges - edges: {0}", edges)
        KIT_ASSERT_ERROR(edges <= MaxEdges, "Cannot make polygon with more than {0} edges - edges: {1}", MaxEdges,
                         edges)
        kit::dynarray<glm::vec2, MaxEdges> vertices{edges};
        const float dangle = 2.f * glm::pi<float>() / edges;
        for (std::size_t i = 0; i < edges; i++)
        {
            const float rotation = i * dangle;
            vertices[i] = {radius * sinf(rotation), radius * cosf(rotation)};
        }
        return vertices;
    }

#ifdef KIT_USE_YAML_CPP
    YAML::Node encode() const override
    {
        return kit::yaml::codec<polygon>::encode(*this);
    }
    bool decode(const YAML::Node &node) override
    {
        return kit::yaml::codec<polygon>::decode(node, *this);
    }
#endif

  private:
    void on_shape_transform_update(const glm::mat3 &ltransform, const glm::mat3 &gtransform) override
    {
        shape2D::on_shape_transform_update(ltransform, gtransform);
        for (std::size_t i = 0; i < vertices.size(); i++)
            vertices.locals(i) = ltransform * glm::vec3(vertices.model[i], 1.f);
        if (m_ltransform.parent)
            for (std::size_t i = 0; i < vertices.size(); i++)
                vertices.globals(i) = gtransform * glm::vec3(vertices.model[i], 1.f);
        else
            for (std::size_t i = 0; i < vertices.size(); i++)
                vertices.globals(i) = vertices.locals[i];

        for (std::size_t i = 0; i < vertices.size(); i++)
        {
            vertices.edges(i) = vertices.globals[i + 1] - vertices.globals[i];
            vertices.normals(i) = glm::normalize(glm::vec2(vertices.edges[i].y, -vertices.edges[i].x));
        }
    }

    void sort_local_vertices()
    {
        const glm::vec2 center = center_of_local_vertices();
        const glm::vec2 reference = vertices.locals[0] - center;

        const auto cmp = [&center, &reference](const glm::vec2 &v1, const glm::vec2 &v2) {
            const glm::vec2 dir1 = v1 - center, dir2 = v2 - center;

            const float det2 = kit::cross2D(reference, dir2);
            if (kit::approaches_zero(det2) && glm::dot(reference, dir2) >= 0.f)
                return false;
            const float det1 = kit::cross2D(reference, dir1);
            if (kit::approaches_zero(det1) && glm::dot(reference, dir1) >= 0.f)
                return true;

            if (det1 * det2 >= 0.f)
                return kit::cross2D(dir1, dir2) > 0.f;
            return det1 > 0.f;
        };
        std::sort(vertices.locals.mbegin(), vertices.locals.mend(), cmp);
    }

    glm::vec2 initialize_properties_and_vertices()
    {
        sort_local_vertices();
        const glm::vec2 current_lcentroid = compute_center_of_mass();

        for (std::size_t i = 0; i < vertices.size(); i++)
            vertices.model(i) = vertices.locals[i] - current_lcentroid;

        m_area = compute_area();
        m_inertia = compute_inertia();
        m_convex = compute_convexity();
        m_radius = compute_radius();
        return current_lcentroid;
    }

    glm::vec2 center_of_local_vertices()
    {
        glm::vec2 center(0.f);
        for (const glm::vec2 &v : vertices.locals)
            center += v;
        return center / (float)vertices.locals.size();
    }

    glm::vec2 compute_center_of_mass()
    {
        const glm::vec2 &p1 = vertices.locals[0]; // No model available yet
        glm::vec2 num(0.f), den(0.f);
        for (std::size_t i = 1; i < vertices.locals.size() - 1; i++)
        {
            const glm::vec2 e1 = vertices.locals[i] - p1;
            const glm::vec2 e2 = vertices.locals[i + 1] - p1;

            const float crs = std::abs(kit::cross2D(e1, e2));
            num += (e1 + e2) * crs;
            den += crs;
        }
        return p1 + num / (3.f * den);
    }

    float compute_area()
    {
        float area = 0.f;
        const glm::vec2 &p1 = vertices.model[0];

        for (std::size_t i = 1; i < vertices.size() - 1; i++)
        {
            const glm::vec2 e1 = vertices.model[i] - p1;
            const glm::vec2 e2 = vertices.model[i + 1] - p1;
            area += std::abs(kit::cross2D(e1, e2));
        }
        return area * 0.5f;
    }

    float compute_inertia()
    {
        const glm::vec2 &p1 = vertices.model[0];
        float inertia = 0.f;

        for (std::size_t i = 1; i < vertices.size() - 1; i++)
        {
            const glm::vec2 &p2 = vertices.model[i];
            const glm::vec2 &p3 = vertices.model[i + 1];
            const glm::vec2 e1 = p1 - p2;
            const glm::vec2 e2 = p3 - p2;

            const float w = glm::length(e1);

            const float w1 = std::abs(glm::dot(e1, e2) / w);
            const float w2 = std::abs(w - w1);

            const float h = std::abs(kit::cross2D(e2, e1)) / w;
            const glm::vec2 p4 = p2 + e1 * w1 / w;

            const float i1 = w1 * h * (w1 * w1 / 3.f + h * h) / 4.f;
            const float i2 = w2 * h * (w2 * w2 / 3.f + h * h) / 4.f;

            const float m1 = 0.5f * w1 * h;
            const float m2 = 0.5f * w2 * h;

            const glm::vec2 cm1 = (p2 + p3 + p4) / 3.f;
            const glm::vec2 cm2 = (p1 + p3 + p4) / 3.f;

            const float icm1 = i1 + m1 * (glm::length2(cm1) - glm::distance2(cm1, p3));
            const float icm2 = i2 + m2 * (glm::length2(cm2) - glm::distance2(cm2, p3));

            const glm::vec2 p13 = p1 - p3;
            const glm::vec2 p43 = p4 - p3;
            const glm::vec2 p23 = -e2;

            if (kit::cross2D(p13, p43) < 0.f)
                inertia += icm1;
            else
                inertia -= icm1;
            if (kit::cross2D(p43, p23) < 0.f)
                inertia += icm2;
            else
                inertia -= icm2;
        }
        return std::abs(inertia) / m_area;
    }

    bool compute_convexity() const
    {
        for (std::size_t i = 0; i < vertices.size(); i++) // No edges available yet
            if (kit::cross2D(vertices.model[i + 1] - vertices.model[i], vertices.model[i + 2] - vertices.model[i + 1]) <
                0.f)
                return false;

        return true;
    }
    float compute_radius() const
    {
        float max_dist = 0.f;
        for (std::size_t i = 0; i < vertices.size(); i++)
        {
            const float dist = glm::length2(vertices.model[i]);
            if (max_dist < dist)
                max_dist = dist;
        }
        return glm::sqrt(max_dist);
    }
};

template <typename T>
concept Polygon = Shape2D<T> && std::same_as<typename T::vertex_container, typename polygon<3>::vertex_container>;

} // namespace geo
