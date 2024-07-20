#pragma once

#include "kit/container/dynarray.hpp"
#include "kit/utility/type_constraints.hpp"
#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
// #define GLM_FORCE_AVX2
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>

namespace geo
{
template <std::size_t N> class polygon;
template <std::size_t N>
    requires(N >= 3)
class vertices2D
{
  public:
    template <class... VectorArgs>
        requires kit::NoCopyCtorOverride<vertices2D, VectorArgs...>
    vertices2D(VectorArgs &&...args) : m_vertices(std::forward<VectorArgs>(args)...)
    {
    }

    const glm::vec2 &operator[](const std::size_t index) const
    {
        return m_vertices[index % m_vertices.size()];
    }

    auto begin() const
    {
        return m_vertices.begin();
    }
    auto end() const
    {
        return m_vertices.end();
    }

    std::size_t size() const
    {
        return m_vertices.size();
    }

    operator const kit::dynarray<glm::vec2, N> &() const
    {
        return m_vertices;
    }

  private:
    kit::dynarray<glm::vec2, N> m_vertices;

    vertices2D(const vertices2D &) = default;
    vertices2D &operator=(const vertices2D &) = default;
    vertices2D(vertices2D &&) = default;
    vertices2D &operator=(vertices2D &&) = default;

    glm::vec2 &operator()(const std::size_t index)
    {
        return m_vertices[index % m_vertices.size()];
    }

    auto mbegin()
    {
        return m_vertices.begin();
    }
    auto mend()
    {
        return m_vertices.end();
    }

    friend class polygon<N>;
};
} // namespace geo