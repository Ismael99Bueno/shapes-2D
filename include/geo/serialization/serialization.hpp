#pragma once
#ifdef KIT_USE_YAML_CPP

#include "geo/shapes2D/circle.hpp"
#include "geo/shapes2D/vertices2D.hpp"
#include "kit/serialization/yaml/codec.hpp"
#include "kit/serialization/yaml/glm.hpp"
#include "kit/serialization/yaml/transform.hpp"

namespace geo
{
template <std::size_t Capacity> class polygon;
}

template <> struct kit::yaml::codec<geo::aabb2D>
{
    static YAML::Node encode(const geo::aabb2D &aabb)
    {
        YAML::Node node;
        node["Min"] = aabb.min;
        node["Max"] = aabb.max;
        node["Min"].SetStyle(YAML::EmitterStyle::Flow);
        node["Max"].SetStyle(YAML::EmitterStyle::Flow);
        return node;
    }
    static bool decode(const YAML::Node &node, geo::aabb2D &aabb)
    {
        if (!node.IsMap() || node.size() != 2)
            return false;

        aabb = {node["Min"].as<glm::vec2>(), node["Max"].as<glm::vec2>()};
        return true;
    }
};

template <> struct kit::yaml::codec<geo::circle>
{
    static YAML::Node encode(const geo::circle &circ)
    {
        YAML::Node node;
        node["Transform"] = circ.ltransform().get();
        node["Radius"] = circ.radius();

        return node;
    }
    static bool decode(const YAML::Node &node, geo::circle &circ)
    {
        if (!node.IsMap() || node.size() != 2)
            return false;
        circ = {geo::transform2D{node["Transform"].as<kit::transform2D<float>>()}, node["Radius"].as<float>()};
        return true;
    }
};

template <std::size_t Capacity> struct kit::yaml::codec<geo::polygon<Capacity>>
{
    static YAML::Node encode(const geo::polygon<Capacity> &poly)
    {
        YAML::Node node;
        node["Transform"] = poly.ltransform().get();

        for (std::size_t i = 0; i < poly.vertices.size(); i++)
        {
            node["Vertices"].push_back(poly.vertices.locals[i]);
            node["Vertices"][i].SetStyle(YAML::EmitterStyle::Flow);
        }
        return node;
    }
    static bool decode(const YAML::Node &node, geo::polygon<Capacity> &poly)
    {
        if (!node.IsMap() || node.size() != 2)
            return false;
        YAML::Node node_v = node["Vertices"];

        kit::dynarray<glm::vec2, Capacity> vertices{node_v.size()};
        for (std::size_t i = 0; i < node_v.size(); i++)
            vertices[i] = node_v[i].as<glm::vec2>();

        const kit::transform2D<float> transform = node["Transform"].as<kit::transform2D<float>>();
        poly = {geo::transform2D{transform}, vertices};
        return true;
    }
};
#endif