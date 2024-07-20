#pragma once
#include <vector>
#include <utility>
#include <array>
#include <limits>

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
// #define GLM_FORCE_AVX2
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <optional>
#ifdef KIT_USE_YAML_CPP
#include <yaml-cpp/yaml.h>
#endif
#include "kit/debug/log.hpp"
#include "kit/profiling/perf.hpp"
