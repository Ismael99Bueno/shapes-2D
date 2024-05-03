#include <glm/vec2.hpp>
#include <vector>
#include <utility>
#include <array>
#include <limits>

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#include <glm/geometric.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/mat2x2.hpp>

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <optional>
#ifdef KIT_USE_YAML_CPP
#include <yaml-cpp/yaml.h>
#endif
#include "kit/debug/log.hpp"
#include "kit/profiling/perf.hpp"
