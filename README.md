# geometry

This project is a library offering useful geometry utilities, including implementations of various shapes (polygons, circles, AABBs) and pairwise functions for shape intersections and ray casts. It primarily serves as the geometry module for my 2D physics engine, [poly-physx](https://github.com/ismawno/poly-physx).

## Features

- Implementation of convex polygon shapes
- Implementation of circle shapes
- AABB (Axis-Aligned Bounding Box) implementation for broad-phase collision detection
- Ray casting functionalities

## Dependencies

- [glm](https://github.com/g-truc/glm)
- [cpp-kit](https://github.com/ismawno/cpp-kit)
- [yaml-cpp](https://github.com/ismawno/yaml-cpp) (optional)
- [spdlog](https://github.com/gabime/spdlog) (optional)

## Building and Usage

This project is intended to be used as a git submodule within another project (parent repo). A premake file is provided for building and linking geometry.

While these build instructions are minimal, this project is primarily for personal use. Although it has been built and tested on multiple machines (MacOS and Windows), it is not necessarily fully cross-platform or easy to build.

## License

geometry is licensed under the MIT License. See LICENSE for more details.
