project "geometry"
language "C++"
cppdialect "c++20"

filter "system:macosx or linux"
   buildoptions {
      "-Wall",
      "-Wextra",
      "-Wpedantic",
      "-Wconversion",
      "-Wno-unused-parameter",
      "-Wno-sign-conversion",
      "-Wno-gnu-anonymous-struct",
      "-Wno-nested-anon-types",
      "-Wno-string-conversion"
   }
filter {}

pchheader "geo/internal/pch.hpp"
pchsource "src/internal/pch.cpp"

staticruntime "off"
kind "StaticLib"

targetdir("bin/" .. outputdir)
objdir("build/" .. outputdir)

files {
   "src/**.cpp",
   "include/**.hpp"
}

includedirs {
   "include",
   "%{wks.location}/cpp-kit/include",
   "%{wks.location}/vendor/yaml-cpp/include",
   "%{wks.location}/vendor/glm",
   "%{wks.location}/vendor/spdlog/include"
}
