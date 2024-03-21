-- project
set_project("muni-rendering-toolchain")

-- version
set_version("0.1")
set_xmakever("2.8.0")

-- set warning all as error
set_warnings("all")

-- rules
add_rules("mode.debug", "mode.release")
add_rules("plugin.compile_commands.autoupdate", {outputdir = "build"})
set_languages("cxx20")

-- ext
add_requires("spdlog 1.13.0")
add_requires("stb 2023.01.30")
add_requires("linalg 2.2")
add_requires("tinyobjloader v2.0.0rc13")
-- targets
target("muni-rendering-toolchain")
    set_kind("headeronly")
    add_includedirs("src/muni", {public = true})
    add_packages("spdlog", {public = true})
    add_packages("stb", {public = true})
    add_packages("linalg", {public = true})
    add_packages("tinyobjloader", {public = true})

target("assignment-4")
    set_kind("binary")
    add_files("src/main.cpp")
    add_deps("muni-rendering-toolchain")
    add_cxflags("-Xclang -fopenmp -I/usr/local/include")
    add_ldflags("-L/opt/homebrew/opt/llvm/lib -lomp")
