set shell := ["bash", "-c"]

default: help

alias h := help
help:
    just -l -f {{justfile()}}

prefix := justfile_directory()
install_dir := (prefix) / "install"

alias s := setup
setup:
    cd {{prefix}}; meson setup --reconfigure --prefix {{install_dir}} build

alias b := build
build:
    cd {{prefix}}; meson compile -C build
    cd {{prefix}}; meson install -C build

alias f := format
# format c++ codes to Mozilla styles
format:
    cd {{prefix}}; clang-format -style=Mozilla -i src/*.cc
