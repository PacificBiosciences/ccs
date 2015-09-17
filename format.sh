#!/usr/bin/env bash
find include src tests -iname '*.h' -o -iname '*.cpp' -exec clang-format -i {} \;
