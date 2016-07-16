#!/usr/bin/env bash -x
for f in `find include src tests -iname '*.h' -o -iname '*.cpp'`
do
    clang-format -i $f
done
