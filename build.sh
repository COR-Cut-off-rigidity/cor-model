#!/bin/sh

BASEDIR=$(dirname "$0")

cmake -S "$BASEDIR" -B "$BASEDIR/build"
cmake --build "$BASEDIR/build"
