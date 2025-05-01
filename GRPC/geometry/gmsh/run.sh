#!/bin/bash

export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/usr/local/lib

/usr/local/bin/gmsh multigap.geo
