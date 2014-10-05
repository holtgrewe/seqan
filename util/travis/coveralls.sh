#!/bin/bash

# Note that this only works if the tests were built and linked using the
# profiling flags.
if [ "$CXX" == "g++" ];
then
  sudo pip install cpp-coveralls
  # redirect to /dev/null to stay below 4MB stdout?
  coveralls -b . -r . -E '.*CMake.*\.cpp' -E '.*CMake.*\.c' -t "${COVERALLS_TOKEN}" #&>/dev/null
fi
