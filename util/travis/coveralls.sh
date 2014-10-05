#!/bin/bash

# Note that this only works if the tests were built and linked using the
# profiling flags.
if [ "$CXX" == "g++" ];
then
  sudo pip install cpp-coveralls
  # redirect to /dev/null to stay below 4MB of stderr/stdout output
  coveralls -b _build -r _build -e core/tests -e extras/tests -t ${COVERALLS_TOKEN} &>/dev/null
fi
