# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /manual/CMakeLists.txt
#
# CMakeLists.txt file for manual.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_manual)

# ---------------------------------------------------------------------------
# Look for Sphinx and balk out if we cannot find it.
find_package (Sphinx)

if (NOT SPHINX_FOUND)
    message (STATUS "  Sphinx not found, not building sphinx manual as a test.")
    return ()
endif (NOT SPHINX_FOUND)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Look for Python and stop if it could not be found.

find_package (PythonInterp)

if (NOT PYTHONINTERP_FOUND)
  message (STATUS "  Python not found, not building sphinx manual as a test.")
  return ()
endif (NOT PYTHONINTERP_FOUND)
# ---------------------------------------------------------------------------

# If this worked, we can proceed.
message (STATUS "  Adding Sphinx build of manual as test.")

# ---------------------------------------------------------------------------
# Create test that builds the documentation.
add_test (NAME build_dox #${CMAKE_BUILD_DIR}/dox/html/js/search.data.js
          COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/util/bin/dox.py
                                       -b ${CMAKE_SOURCE_DIR}/core
                                       -b ${CMAKE_SOURCE_DIR}/extras
                                       -i ${CMAKE_SOURCE_DIR}/core/include/seqan
                                       -i ${CMAKE_SOURCE_DIR}/extras/include/seqan
                                       -i pages
                                       --image-dir ${CMAKE_SOURCE_DIR}/docs2/images)

# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Add building Sphinx documentation as a test.
add_test (NAME build_manual
                  COMMAND ${SPHINX_EXECUTABLE} -W -n -E
                  -D doxlinks_dox_json=${CMAKE_BINARY_DIR}/manual/html/js/search.data.js
                  -b html -d ${CMAKE_BINARY_DIR}/manual/doctrees
                  ${CMAKE_CURRENT_SOURCE_DIR}/source
                  ${CMAKE_BINARY_DIR}/manual/html
                  )
# ---------------------------------------------------------------------------
