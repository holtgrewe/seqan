# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/bam_realigner/CMakeLists.txt
#
# CMakeLists.txt file for bam_realigner.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_apps_bam_realigner)
message (STATUS "Configuring extras/apps/bam_realigner")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# enable C++11 support
find_package (CXX11)
if (NOT CXX11_FOUND)
  message (STATUS "  C++11 not found, not building bam_realigner.")
  return ()
endif (NOT CXX11_FOUND)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS_MAIN})
include_directories (SYSTEM ${SEQAN_INCLUDE_DIRS_DEPS})
# Updated C++ flags.
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

add_executable (bam_realigner
                bam_realigner.cpp
                bam_realigner_app.cpp
                bam_realigner_app.h
                bam_realigner_options.h
                bam_realigner_options.cpp
                realigner_step.h
                realigner_step.cpp)
target_link_libraries (bam_realigner ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (bam_realigner)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install bam_realigner in ${PREFIX}/bin directory
install (TARGETS bam_realigner
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/bam_realigner for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (bam_realigner)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:bam_realigner")
  set (CPACK_PACKAGE_NAME "bam_realigner")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "bam_realigner")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")

  seqan_configure_cpack_app (bam_realigner "bam_realigner")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:bam_realigner")

