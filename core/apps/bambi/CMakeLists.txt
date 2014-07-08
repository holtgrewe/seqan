# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/bambi/CMakeLists.txt
#
# CMakeLists.txt file for bambi.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_apps_bambi)
message (STATUS "Configuring core/apps/bambi")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# Require C++11
find_package (CXX11)
if (NOT CXX11_FOUND)
  message (STATUS "  C++11 not found, skipping.")
  return ()
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (bam2fastq
                bam2fastq/bam2fastq.cpp
                bam2fastq/converter_thread.h
                bam2fastq/job_queue.h
                bam2fastq/sequence_sink.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (bam2fastq ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (bambi)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install bambi in ${PREFIX}/bin directory
install (TARGETS bam2fastq
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/bambi for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (bambi)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:bambi")
  set (CPACK_PACKAGE_NAME "bambi")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "bambi")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")

  seqan_configure_cpack_app (bambi "bambi")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:bambi")

