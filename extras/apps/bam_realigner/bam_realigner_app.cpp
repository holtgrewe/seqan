// ==========================================================================
//                               BAM Realigner
// ==========================================================================
// Copyright (c) 2014, Manuel Holtgrewe
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "bam_realigner_app.h"

#include <iostream>
#include <vector>
#include <string>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/simple_intervals_io.h>

#include "bam_realigner_options.h"
#include "realigner_step.h"

namespace {  // anonymous namespace

}  // anonymous namespace

// ---------------------------------------------------------------------------
// Class BamRealignerAppImpl
// ---------------------------------------------------------------------------

class BamRealignerAppImpl
{
public:
    BamRealignerAppImpl(BamRealignerOptions const & options) :
            options(options), bamFileIn(bamFileOut)
    {}

    void run();

private:

    // Open FASTA index file.
    void openFai();
    // Open input BAM file and bai index.
    void openBamIn();
    // Open intervals file.
    void openIntervals();

    // Open output BAM file.
    void openBamOut();
    // Open output MSA txt file.
    void openMsasTxtOut();

    // Process regions one-by-one.
    void processAllRegions();
    void processOneRegion(seqan::GenomicRegion const & region);

    // Program configuration.
    BamRealignerOptions options;

    // Objects used for I/O.
    seqan::BamFileOut bamFileOut;
    seqan::VirtualStream<char, seqan::Output> msasTxtOut;
    seqan::FaiIndex faiIndex;
    seqan::BamFileIn bamFileIn;
    seqan::BamIndex<seqan::Bai> baiIndex;
    seqan::SimpleIntervalsFileIn intervalsFileIn;

    // BAM header is read into this variable.
    seqan::BamHeader bamHeader;
};

void BamRealignerAppImpl::run()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "BAM Realigner\n"
                  << "=============\n\n";
        options.print(std::cerr);
    }

    // Open Input Files

    if (options.verbosity >= 1)
    {
        std::cerr << "\n"
                  << "__OPENING INPUT FILES____________________________________________\n"
                  << "\n";
    }

    openFai();
    openBamIn();
    openIntervals();

    // Open Input Files

    if (options.verbosity >= 1)
    {
        std::cerr << "\n"
                  << "__OPENING OUTPUT FILES____________________________________________\n"
                  << "\n";
    }

    openBamOut();
    openMsasTxtOut();

    // Process Intervals

    processAllRegions();

    // Writing Output
}

void BamRealignerAppImpl::processAllRegions()
{
    std::cerr << "\n"
              << "__PROCESSING REGIONS_____________________________________________\n"
              << "\n";

    seqan::GenomicRegion region;
    for (unsigned no = 1; !atEnd(intervalsFileIn); ++no)
    {
        readRecord(region, intervalsFileIn);
        seqan::CharString buffer;
        region.toString(buffer);
        // std::cerr << "\r                                                   "
        //           << "\rProcessing (#" << no << ") " << buffer << std::flush;
        std::cerr << "Processing (#" << no << ") " << buffer << "\n";

        processOneRegion(region);
    }

    std::cerr << " DONE\n";
}

void BamRealignerAppImpl::processOneRegion(seqan::GenomicRegion const & region)
{
    RealignerStep worker(bamFileOut, msasTxtOut, bamFileIn, baiIndex, faiIndex, region, options);
    worker.run();
}

void BamRealignerAppImpl::openFai()
{
    if (options.verbosity >= 1)
        std::cerr << "    Opening " << options.inReferencePath << " (using FAI index) ...";
    if (!open(faiIndex, options.inReferencePath.c_str()))
    {
        std::cerr << " (building .fai)";
        if (!build(faiIndex, options.inReferencePath.c_str()))
            throw seqan::IOError("Could not build .fai index.");
        if (!save(faiIndex))
            throw seqan::IOError("Could not save .fai index.");
    }
    if (options.verbosity >= 1)
        std::cerr << "OK\n";
}

void BamRealignerAppImpl::openBamIn()
{
    if (options.verbosity >= 1)
        std::cerr << "    Opening " << options.inAlignmentPath << " ...";
    if (!open(bamFileIn, options.inAlignmentPath.c_str()))
        throw seqan::IOError("Could not open BAM file.");
    if (options.verbosity >= 1)
        std::cerr << " OK\n";

    if (options.verbosity >= 1)
        std::cerr << "        Reading header ...";
    readRecord(bamHeader, bamFileIn);
    if (options.verbosity >= 1)
        std::cerr << " OK\n";

    std::string baiPath = options.inAlignmentPath + ".bai";
    if (options.verbosity >= 1)
        std::cerr << "    Opening " << baiPath << " ...";
    if (!open(baiIndex, baiPath.c_str()))
        throw seqan::IOError("Could not open BAI file.");
    if (options.verbosity >= 1)
        std::cerr << "OK\n";
}

void BamRealignerAppImpl::openIntervals()
{
    if (options.verbosity >= 1)
        std::cerr << "    Opening " << options.inIntervalsPath << " ...";
    if (!open(intervalsFileIn, options.inIntervalsPath.c_str()))
        throw seqan::IOError("Could not open intervals file.");
    if (options.verbosity >= 1)
        std::cerr << "OK\n";
}

void BamRealignerAppImpl::openBamOut()
{
    if (options.verbosity >= 1)
        std::cerr << "    Opening " << options.outAlignmentPath << " ...";
    if (!open(bamFileOut, options.outAlignmentPath.c_str()))
        throw seqan::IOError("Could not open output BAM file.");
    if (options.verbosity >= 1)
        std::cerr << "OK\n";
    writeRecord(bamFileOut, bamHeader);
}

void BamRealignerAppImpl::openMsasTxtOut()
{
    if (options.outMsasPath.empty())
        return;

    if (options.verbosity >= 1)
        std::cerr << "    Opening " << options.outMsasPath << " ...";
    if (!open(msasTxtOut, options.outMsasPath.c_str()))
        throw seqan::IOError("Could not open output MSA output txt file.");
    if (options.verbosity >= 1)
        std::cerr << "OK\n";
}


// ---------------------------------------------------------------------------
// Class BamRealignerApp
// ---------------------------------------------------------------------------

BamRealignerApp::BamRealignerApp(BamRealignerOptions const & options) :
        impl(new BamRealignerAppImpl(options))
{}

BamRealignerApp::~BamRealignerApp()  // for pimpl
{}

void BamRealignerApp::run()
{
    impl->run();
}
