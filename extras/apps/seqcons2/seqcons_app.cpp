// ==========================================================================
//                                 SeqCons
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#include "seqcons_app.h"

#include <cctype>
#include <iostream>

#include <seqan/realign.h>
#include <seqan/sequence.h>
#include <seqan/store.h>

#include "seqcons_options.h"

namespace {

// ---------------------------------------------------------------------------
// Function endsWithIgnoreCase()
// ---------------------------------------------------------------------------

bool endsWithIgnoreCase(std::string str, std::string suffix)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(), ::tolower);
    return seqan::endsWith(str, suffix);
}

}

// ----------------------------------------------------------------------------
// Class SeqConsAppImpl
// ----------------------------------------------------------------------------

class SeqConsAppImpl
{
public:
    SeqConsAppImpl(SeqConsOptions const & options) : options(options)
    {}

    void run();

private:
    // Load reads without alignment information from FASTA file.  Will create one pseudo-contig and put all reads at the
    // first position.
    void loadReads(char const * fileName);
    // Load alignments from SAM file, will set reference to string of N of sufficient length.
    void loadAlignments(char const * fileName);

    // Perform the consensus alignment, optionally interpreting coordinates.
    void performConsensusAlignment(bool useContigID, bool usePositions);
    // Perform realignment on store only.
    void performRealignment();

    // The fragment store for the data.
    seqan::FragmentStore<> store;

    // Configuration.
    SeqConsOptions options;
};

void SeqConsAppImpl::run()
{
    // Load read or alignment data.
    if (options.verbosity >= 1)
        std::cout << "\n__LOADING DATA_______________________________________________________________\n"
                  << '\n';
    if (endsWithIgnoreCase(options.inputFile, ".sam"))
        loadAlignments(options.inputFile.c_str());
    else  // sequence file
        loadReads(options.inputFile.c_str());

    if (options.verbosity >= 1)
        std::cout << "\n__COMPUTATION________________________________________________________________\n"
                  << '\n';
    // Perform the consensus or realignment computation.
    switch (options.operation)
    {
        case SeqConsOptions::MSA_CONSENSUS:
            performConsensusAlignment(false, false);
            break;

        case SeqConsOptions::CTG_CONSENSUS:
            performConsensusAlignment(true, false);
            break;

        case SeqConsOptions::POS_CONSENSUS:
            performConsensusAlignment(true, true);
            break;

        case SeqConsOptions::REALIGN:
            performRealignment();
            break;

        case SeqConsOptions::NOP:
        default:
            // do nothing, will just write out store
            break;
    }
}

void SeqConsAppImpl::loadReads(char const * fileName)
{
}

void SeqConsAppImpl::loadAlignments(char const * fileName)
{
}

void SeqConsAppImpl::performConsensusAlignment(bool useContigID, bool usePositions)
{
    // Setup the consensus alignment options.
    seqan::ConsensusAlignmentOptions caOptions;
    caOptions.useContigID = useContigID;
    caOptions.usePositions = usePositions;
    caOptions.runRealignment = false;  // will run manually.

    caOptions.overlapMaxErrorRate = options.overlapMaxErrorRate;
    caOptions.overlapMinLength = options.overlapMinLength;
    caOptions.overlapMinCount = options.overlapMinCount;
    caOptions.posDelta = options.overlapWindowSize;

    caOptions.kMerSize = options.kMerSize;
    caOptions.kMerMaxOcc = options.kMerMaxOcc;

    // Perform the consensus alignment.
    double startTime = seqan::sysTime();
    if (options.verbosity >= 1)
        std::cerr << "Performing consensus computation...\n";
    consensusAlignment(store, caOptions);
    if (options.verbosity >= 2)
        std::cerr << "\t=> consensus step took " << seqan::sysTime() - startTime << "s\n";

    // Finally, compute realignment of the resulting multi-read alignment.
    performRealignment();
}

void SeqConsAppImpl::performRealignment()
{
    double startTime = seqan::sysTime();
    if (options.verbosity >= 1)
        std::cerr << "Performing realignment...\n";
    for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
        reAlignment(store, contigID, /*method=*/2, options.reAlignmentBandwidth, /*includeReference=*/false);
    if (options.verbosity >= 2)
        std::cerr << "\t=> realignment took " << seqan::sysTime() - startTime << "s\n";
}

// ----------------------------------------------------------------------------
// Class SeqConsApp
// ----------------------------------------------------------------------------

SeqConsApp::SeqConsApp(SeqConsOptions const & options) : impl(new SeqConsAppImpl(options))
{}

SeqConsApp::~SeqConsApp()
{}

void SeqConsApp::run()
{
    impl->run();
}