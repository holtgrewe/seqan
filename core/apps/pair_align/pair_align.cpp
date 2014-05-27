/*==========================================================================
               SeqAn - The Library for Sequence Analysis
                         http://www.seqan.de 
============================================================================
Copyright (C) 2007-2012

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
==========================================================================*/

#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>

#include <iostream>
#include <fstream>


using namespace seqan;

// --------------------------------------------------------------------------
// Class Options
// --------------------------------------------------------------------------

struct Options
{
    static int const INVALID_DIAGONAL;

    seqan::CharString inputFile;
    seqan::CharString outputFile;
    seqan::CharString alphabet;
    seqan::CharString method;
    int outputFormat;
    int gop;
    int gex;
    seqan::CharString matrix;
    int msc;
    int mmsc;
    int low;
    int high;
    seqan::CharString config;

    Options() : gop(0), gex(0), msc(0), mmsc(0), low(INVALID_DIAGONAL), high(INVALID_DIAGONAL)
    {}
};

int const Options::INVALID_DIAGONAL = seqan::MaxValue<int>::VALUE;

//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////

template<typename TScore, typename TSc>
inline void
_setMatchScore(TScore&, TSc) {
    // No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TScore, typename TSc>
inline void
_setMismatchScore(TScore&, TSc) {
    // No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TSc>
inline void
_setMatchScore(Score<int, Simple>& sc, TSc msc) {
    sc.data_match = msc;
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TSc>
inline void
_setMismatchScore(Score<int, Simple>& sc, TSc mmsc) {
    sc.data_mismatch = mmsc;
}


//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore>
inline void
_initAlignParams(Options const & options, TScore& sc) {
    // Set options
    sc.data_gap_open = options.gop;
    sc.data_gap_extend = options.gex;
    int msc = options.msc;
    _setMatchScore(sc, msc);
    int mmsc = options.mmsc;
    _setMismatchScore(sc, mmsc);
    ::std::string seqfile = toCString(options.inputFile);
    ::std::string outfile = toCString(options.outputFile);
    unsigned int method = 0;
    String<char> meth = options.method;
    if (meth == "nw") method = 0;
    else if (meth == "gotoh") method = 1;
    else if (meth == "sw") method = 2;
    else if (meth == "lcs") method = 3;
    int low = 0;
    int high = 0;
    bool banded = false;
    if (options.low != Options::INVALID_DIAGONAL)
    {
        low = options.low;
        banded = true;
    }
    if (options.high != Options::INVALID_DIAGONAL)
    {
        high = options.high;
        banded = true;
    }

    // Check options
    if (low > high) banded = false;
    
    // Do pairwise alignment
    String<char> config = options.config;
    if (!empty(config))
    {
        if (config == "tttt")
            performPairwiseAlignment<TAlphabet, AlignConfig<true, true, true, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "tttf")
            performPairwiseAlignment<TAlphabet, AlignConfig<true, true, true, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ttft")
            performPairwiseAlignment<TAlphabet, AlignConfig<true, true, false, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ttff")
            performPairwiseAlignment<TAlphabet, AlignConfig<true, true, false, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "tftt")
            performPairwiseAlignment<TAlphabet, AlignConfig<true, false, true, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "tftf")
            performPairwiseAlignment<TAlphabet, AlignConfig<true, false, true, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "tfft")
            performPairwiseAlignment<TAlphabet, AlignConfig<true, false, false, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "tfff")
            performPairwiseAlignment<TAlphabet, AlignConfig<true, false, false, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "fttt")
            performPairwiseAlignment<TAlphabet, AlignConfig<false, true, true, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "fttf")
            performPairwiseAlignment<TAlphabet, AlignConfig<false, true, true, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ftft")
            performPairwiseAlignment<TAlphabet, AlignConfig<false, true, false, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ftff")
            performPairwiseAlignment<TAlphabet, AlignConfig<false, true, false, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "fftt")
            performPairwiseAlignment<TAlphabet, AlignConfig<false, false, true, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "fftf")
            performPairwiseAlignment<TAlphabet, AlignConfig<false, false, true, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ffft")
            performPairwiseAlignment<TAlphabet, AlignConfig<false, false, false, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ffff")
            performPairwiseAlignment<TAlphabet, AlignConfig<false, false, false, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
    }
    else
    {
        performPairwiseAlignment<TAlphabet, AlignConfig<false, false, false, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
    }
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(Options const & options, Dna5 const) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, options.matrix);
        _initAlignParams<Dna5>(options, sc);
    }
    else
    {
        Score<int> sc;
        _initAlignParams<Dna5>(options, sc);
    }
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(Options const & options, char const) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, options.matrix);
        _initAlignParams<char>(options, sc);
    }
    else
    {
        Score<int> sc;
        _initAlignParams<char>(options, sc);
    }
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(Options const & options, Rna5 const) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {    
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, options.matrix);
        _initAlignParams<Rna5>(options, sc);
    }
    else
    {
        Score<int> sc;
        _initAlignParams<Rna5>(options, sc);
    }
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(Options const & options, AminoAcid const) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, options.matrix);
        _initAlignParams<AminoAcid>(options, sc);
    }
    else
    {
        Blosum62 sc;
        _initAlignParams<AminoAcid>(options, sc);
    }
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("pair_align");
    // Set short description, version, and date.
    setShortDescription(parser, "Pairwise alignment");
    setVersion(parser, "1.1");
    setDate(parser, "November 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-s\\fP \\fIIN.fa\\fP");
	setCategory(parser, "Sequence Alignment");
    addDescription(parser,
                   "The program allows to align two sequences using dyamic programming alignment algorithms while "
                   "tweaking various parameters.");

    addSection(parser, "Main Options");
    addOption(parser, seqan::ArgParseOption("s", "seq", "FASTA file with two sequences.", seqan::ArgParseOption::INPUTFILE, "IN.fa"));
    setRequired(parser, "seq");
    setValidValues(parser, "seq", "fasta fa");
    addOption(parser, seqan::ArgParseOption("a", "alphabet", "Sequence alphabet.", seqan::ArgParseOption::STRING, "ALPHABET"));
    setValidValues(parser, "alphabet", "protein dna rna text");
    setDefaultValue(parser, "alphabet", "protein");
    addOption(parser, seqan::ArgParseOption("m", "method",
                                            "DP alignment method: Needleman-Wunsch, Gotoh, Smith-Waterman, "
                                            "Longest Common Subsequence",
                                            seqan::ArgParseOption::STRING, "METHOD"));
    setValidValues(parser, "method", "nw gotoh sw lcs");
    setDefaultValue(parser, "method", "gotoh");
    addOption(parser, seqan::ArgParseOption("o", "outfile", "Output filename.", seqan::ArgParseOption::OUTPUTFILE, "OUT"));
    setDefaultValue(parser, "outfile", "out.fasta");
	setValidValues(parser, "outfile", "fa fasta msf");
	//TODO(rmaerker): We removed this option. The file format is derived from the outfile format.
    //addOption(parser, seqan::ArgParseOption("f", "format", "Output format.", seqan::ArgParseOption::STRING));
    //setValidValues(parser, "format", "fa fasta msf");
    //setDefaultValue(parser, "format", "fasta");

    addSection(parser, "Scoring Options");
    addOption(parser, seqan::ArgParseOption("g", "gop", "Gap open penalty.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "gop", "-11");
    addOption(parser, seqan::ArgParseOption("e", "gex", "Gap extension penalty.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "gex", "-1");
    addOption(parser, seqan::ArgParseOption("ma", "matrix", "Score matrix.", seqan::ArgParseOption::STRING, "MATRIX_FILE"));
    addOption(parser, seqan::ArgParseOption("ms", "msc", "Match score.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "msc", "5");
    addOption(parser, seqan::ArgParseOption("mm", "mmsc", "Mismatch penalty.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "mmsc", "-4");

    addSection(parser, "Banded Alignment Options");
    addOption(parser, seqan::ArgParseOption("lo", "low", "Lower diagonal.", seqan::ArgParseOption::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("hi", "high", "Upper diagonal.", seqan::ArgParseOption::INTEGER, "INT"));

    addSection(parser, "DP Matrix Configuration Options");
    addOption(parser, seqan::ArgParseOption("c", "config", "Alignment configuration.", seqan::ArgParseOption::STRING, "CONF"));
    setValidValues(parser, "config", "ffff ffft fftf fftt ftff ftft fttf fttt tfff tfft tftf tftt ttff ttft tttf tttt");

    addTextSection(parser, "Alignment configuration");
    addText(parser,
            "The alignment configuration is a string of four characters, each being either t or f. All "
            "combinations are allowed. The meaning is as follows.");
    addListItem(parser, "tfff", "First row initialized with 0s.");
    addListItem(parser, "ftff", "First column initialized with 0s.");
    addListItem(parser, "fftf", "Search last column for maximum.");
    addListItem(parser, "ffft", "Search last row for maximum.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getOptionValue(options.inputFile, parser, "seq");
    getOptionValue(options.outputFile, parser, "outfile");
	// Guess file format based on extension of file.
	CharString tmp = options.outputFile;
	if (endsWith(tmp, ".fa") || endsWith(tmp, "fasta"))
	    options.outputFormat = 0;
	else if (endsWith(tmp, ".msf"))
	    options.outputFormat = 1;

    getOptionValue(options.alphabet, parser, "alphabet");
    getOptionValue(options.method, parser, "method");

    //getOptionValue(options.format, parser, "format");
    getOptionValue(options.gop, parser, "gop");
    getOptionValue(options.gex, parser, "gex");
    getOptionValue(options.matrix, parser, "matrix");
    getOptionValue(options.msc, parser, "msc");
    getOptionValue(options.mmsc, parser, "mmsc");
    getOptionValue(options.low, parser, "low");
    getOptionValue(options.high, parser, "high");
    getOptionValue(options.config, parser, "config");

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, const char *argv[])
{
    // Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    
    // Initialize scoring matrices
    if (options.alphabet == "dna") _initScoreMatrix(options, Dna5());
    else if (options.alphabet == "rna") _initScoreMatrix(options, Rna5());
    else if (options.alphabet == "protein") _initScoreMatrix(options, AminoAcid());
    else _initScoreMatrix(options, char());

    return 0;
}
