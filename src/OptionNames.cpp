// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.


#include <pacbio/ccs/OptionNames.h>

using namespace PacBio::CCS;

namespace PacBio {
namespace CCS {

	const std::string OPTION_NAMES::ZMWS = "zmws";
	const std::string OPTION_NAMES::MIN_SNR = "minSnr";
	const std::string OPTION_NAMES::MIN_READ_SCORE ="minReadScore";
	const std::string OPTION_NAMES::REPORT_FILE = "reportFile";
	const std::string OPTION_NAMES::NUM_THREADS = "numThreads";
	const std::string OPTION_NAMES::LOG_FILE = "logFile";
	const std::string OPTION_NAMES::LOG_LEVEL = "logLevel";
	const std::string OPTION_NAMES::MIN_LENGTH = "minLength";
	const std::string OPTION_NAMES::MIN_PASSES = "minPasses";
	const std::string OPTION_NAMES::MIN_PREDICTED_ACCURACY = "minPredictedAccuracy";
	const std::string OPTION_NAMES::MIN_Z_SCORE = "minZScore";
	const std::string OPTION_NAMES::MAX_DROP_FRAC = "maxDropFrac";


} // namespace CCS
} // namespace PacBio
