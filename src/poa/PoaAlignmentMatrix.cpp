// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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

// Author: David Alexander

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "PoaAlignmentMatrix.h"
#include "PoaGraphImpl.h"

namespace PacBio {
namespace Consensus {
namespace detail {

void PoaAlignmentMatrixImpl::Print() const
{
    const int COL_WIDTH = 6;
    const int ROW_HEIGHT = 2;

    const std::map<MoveType, char> moveCode = {
        {InvalidMove, 'Z'},  {StartMove, 'S'},  {EndMove, 'E'},  {MatchMove, '='},
        {MismatchMove, 'X'}, {DeleteMove, 'D'}, {ExtraMove, 'I'}};

    // Display something like this:
    //
    //           0     2     3     4     5     6     7     8     9     1
    //           ^     A     C     G     T     A     C     G     T     $
    //
    //     -    0Z    0S    0S    0S                               -infZ
    //
    //     A    0S    3=    0S    0S    0S                         -infZ
    //
    //     C    0S    0S    6=    2D    0S    0S                   -infZ
    //
    //     G    0S    0S    2I    9=    5D    1D    0S             -infZ
    //
    //     T    0S    0S    0S    5I   12=    8D    4D    0S       -infZ
    //
    //     A    0S    3=    0S    1I    8I   15=   11D    7D    3D -infZ
    //
    //     C    0S    0S    6=    2D    4I   11I   18=   14D   10D -infZ
    //
    //     G          0S    2I    9=    5D    7I   14I   21=   17D -infZ
    //
    std::stringstream header0;  // Vertex IDs
    std::stringstream header1;  // Vertex labels
    std::vector<std::stringstream> outputRows(NumRows());

    header0 << std::setw(COL_WIDTH) << std::right << "";
    header1 << std::setw(COL_WIDTH) << std::right << "";
    for (size_t row = 0; row < NumRows(); row++) {
        char readBase;
        if (row == 0) {
            readBase = '-';
        } else {
            readBase = readSequence_[row - 1];
        }
        outputRows[row] << std::setw(COL_WIDTH) << std::right << readBase;
    }

    for (const auto& v : graph_->sortedVertices()) {
        const AlignmentColumn* col = columns_.at(v);
        const PoaNode& node = graph_->getPoaNode(v);

        header0 << std::setw(COL_WIDTH) << std::right << node.Id;
        header1 << std::setw(COL_WIDTH) << std::right << node.Base;

        for (size_t j = 0; j < NumRows(); j++) {
            if (col->HasRow(j)) {
                float score = col->Score[j];
                string scoreFmt = score == -FLT_MAX ? "-inf" : std::to_string(int(score));
                string cellFmt = scoreFmt + moveCode.at(col->ReachingMove[j]);
                outputRows[j] << std::setw(COL_WIDTH) << std::right << cellFmt;
            } else {
                outputRows[j] << std::setw(COL_WIDTH) << std::right << "";
            }
        }
    }

    std::cout << header0.str() << std::endl;
    std::cout << header1.str() << std::endl;
    for (size_t row = 0; row < NumRows(); row++) {
        for (size_t skip = ROW_HEIGHT; skip > 1; skip--) {
            std::cout << std::endl;
        }
        std::cout << outputRows[row].str() << std::endl;
    }
}
}
}
}  // PacBio::Consensus::detail
