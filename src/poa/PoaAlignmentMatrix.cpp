// Author: David Alexander

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "PoaAlignmentMatrix.h"
#include "PoaGraphImpl.h"

using namespace PacBio::Align;

namespace PacBio {
namespace Poa {
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
    std::ostringstream header0;  // Vertex IDs
    std::ostringstream header1;  // Vertex labels
    std::vector<std::ostringstream> outputRows(NumRows());

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
}  // PacBio::Poa::detail
