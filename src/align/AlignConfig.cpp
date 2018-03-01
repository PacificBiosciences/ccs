// Authors: David Alexander, Lance Hepler

#include <pacbio/align/AlignConfig.h>

namespace PacBio {
namespace Align {

AlignParams::AlignParams(int match, int mismatch, int insert, int delete_)
    : Match(match), Mismatch(mismatch), Insert(insert), Delete(delete_)
{
}

AlignParams AlignParams::Default() { return {0, -1, -1, -1}; }

AlignConfig::AlignConfig(AlignParams params, AlignMode mode) : Params(params), Mode(mode) {}

AlignConfig AlignConfig::Default() { return {AlignParams::Default(), AlignMode::GLOBAL}; }

}  // namespace Align
}  // namespace PacBio
