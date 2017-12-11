// Copyright (c) 2017, Pacific Biosciences of California, Inc.
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

// Authors: David Seifert

#pragma once

namespace PacBio {
namespace Consensus {

// This is a small trick to initialize the internal data structures.
//
// Without this artifical dummy dependency chain of static constructors
// into the main unanimity library, most Unix ld variants would remove
// a bunch of internal initializers. These can be forcibly included by
// using such constructs as `-Wl,--whole-archive` (GNU Binutils) or
// `-force_load` (macOS ld).
//
// Unfortunately, the previous flags make deployment unnecessarily hard,
// as consuming code needs to be aware of these static initializers,
// in order to use the right flag. By using a chain of constructors
// into the library, we can cleanly circumvent this issue, as `ld` will
// be forced to pull in the library's static initializers without
// requiring non-portable toolchain hacks.
//
// See also
// http://www.lysium.de/blog/index.php?/archives/222-Lost-static-objects-in-static-libraries-with-GNU-linker-ld.html

static struct FactoryInit
{
    FactoryInit();
} initFactory;

}  // namespace Consensus
}  // namespace PacBio
