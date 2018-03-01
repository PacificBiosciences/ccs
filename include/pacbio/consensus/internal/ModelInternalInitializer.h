// Author: David Seifert

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
