#include <string>

#define CHIMERA_GIT_SHA1 "46614a9*"
#define CHIMERA_VERSION "0.0.1"

namespace PacBio {
namespace Chimera {

struct Version
{
    static const size_t Major = @PacBioChimera_VERSION_MAJOR @;
    static const size_t Minor = @PacBioChimera_VERSION_MINOR @;
    static const size_t Patch = @PacBioChimera_VERSION_PATCH @;

    static std::string ToString() { return CHIMERA_VERSION; }
};

}  // namespace Chimera
}  // namespace PacBio
