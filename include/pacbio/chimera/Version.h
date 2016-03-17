#include <string>

#define CHIMERA_GIT_SHA1 "c07a339"
#define CHIMERA_VERSION  "0.0.1"

namespace PacBio {
namespace Chimera {

struct Version
{
    static const size_t Major = 0;
    static const size_t Minor = 0;
    static const size_t Patch = 1;

    static std::string ToString()
    { return CHIMERA_VERSION; }
};

}  // namespace Chimera
}  // namespace PacBio
