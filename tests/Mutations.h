
#include <string>
#include <vector>

#include <pacbio/consensus/Mutation.h>

std::vector<PacBio::Consensus::Mutation> Mutations(const std::string& tpl, const size_t start,
                                                   const size_t end);

std::vector<PacBio::Consensus::Mutation> Mutations(const std::string& tpl);
