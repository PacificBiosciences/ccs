// Author: David Alexander

#include "ModelFactory.h"

#include <pacbio/consensus/EasyReadScorer.h>
#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Template.h>
#include <pacbio/data/Read.h>

#include <memory>
#include <string>

namespace PacBio {
namespace Consensus {

class ModelConfig;

Evaluator* EasyReadScorer::MakeEvaluator(const std::string& tplString,
                                         const PacBio::Data::MappedRead& mappedRead,
                                         double minZScore, double scoreDiff)
{
    // Make a "ModelConfig"
    std::unique_ptr<ModelConfig> cfg =
        ModelFactory::Create(mappedRead.Model, mappedRead.SignalToNoise);
    // make a "template" object
    auto tpl = std::make_unique<Template>(tplString, std::move(cfg));

    return new Evaluator(std::move(tpl), mappedRead, minZScore, scoreDiff);
}

}  // PacBio::Consensus
}  // PacBio
