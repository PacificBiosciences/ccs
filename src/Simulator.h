// Author: David Seifert

#pragma once

#include <cassert>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>

namespace PacBio {
namespace Consensus {
namespace {

/// The simulator completes the cycle
///
///   .--> parameter inference from data ---.
///   |                                     |
///   '-- data simulation from parameters <-'
///
/// The simulator model has been implemented in its basic form by
/// Nigel Delaney and Yuan Tian in unitem. The implementation here
/// is in order to consolidate the code base.

/// SimulateReadImpl is the main function for generating reads from the
/// different chemistry models. In order to draw reads, the factory
/// function of the ModelConfig class can be used. For setting up
/// simulation of reads, every chemistry models needs to provide an
/// explicit specialization of
///   1. InitialiseModel
///   2. GenerateReadData

struct BaseData
{
    char base;
    uint8_t pw;
    uint8_t ipd;
};

template <typename TInitialiseModel, typename TGenerateReadData>
std::pair<Data::Read, std::vector<MoveType>> SimulateReadImpl(std::default_random_engine* const rng,
                                                              const std::string& tpl,
                                                              const std::string& readname,
                                                              TInitialiseModel initialiseModel,
                                                              TGenerateReadData generateReadData)
{
    if (tpl.empty()) throw std::invalid_argument("Cannot pass empty template to the Simulator!");

    std::string readBases;
    std::vector<uint8_t> readPw, readIpd;
    std::vector<MoveType> statePath;

    MoveType state = MoveType::MATCH;
    uint8_t pw, ipd;

    const std::pair<Data::SNR, std::vector<TemplatePosition>> snrTransModel{
        initialiseModel(rng, tpl)};
    const Data::SNR& snrs = snrTransModel.first;
    const std::vector<TemplatePosition>& transModel = snrTransModel.second;

    size_t locus = 0;
    while (locus < tpl.size()) {
        const AlleleRep prev = (locus ? transModel[locus - 1].Idx : AlleleRep::FromASCII('A'));
        const AlleleRep curr = transModel[locus].Idx;

        // the first base MUST be a match
        if (locus == 0)
            state = MoveType::MATCH;
        else {
            // 1. generate new state
            std::discrete_distribution<uint8_t> stateDistrib(&transModel[locus - 1].Match,
                                                             &transModel[locus - 1].Deletion + 1);
            state = static_cast<MoveType>(stateDistrib(*rng));
        }
        statePath.push_back(state);

        switch (state) {
            case MoveType::MATCH:
                ++locus;
            case MoveType::BRANCH:
            case MoveType::STICK:
                break;
            case MoveType::DELETION:
                ++locus;
                continue;
        }

        // 2. generate new base data (base + PW + IPD)
        assert(state != MoveType::DELETION);
        const BaseData generatedBase = generateReadData(rng, state, prev, curr);

        readBases.push_back(generatedBase.base);
        readPw.push_back(generatedBase.pw);
        readIpd.push_back(generatedBase.ipd);
    }

    return {{readname, readBases, readIpd, readPw, snrs, "simulate"}, std::move(statePath)};
}
}
}  // namespace Consensus
}  // namespace PacBio
