// Author: Armin Töpfer

#pragma once

#include <stdexcept>
#include <string>

#include <pacbio/data/State.h>

namespace PacBio {
namespace Exception {

class StateError : public std::runtime_error
{
public:
    StateError(PacBio::Data::State state, const std::string& msg)
        : std::runtime_error(msg), state_(state)
    {
    }

    PacBio::Data::State WhatState() const { return state_; }
    virtual const char* what() const noexcept override { return std::runtime_error::what(); }
private:
    PacBio::Data::State state_;
};

class TemplateTooSmall : public StateError
{
public:
    TemplateTooSmall() : StateError(PacBio::Data::State::TEMPLATE_TOO_SMALL, "Template too short!")
    {
    }
};

class AlphaBetaMismatch : public StateError
{
public:
    AlphaBetaMismatch()
        : StateError(PacBio::Data::State::ALPHA_BETA_MISMATCH, "Alpha/beta mismatch!")
    {
    }
};

}  // namespace Exception
}  // namespace PacBio
