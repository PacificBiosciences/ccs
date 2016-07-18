// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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
        : std::runtime_error(msg)
        , state_(state) 
    {}

    PacBio::Data::State WhatState() const { return state_; }
    virtual const char* what() const noexcept { return std::runtime_error::what(); }
private:
    PacBio::Data::State state_;
};

class TemplateTooSmall : public StateError
{
public:
    TemplateTooSmall() 
        : StateError(PacBio::Data::State::TEMPLATE_TOO_SMALL, "Template too short!") 
    {}
};

class AlphaBetaMismatch : public StateError
{
public:
    AlphaBetaMismatch() 
        : StateError(PacBio::Data::State::ALPHA_BETA_MISMATCH, "Alpha/beta mismatch!") 
    {}
};

class ChemistryNotFound : public std::runtime_error
{
public:
    ChemistryNotFound(const std::string& name)
        : std::runtime_error(std::string("chemistry not found: '") + name + "'")
    {
    }
};

}  // namespace Exception
}  // namespace PacBio
