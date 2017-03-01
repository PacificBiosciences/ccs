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

// Author: Lance Hepler

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace PacBio {
namespace Exception {

class ModelNamingError : public std::runtime_error
{
public:
    ModelNamingError(const std::string& msg) : std::runtime_error(msg) {}
};
}

namespace Consensus {

struct ModelForm
{
    enum : uint8_t
    {
        MARGINAL = 0,
        SNR = 1,
        PWSNRA = 2,
        PWSNR = 3
    };

    ModelForm(const uint8_t u) : v{u}
    {
        if (v > PWSNR) throw std::runtime_error("invalid model form!");
    }

    ModelForm(const std::string& s)
    {
        static const std::map<std::string, uint8_t> m{
            {"Marginal", MARGINAL}, {"Snr", SNR}, {"PwSnrA", PWSNRA}, {"PwSnr", PWSNR}};
        const auto it = m.find(s);
        if (it == m.end()) throw Exception::ModelNamingError("invalid model form: " + s);
        v = it->second;
    }

    bool operator<(const ModelForm o) const { return v < o.v; }

    operator std::string() const
    {
        switch (v) {
            case MARGINAL:
                return "Marginal";
            case SNR:
                return "Snr";
            case PWSNRA:
                return "PwSnrA";
            case PWSNR:
                return "PwSnr";
            default:
                throw std::runtime_error("invalid model form!");
        }
    }

    static std::vector<std::string> Preferences()
    {
        std::vector<std::string> prefs;
        for (uint8_t i = 4; i > 0; --i)
            prefs.emplace_back(ModelForm(i - 1));
        return prefs;
    }

private:
    uint8_t v;
};

struct ModelOrigin
{
    enum : uint8_t
    {
        COMPILED = 0,
        BUNDLED = 1,
        PROVIDED = 2
    };

    ModelOrigin(const uint8_t u) : v{u}
    {
        if (v > PROVIDED) throw std::runtime_error("invalid model origin!");
    }

    ModelOrigin(const std::string& s)
    {
        static const std::map<std::string, uint8_t> m{
            {"Compiled", COMPILED}, {"Bundled", BUNDLED}, {"FromFile", PROVIDED}};
        const auto it = m.find(s);
        if (it == m.end()) throw Exception::ModelNamingError("invalid model origin: " + s);
        v = it->second;
    }

    bool operator<(const ModelOrigin o) const { return v < o.v; }

    operator std::string() const
    {
        switch (v) {
            case COMPILED:
                return "Compiled";
            case BUNDLED:
                return "Bundled";
            case PROVIDED:
                return "FromFile";
            default:
                throw std::runtime_error("invalid model origin!");
        }
    }

    static std::vector<std::string> Preferences()
    {
        std::vector<std::string> prefs;
        for (uint8_t i = 3; i > 0; --i)
            prefs.emplace_back(ModelOrigin(i - 1));
        return prefs;
    }

private:
    uint8_t v;
};

struct ModelName : public std::tuple<std::string, ModelForm, ModelOrigin>
{
    ModelName(const std::string& chemistry, const ModelForm form, const ModelOrigin origin)
        : std::tuple<std::string, ModelForm, ModelOrigin>{chemistry, form, origin}
    {
    }
    ModelName(const std::string& s) : std::tuple<std::string, ModelForm, ModelOrigin>(FromString(s))
    {
    }

    operator std::string() const
    {
        static const auto d = "::";
        return std::get<0>(*this) + d + std::string(std::get<1>(*this)) + d +
               std::string(std::get<2>(*this));
    }

private:
    static std::tuple<std::string, ModelForm, ModelOrigin> FromString(const std::string& s)
    {
        static const auto d = "::";
        const size_t fst = s.find(d);
        if (fst == std::string::npos)
            throw Exception::ModelNamingError("invalid model name: " + s);
        std::string chemistry = s.substr(0, fst);
        const size_t snd = s.find(d, fst + 2);
        if (snd == std::string::npos)
            throw Exception::ModelNamingError("invalid model name: " + s);
        ModelForm form(s.substr(fst + 2, snd - (fst + 2)));
        ModelOrigin origin(s.substr(snd + 2));
        return std::make_tuple(std::move(chemistry), std::move(form), std::move(origin));
    }
};
}
}
