// Author: Armin Töpfer

#pragma once

#include <string>
#include <vector>

#include <pbcopper/cli/CLI.h>

namespace PacBio {
namespace Data {
struct PlainOption
{
    std::string id;
    std::vector<std::string> cliOptions;
    std::string name;
    std::string description;
    JSON::Json defaultValue;
    JSON::Json choices = JSON::Json(nullptr);
    CLI::OptionFlags flags;

    PlainOption(const std::string& id, const std::vector<std::string>& cliOptions,
                const std::string& name, const std::string& description,
                const JSON::Json& defaultValue, const JSON::Json& choices = JSON::Json(nullptr),
                const CLI::OptionFlags& flags = CLI::OptionFlags::DEFAULT)
        : id(id)
        , cliOptions(cliOptions)
        , name(name)
        , description(description)
        , defaultValue(defaultValue)
        , choices(choices)
        , flags(flags)
    {
    }

    operator CLI::Option() const
    {
        return {id, cliOptions, description, defaultValue, choices, flags};
    }
    operator std::pair<std::string, std::string>() const { return std::make_pair(id, name); }
    operator std::string() const { return id; }
};
}
}  // :: PacBio::CLI