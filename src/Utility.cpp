// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

#include <climits>
#include <cstdlib>
#include <fstream>

#include <sys/stat.h>

#include <boost/algorithm/string.hpp>

#include <pacbio/ccs/Utility.h>

using namespace std;

namespace PacBio {
namespace CCS {

string AbsolutePath(const string& path)
{
    char buf[PATH_MAX + 1];
    char *res = realpath(path.c_str(), buf);

    if (res != nullptr)
    {
        string abspath(res);
        return abspath;
    }

    return path;
}

bool FileExists(const string& path)
{
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
}

void FlattenFofn(vector<string>& res, const string& file)
{
    using boost::algorithm::iends_with;
    using boost::algorithm::trim;

    if (iends_with(file, ".fofn"))
    {
        ifstream fofn(file);
        string line;
        while (getline(fofn, line))
        {
            trim(line);
            FlattenFofn(res, line);
        }
    }
    else if (iends_with(file, ".bam"))
        res.push_back(file);
    else
        throw invalid_argument("not a .fofn or .bam file!");
}

vector<string> FlattenFofn(const vector<string>& files)
{
    vector<string> res;

    for (const auto& file : files)
    {
        FlattenFofn(res, file);
    }

    return res;
}

} // namespace CCS
} // namespace PacBio
