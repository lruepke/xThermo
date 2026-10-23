#include "stdfunc.h"
namespace xThermal{
    std::vector<std::string> string_split (const std::string& s, const std::string& delimiter)
    {
        size_t pos_start = 0, pos_end, delim_len = delimiter.length();
        std::string token;
        std::vector<std::string> res;
        while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos)
        {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
        }
        res.push_back (s.substr (pos_start));
        return res;
    }

    std::string extname_file(const std::string& filepath)
    {
        std::string extname;
        std::vector<std::string> tmp=string_split(filepath,".");
        if(!tmp.empty())
        {
            extname=tmp[tmp.size()-1];
        }
        return extname;
    }

    std::string filename_without_ext(const std::string& filepath)
    {
        std::string basename;
        std::vector<std::string> tmp=string_split(filepath,".");
        if(!tmp.empty())
        {
            basename = tmp[0];
            for (int i = 1; i < tmp.size()-1; ++i) {
                basename +="."+tmp[i];
            }
        } else
        {
            basename = filepath;
        }
        return basename;
    }

}