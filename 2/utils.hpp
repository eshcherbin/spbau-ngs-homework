#ifndef __UTILS_H__
#define __UTILS_H__

#include <algorithm>
#include <string>
#include <unordered_map>

using namespace std;

unordered_map<char, char> COMP = {
    {'A', 'T'}, 
    {'T', 'A'}, 
    {'G', 'C'}, 
    {'C', 'G'}
};

inline string rev_comp(const string& s) {
    string res(s);
    transform(res.begin(), res.end(), res.begin(), [&](char c) { return COMP[c]; });
    reverse(res.begin(), res.end());
    return res;
}

inline int comp_node(int v) {
    return v ^ 1;
}

#endif
