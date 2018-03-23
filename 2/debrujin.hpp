#ifndef __DEBRUJIN_H__
#define __DEBRUJIN_H__

#include <vector>
#include <string>

using namespace std;

struct node {
    string kmer;
    int cnt;
    vector<int> in;
    vector<int> out;

    node(const string &kmer_) : kmer(kmer_), cnt(0), in(), out() {}
};

struct edge {
    int v;
    int u;
    string seq;
    long long sum_kmer_cnt;
};

int add_node(const string&);

int add_edge(int, int, const string&);

void process_one_read(const string&);

void remove_redundant_node(int);

#endif
