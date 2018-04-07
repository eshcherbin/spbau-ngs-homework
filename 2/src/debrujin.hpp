#ifndef __DEBRUJIN_H__
#define __DEBRUJIN_H__

#include <vector>
#include <string>

using namespace std;

// --------------
// | structures |
// --------------

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

// ------------------
// | building graph |
// ------------------

int add_node(const string&);

void add_edge_pair(int, int, const string &);

void process_one_read(const string &);

// ----------------------
// | contracting graph  |
// ----------------------

void remove_redundant_nodes();

// ------------------
// | simplification |
// ------------------

void remove_tips();

void remove_bad_qual_edges();

void remove_bulges();

// ----------
// | output |
// ----------

void write_fasta();

void write_dot();

#endif
