#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <cassert>
#include <sstream>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <seqan/seq_io.h>

#include "debrujin.hpp"
#include "config.hpp"

using namespace std;

// --------------------
// | global variables |
// --------------------

int k;
vector<node> nodes;
vector<edge> edges;
unordered_map<string, int> kmer_cnt;
unordered_map<string, int> kmer2node;
set<pair<int, int>> edges_set;

// -----------
// | utility |
// -----------

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

inline int get_node(const string &kmer) {
    if (!kmer2node.count(kmer)) {
        return add_node(kmer);
    }
    return kmer2node[kmer];
}

inline size_t n_nodes() {
    return nodes.size();
}

inline size_t n_edges() {
    return edges.size();
}

inline double get_cov(const edge &e) {
    return ((double) e.sum_kmer_cnt / (e.seq.size() - k));
}

// ------------------
// | building graph |
// ------------------

int add_node(const string& kmer) {
    string rev = rev_comp(kmer);
    nodes.push_back(node(kmer));
    kmer2node[kmer] = n_nodes() - 1;
    nodes.push_back(node(rev));
    kmer2node[rev] = n_nodes() - 1;
    return n_nodes() - 2;
}

void add_edge_pair(int v, int u, const string &seq) {
    if (edges_set.count({v, u})) {
        return;
    }
    string rev = rev_comp(seq);
    if (seq != rev) {
        int u_comp = comp_node(u);
        int v_comp = comp_node(v);
        edges.push_back({u_comp, v_comp, rev, 0LL});
        nodes[u_comp].out.push_back(n_edges() - 1);
        nodes[v_comp].in.push_back(n_edges() - 1);
        edges_set.insert({u_comp, v_comp});
    }
    edges.push_back({v, u, seq, 0LL});
    nodes[v].out.push_back(n_edges() - 1);
    nodes[u].in.push_back(n_edges() - 1);
    edges_set.insert({v, u});
}

void process_one_read(const string &read) {
    int prev = -1;
    for (int i = 0; i <= (int) read.size() - k; ++i) {
        string kmer = read.substr(i, k);
        int cur = get_node(kmer);
        if (prev != -1) {
            add_edge_pair(prev, cur, nodes[prev].kmer + kmer[kmer.size() - 1]);
        }
        //++nodes[cur].cnt;
        //++nodes[comp_node(cur)].cnt;
        prev = cur;
    }
    for (int i = 0; i <= (int) read.size() - k - 1; ++i) {
        string kmer = read.substr(i, k + 1);
        string rev = rev_comp(kmer);
        ++kmer_cnt[kmer];
        ++kmer_cnt[rev];
    }
}

// ----------------------
// | contracting graph  |
// ----------------------

void remove_redundant_node(int v) {
    int e_in = nodes[v].in[0];
    int e_out = nodes[v].out[0];
    assert(nodes[v].in.size() == 1 && nodes[v].out.size() == 1 && edges[e_in].v != v);
    edges[e_in].u = edges[e_out].u;
    edges[e_in].sum_kmer_cnt += edges[e_out].sum_kmer_cnt;
    edges[e_in].seq += edges[e_out].seq.substr(k);
    nodes[v] = node("");
    edges[e_out] = {-1, -1, "", 0};
    replace(nodes[edges[e_in].u].in.begin(), nodes[edges[e_in].u].in.end(), e_out, e_in);
}

void remove_redundant_nodes() {
    for (int i = 0; i < (int) n_nodes(); ++i) {
        if (nodes[i].out.size() == 1 && nodes[i].in.size() == 1 && nodes[i].in[0] != nodes[i].out[0]) {
            remove_redundant_node(i);
        }
    }
}

// ------------------
// | simplification |
// ------------------

// an edge E, which is the only incident edge for some vertex V, is considered to be a tip iff its length is 
// less than 3*K and its coverage is less than 0.5*max_cov, where max_cov is the maximum coverage among 
// all edges **incident to the E's other end**
// (reasoning: the edge with maximum coverage is likely to be from the "main path")
bool try_remove_tip(int v) {
    if (nodes[v].in.size() + nodes[v].out.size() != 1) {
        return false;
    }
    int e_id = nodes[v].in.size() ? nodes[v].in[0] : nodes[v].out[0];
    edge &e = edges[e_id];
    int u = e.v == v ? e.u : e.v;
    double max_cov = get_cov(e);
    for (int i : nodes[u].in) {
        max_cov = max(max_cov, get_cov(edges[i]));
    }
    for (int i : nodes[u].out) {
        max_cov = max(max_cov, get_cov(edges[i]));
    }
    if ((int) e.seq.size() < 3 * k && get_cov(e) < 0.75 * max_cov) { // check if e is a tip
        nodes[v] = node("");
        edges[e_id] = {-1, -1, "", 0};
        nodes[u].in.erase(remove(nodes[u].in.begin(), nodes[u].in.end(), e_id), nodes[u].in.end());
        nodes[u].out.erase(remove(nodes[u].out.begin(), nodes[u].out.end(), e_id), nodes[u].out.end());
        return true;
    }
    return false;
}

void remove_tips() {
    bool at_least_one_removed;
    do {
        at_least_one_removed = false;
        for (int v = 0; v < (int) n_nodes(); ++v) {
            if (try_remove_tip(v)) {
                at_least_one_removed = true;
            }
        }
    } while (at_least_one_removed);
}

// an edge is considered low quality if its length is less than 3*K and its coverage is less than a given threshold 
void remove_bad_qual_edges() {
    for (int i = 0; i < (int) n_edges(); ++i) {
#ifndef BAD_COVERAGE_THRESHOLD
#define BAD_COVERAGE_THRESHOLD 0
#endif
        if (edges[i].seq.size() && (int) edges[i].seq.size() < 3 * k && get_cov(edges[i]) < BAD_COVERAGE_THRESHOLD) {
            int v = edges[i].v;
            int u = edges[i].u;
            edges[i] = {-1, -1, "", 0};
            nodes[v].out.erase(remove(nodes[v].out.begin(), nodes[v].out.end(), i), nodes[v].out.end());
            nodes[u].in.erase(remove(nodes[u].in.begin(), nodes[u].in.end(), i), nodes[u].in.end());
        }
    }
}

// a bulge is a pair of edges between two nodes V and U, which are the only in edges for U, the only out edges
// for V, which satisfy the following quality criteria: their lengths are less than 3*k and
// <smaller coverage> is less than 0.75*<bigger coverage>
void remove_bulges() {
    for (int v = 0; v < (int) n_nodes(); ++v) {
        if (nodes[v].out.size() != 2) {
            continue;
        }
        int e_a = nodes[v].out[0];
        int e_b = nodes[v].out[1];
        if (edges[e_a].u != edges[e_b].u) {
            continue;
        }
        int u = edges[e_a].u;
        if (nodes[u].in.size() > 2) {
            continue;
        }
        if (get_cov(edges[e_a]) < get_cov(edges[e_b])) {
            swap(e_a, e_b);
        }
        if ((int) edges[e_a].seq.size() < 3 * k && (int) edges[e_b].seq.size() < 3 * k 
                && get_cov(edges[e_b]) < 0.75 * get_cov(edges[e_a])) {
            edges[e_b] = {-1, -1, "", 0};
            nodes[v].out.erase(remove(nodes[v].out.begin(), nodes[v].out.end(), e_b), nodes[v].out.end());
            nodes[u].in.erase(remove(nodes[u].in.begin(), nodes[u].in.end(), e_b), nodes[u].in.end());
        }
    }
}

// ----------
// | output |
// ----------

void write_fasta() {
    seqan::SeqFileOut fasta("edges.fasta");
    stringstream name;
    name << fixed << setprecision(6);
    for (int i = 0; i < (int) n_edges(); ++i) {
        edge &e = edges[i];
        if (!e.seq.size()) {
            continue;
        }
        name << "EDGE_" 
            << i 
            << "_length_" 
            << e.seq.size() - k
            << "_cov_"
            << get_cov(e);
        seqan::writeRecord(fasta, name.str(), e.seq);
        name.str("");
    }
}

void write_dot() {
    ofstream dot("graph.dot");
    dot << fixed << setprecision(6);
    dot << "digraph debrujin_graph {\n";
    for (edge &e : edges) {
        if (!e.seq.size()) {
            continue;
        }
        dot << e.v 
            << " -> " 
            << e.u
            << " [label=\"length = " 
            << e.seq.size() - k
            << "; cov = " 
            << get_cov(e)
            << "\"]\n";
    }
    dot << "}\n";
}

// --------
// | main |
// --------

int main(int argc, char **argv) {
    // get arguments
    assert(argc == 3);
    char* fn = argv[1];

    k = stoi(argv[2]);

    // read and process input
    seqan::SeqFileIn fin(fn);
    while (!seqan::atEnd(fin)) {
        string name, seq;
        seqan::readRecord(name, seq, fin);
        process_one_read(seq);
    }
    edges_set.clear();

    // calculate initial kmer counts on edges
    for (edge &e : edges) {
        e.sum_kmer_cnt = kmer_cnt[e.seq];
    }

    // contract graph
    remove_redundant_nodes();

#ifdef REMOVE_BAD_QUAL_EDGES_ENABLED
#ifdef BAD_COVERAGE_THRESHOLD
    // remove short and low coverage edges
    remove_bad_qual_edges();
    remove_redundant_nodes();
#endif
#endif

#ifdef REMOVE_TIPS_ENABLED
    // remove tips
    remove_tips();
    remove_redundant_nodes();
#endif

#ifdef REMOVE_BULGES_ENABLED
    // remove bulges
    remove_bulges();
    remove_redundant_nodes();
#endif

    // write to FASTA
    write_fasta();

    // write to DOT
    write_dot();
}
