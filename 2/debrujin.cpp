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
#include "utils.hpp"

using namespace std;

int k;
vector<node> nodes;
vector<edge> edges;
unordered_map<string, int> kmer2node;
set<pair<int, int>> edges_set;

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
    return ((double) e.sum_kmer_cnt / (e.seq.size() - k + 1));
}

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
        ++nodes[cur].cnt;
        ++nodes[comp_node(cur)].cnt;
        prev = cur;
    }
}

void remove_redundant_node(int v) {
    int e_in = nodes[v].in[0];
    int e_out = nodes[v].out[0];
    assert(nodes[v].in.size() == 1 && nodes[v].out.size() == 1 && edges[e_in].v != v);
    edges[e_in].u = edges[e_out].u;
    edges[e_in].sum_kmer_cnt += edges[e_out].sum_kmer_cnt - nodes[v].cnt;
    edges[e_in].seq += edges[e_out].seq.substr(k);
    nodes[v] = node("");
    edges_set.erase({edges[e_out].v, edges[e_out].u});
    edges[e_out] = {-1, -1, "", 0};
    replace(nodes[edges[e_in].u].in.begin(), nodes[edges[e_in].u].in.end(), e_out, e_in);
}

// an edge E is considered to be a tip iff its length is less than 2*k and its coverage is less than 0.5*max_cov,
// where max_cov is the maximum coverage among all edges **incident to the E's other end**.
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
    if ((int) e.seq.size() < 2 * k && get_cov(e) < 0.5 * max_cov) { // check if e is a tip
        nodes[v] = node("");
        edges_set.erase({e.v, e.u});
        edges[e_id] = {-1, -1, "", 0};
        nodes[u].in.erase(remove(nodes[u].in.begin(), nodes[u].in.end(), e_id), nodes[u].in.end());
        nodes[u].out.erase(remove(nodes[u].out.begin(), nodes[u].out.end(), e_id), nodes[u].out.end());
    }
    return false;
}

void remove_redundant_nodes() {
    for (int i = 0; i < (int) n_nodes(); ++i) {
        if (nodes[i].out.size() == 1 && nodes[i].in.size() == 1 && nodes[i].in[0] != nodes[i].out[0]) {
            remove_redundant_node(i);
        }
    }
}

int main(int argc, char **argv) {
    // get arguments
    assert(argc == 3);
    char* fn = argv[1];

    k = stoi(argv[2]);

    bool REMOVE_TIPS_ENABLED = true;

    // read and process input
    seqan::SeqFileIn fin(fn);
    while (!seqan::atEnd(fin)) {
        string name, seq;
        seqan::readRecord(name, seq, fin);
        process_one_read(seq);
    }

    // calculate initial kmer counts on edges
    for (edge &e : edges) {
        e.sum_kmer_cnt = nodes[e.v].cnt + nodes[e.u].cnt;
    }

    // remove redundant nodes
    remove_redundant_nodes();

    // remove tips
    if (REMOVE_TIPS_ENABLED) {
        bool at_least_one_removed;
        do {
            at_least_one_removed = false;
            for (int v = 0; v < (int) n_nodes(); ++v) {
                if (try_remove_tip(v)) {
                    at_least_one_removed = true;
                }
            }
        } while (at_least_one_removed);
        remove_redundant_nodes();
    }

    // write to FASTA
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

    // write to DOT
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

    // debug output
    /*cout << n_nodes() << '\n';
    for (node &v : nodes) {
        cout << v.kmer << ' ' << v.in.size() << ' ' << v.out.size() << '\n';
    }
    cout << '\n';
    cout << n_edges() << '\n';
    for (edge &e : edges) {
        cout << e.v << ' ' << e.u << ' ' << e.seq << '\n';
    }*/
}
