// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "graph.h"
#include "pvector.h"
#include "util.h"

/*
GAP Benchmark Suite
File:   Benchmark
Author: Scott Beamer

Various helper functions to ease writing of kernels
*/

// Default type signatures for commonly used types
typedef int32_t NodeID;
typedef int32_t WeightT;
typedef NodeWeight<NodeID, WeightT> WNode;

typedef CSRGraph<NodeID> Graph;
typedef CSRGraph<NodeID, WNode> WGraph;

template <typename NodeID_, typename DestID_ = NodeID_, typename WeightT_ = NodeID_, bool invert = true>
class Reader {
    typedef EdgePair<NodeID_, DestID_> Edge;
    typedef pvector<Edge> EdgeList;
    std::string filename_;

   public:
    explicit Reader(std::string filename) : filename_(filename) {
    }

    std::string GetSuffix() {
        std::size_t suff_pos = filename_.rfind('.');
        if (suff_pos == std::string::npos) {
            std::cout << "Couldn't find suffix of " << filename_ << std::endl;
            std::exit(-1);
        }
        return filename_.substr(suff_pos);
    }

    CSRGraph<NodeID_, DestID_, invert> ReadSerializedGraph() {
        bool weighted = GetSuffix() == ".wsg";
        if (!std::is_same<NodeID_, SGID>::value) {
            std::cout << "serialized graphs only allowed for 32bit" << std::endl;
            std::exit(-5);
        }
        if (!weighted && !std::is_same<NodeID_, DestID_>::value) {
            std::cout << ".sg not allowed for weighted graphs" << std::endl;
            std::exit(-5);
        }
        if (weighted && std::is_same<NodeID_, DestID_>::value) {
            std::cout << ".wsg only allowed for weighted graphs" << std::endl;
            std::exit(-5);
        }
        if (weighted && !std::is_same<WeightT_, SGID>::value) {
            std::cout << ".wsg only allowed for int32_t weights" << std::endl;
            std::exit(-5);
        }
        std::ifstream file(filename_);
        if (!file.is_open()) {
            std::cout << "Couldn't open file " << filename_ << std::endl;
            std::exit(-6);
        }
        Timer t;
        t.Start();
        bool directed;
        SGOffset num_nodes, num_edges;
        DestID_ **index = nullptr, **inv_index = nullptr;
        DestID_ *neighs = nullptr, *inv_neighs = nullptr;
        file.read(reinterpret_cast<char *>(&directed), sizeof(bool));
        file.read(reinterpret_cast<char *>(&num_edges), sizeof(SGOffset));
        file.read(reinterpret_cast<char *>(&num_nodes), sizeof(SGOffset));
        pvector<SGOffset> offsets(num_nodes + 1);
        neighs = new DestID_[num_edges];
        std::streamsize num_index_bytes = (num_nodes + 1) * sizeof(SGOffset);
        std::streamsize num_neigh_bytes = num_edges * sizeof(DestID_);
        file.read(reinterpret_cast<char *>(offsets.data()), num_index_bytes);
        file.read(reinterpret_cast<char *>(neighs), num_neigh_bytes);
        index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        if (directed && invert) {
            inv_neighs = new DestID_[num_edges];
            file.read(reinterpret_cast<char *>(offsets.data()), num_index_bytes);
            file.read(reinterpret_cast<char *>(inv_neighs), num_neigh_bytes);
            inv_index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, inv_neighs);
        }
        file.close();
        t.Stop();
        PrintTime("Read Time", t.Seconds());
        if (directed)
            return CSRGraph<NodeID_, DestID_, invert>(num_nodes, index, neighs, inv_index, inv_neighs);
        else
            return CSRGraph<NodeID_, DestID_, invert>(num_nodes, index, neighs);
    }
};

CSRGraph<NodeID, NodeID, true> MakeUnweightedGraph(std::string filename) {
    Reader<NodeID, NodeID, WeightT, true> r(filename);
    return r.ReadSerializedGraph();
}

CSRGraph<NodeID, WNode, true> MakeWeightedGraph(std::string filename) {
    Reader<NodeID, WNode, WeightT, true> r(filename);
    return r.ReadSerializedGraph();
}

template <typename NodeID_, typename rng_t_, typename uNodeID_ = typename std::make_unsigned<NodeID_>::type>
class UniDist {
   public:
    UniDist(NodeID_ max_value, rng_t_ &rng) : rng_(rng) {
        no_mod_ = rng_.max() == static_cast<uNodeID_>(max_value);
        mod_ = max_value + 1;
        uNodeID_ remainder_sub_1 = rng_.max() % mod_;
        if (remainder_sub_1 == mod_ - 1)
            cutoff_ = 0;
        else
            cutoff_ = rng_.max() - remainder_sub_1;
    }

    NodeID_ operator()() {
        uNodeID_ rand_num = rng_();
        if (no_mod_)
            return rand_num;
        if (cutoff_ != 0) {
            while (rand_num >= cutoff_)
                rand_num = rng_();
        }
        return rand_num % mod_;
    }

   private:
    rng_t_ &rng_;
    bool no_mod_;
    uNodeID_ mod_;
    uNodeID_ cutoff_;
};

// Used to pick random non-zero degree starting points for search algorithms
template <typename GraphT_>
class SourcePicker {
   public:
    explicit SourcePicker(const GraphT_ &g)
        : given_source_(-1), non_changing_(true), rng_(kRandSeed), udist_(g.num_nodes() - 1, rng_), g_(g) {
    }

    NodeID PickNext() {
        if (given_source_ != -1)
            return given_source_;
        if (non_changing_ && last_ != -1)
            return last_;

        NodeID source;
        do {
            source = udist_();
        } while (g_.out_degree(source) == 0);
        last_ = source;

        return source;
    }

   private:
    NodeID given_source_;
    bool non_changing_;
    std::mt19937_64 rng_;
    UniDist<NodeID, std::mt19937_64> udist_;
    const GraphT_ &g_;
    NodeID last_ = -1;
};

#endif  // BENCHMARK_H_
