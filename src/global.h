#pragma once

#include <vector>
#include <iostream>
#include <unordered_map>
#include <list>
#include <algorithm>
#include <cstdint>

#include "kmers.h"
#include "khash.h"
#include "khash_utils.h"
#include "simplitigs.h"

/// Provide possibility to access reverse complements as if they were in the field.
#define accessPrefix(field, index) (((field).size() > (index)) ? simplitig_first(kmer_t(0), (field)[(index)], k) : \
        ReverseComplement(simplitig_last(kmer_t(0), (field)[(index) - (field).size()], k), k))

/// Provide possibility to access reverse complements as if they were in the field.
#define accessSuffix(field, index) (((field).size() > (index)) ? simplitig_last(kmer_t(0), (field)[(index)], k) : \
        ReverseComplement(simplitig_first(kmer_t(0), (field)[(index) - (field).size()], k), k))

/// Provide possibility to access reverse complements as if they were in the field.
#define accessSimplitig(field, index) (((field).size() > (index)) ? (field)[(index)] : \
        simplitig_reverse_complement((field)[(index) - (field).size()]))

/// Provide possibility to know first and last for reverse complements without storing them.
#define accessFirstLast(a, b, index, n) (((n) > (index)) ? (a)[(index)] : \
        (((b)[(index) - (n)]) + (n)) % (2*(n)))


/// Determines which fraction of k-mers store its prefixes at one time.
int MEMORY_REDUCTION_FACTOR = 16;
/// Determines the number of prefix bits based on which the k-mers are presorted.
constexpr int SORT_FIRST_BITS_DEFAULT = 8;

typedef std::pair<std::vector<size_t>, std::vector<unsigned char>> overlapPath;

/// Rearrange the k-mers so that k-mers next to each other in sorted order appear close so that they are in the same bucket.
template <typename kmer_t>
void PartialPreSort(std::vector<kmer_t> &vals, int k) {
    int SORT_FIRST_BITS = std::min(2 * k, SORT_FIRST_BITS_DEFAULT);
    uint64_t DIFFERENT_PREFIXES_COUNT = 1ULL << SORT_FIRST_BITS;
    kmer_t PREFIX_MASK = DIFFERENT_PREFIXES_COUNT - kmer_t(1);
    std::vector<size_t> counts(DIFFERENT_PREFIXES_COUNT, 0);
    int shift = (2 * k) - SORT_FIRST_BITS;
    kmer_t mask = PREFIX_MASK << shift;
    for (auto &&kMer : vals) counts[(uint64_t)((kMer & mask) >> shift)]++;
    std::vector<std::vector<kmer_t>> distributed(DIFFERENT_PREFIXES_COUNT);
    for (uint64_t i = 0; i < DIFFERENT_PREFIXES_COUNT; ++i) distributed[i] = std::vector<kmer_t> (counts[i]);
    for (uint64_t i = 0; i < DIFFERENT_PREFIXES_COUNT; ++i) counts[i] = 0;
    for (auto &&kMer : vals) {
        uint64_t index = (kMer & mask) >> shift;
        distributed[index][counts[index]++] = kMer;
    }
    size_t index = 0;
    for (auto && bucket : distributed) {
        for (auto && kMer : bucket) {
            vals[index++] = kMer;
        }
    }
}

/// Greedily find the approximate Hamiltonian path with longest overlaps.
/// k is the size of one k-mer and n is the number of distinct k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected to contain only one k-mer from a complement pair.
/// Moreover, if so, the resulting Hamiltonian path contains two superstrings which are reverse complements of one another.
/// If lower_bound is set to true, return a shortest cycle cover instead.
template <typename kmer_t, typename kh_wrapper_t>
overlapPath OverlapHamiltonianPath (kh_wrapper_t wrapper, kmer_t kmerType, std::vector<simplitig_t> &simplitigs, int k, bool complements,
                                    bool lower_bound = false) {
    size_t n = simplitigs.size();
    size_t kMersCount = n * (1 + complements);
    size_t batchSize = kMersCount / MEMORY_REDUCTION_FACTOR + 1;
    std::vector<size_t> edgeFrom(kMersCount, -1);
    std::vector<unsigned char> overlaps(kMersCount, -1);
    std::vector<bool> suffixForbidden(kMersCount, false);
    std::vector<bool> prefixForbidden(kMersCount, false);
    // For reverse complements, compute first from last and vice versa.
    auto first = new size_t[n];
    auto last = new size_t[n];
    // Index next relative to the batch.
    auto next = new size_t[batchSize];
    for (size_t i = 0; i < n; ++i) {
        first[i] = last[i] = i;
    }
    auto *prefixes = wrapper.kh_init_map();
    wrapper.kh_resize_map(prefixes, (kMersCount / MEMORY_REDUCTION_FACTOR + 1 ) * 100 / 77 );
    for (int d = k - 1; d >= 0; --d) {
        // In order to reduce memory requirements, the prefixes are not processed at once, but in batches.
        // As a cost, this slows down the algorithm.
        for (int part = 0; part < MEMORY_REDUCTION_FACTOR; part++) {
            wrapper.kh_clear_map(prefixes);
            for (size_t i = 0; i < batchSize; ++i) {
                next[i] = (size_t)-1;
            }
            size_t to = std::min(kMersCount, (part + 1) * batchSize);
            size_t from = part * batchSize;
            for (size_t i = from; i < to; ++i)
                if (!prefixForbidden[i]) {
                    next[i - from] = -1;
                    kmer_t prefix = BitPrefix(accessPrefix(simplitigs,i), k, d);
                    auto prefix_key = wrapper.kh_get_from_map(prefixes, prefix);
                    if (prefix_key != kh_end(prefixes)) {
                        next[i - from] = kh_val(prefixes, prefix_key);
                    } else {
                        int ret;
                        prefix_key = wrapper.kh_put_to_map(prefixes, prefix, &ret);
                    }
                    kh_value(prefixes, prefix_key) = i;
                }
            for (size_t i = 0; i < kMersCount; ++i)
                if (!suffixForbidden[i]) {
                    kmer_t suffix = BitSuffix(accessSuffix(simplitigs, i), d);
                    auto suffix_key = wrapper.kh_get_from_map(prefixes, suffix);
                    if (suffix_key == kh_end(prefixes)) continue;
                    size_t previous, j;
                    previous = j = kh_val(prefixes, suffix_key);
                    while (j != size_t(-1) && \
                            // k-mers are complementary
                           ((i + n) % (2 * n) == j \
                           // forms a cycle
                           || (!lower_bound && accessFirstLast(first, last, i, n) == j) \
                           // k-mer is already used
                           || prefixForbidden[j])) {
                        size_t new_j = next[j - from];
                        // If the k-mer is forbidden, remove it to keep the complexity linear.
                        // This is not done with the first k-mer but that is not a problem.
                        if (prefixForbidden[j]) next[previous - from] = new_j;
                        else previous = j;
                        j = new_j;
                    }
                    if (j == size_t(-1)) {
                        continue;
                    }
                    std::vector<std::pair<size_t, size_t>> new_edges({{i, j}});
                    // Add also the edge between complementary k-mers in the opposite direction.
                    if (complements) new_edges.emplace_back((j + n) % kMersCount, (i + n) % kMersCount);
                    for (auto [x, y]: new_edges) {
                        edgeFrom[x] = y;
                        overlaps[x] = d;
                        prefixForbidden[y] = true;
                        auto lastY =  accessFirstLast(last, first, y, n);
                        auto firstX = accessFirstLast(first, last, x, n);
                        if (lastY < n) first[lastY] = firstX;
                        if (firstX < n) last[firstX] = lastY;
                        suffixForbidden[x] = true;
                    }
                    next[previous - from] = next[j - from];
                }
        }
    }

    wrapper.kh_destroy_map(prefixes);
    delete[](next);
    delete[](first);
    delete[](last);
    return {edgeFrom, overlaps};
}


void PrintSimplitigStart(simplitig_t simplitig, std::ostream& of, std::ostream *maskf, int k) {
    int count = kmers_in_simplitig(simplitig, k);
    for (size_t i = 0; i < count; ++i) {
        char nucleotide = simplitig_at_index(simplitig, i);
        of << nucleotide;
        if (maskf != nullptr) (*maskf) << nucleotide;
    }
}


/// Construct the superstring and its mask from the given overlapPath path in the overlap graph.
/// If reverse complements are considered and the overlapPath path contains two paths which are reverse complements of one another,
/// return only one of them.
template <typename kmer_t, typename kh_wrapper_t >
void SuperstringFromPath(kh_wrapper_t wrapper, kmer_t kmerType, const overlapPath &hamiltonianPath, const std::vector<simplitig_t> &simplitigs, std::ostream& of, std::ostream *maskf, const int k, const bool complements) {
    size_t kMersCount = simplitigs.size() * (1 + complements);
    auto edgeFrom = hamiltonianPath.first;
    auto overlaps = hamiltonianPath.second;

    // Find the vertex in the overlap graph with in-degree 0.
    std::vector<bool> isStart(kMersCount, true);
    for (auto edge : edgeFrom) {
        if (edge != size_t(-1)) isStart[edge] = false;
    }
    size_t start = 0;
    for (; start < kMersCount && !isStart[start]; ++start);

    auto *kMersDict = wrapper.kh_init_set();

    if (maskf != nullptr) {
        for (auto simplitig : simplitigs) {
            auto simplitigKMers = kmers_in_simplitig(simplitig, k);
            for (size_t i = 0; i < simplitigKMers; ++i) {
                auto kMer = kmer_at_simplitig_index(kmerType, simplitig, k, i);
                int ret;
                wrapper.kh_put_to_set(kMersDict, kMer, &ret);
            }
        }
    }


    simplitig_t simplitig = accessSimplitig(simplitigs, start);
    kmer_t last = simplitig_last(kmerType, simplitig, k);

    PrintSimplitigStart(simplitig, of, maskf, k);

    // Move from the first k-mer to the last which has no successor.
    while(edgeFrom[start] != size_t(-1)) {
        int overlapLength = overlaps[start];
        simplitig = accessSimplitig(simplitigs, edgeFrom[start]);
        std::string unmaskedNucleotides = NumberToKMer(BitPrefix(last, k-1, k-1-overlapLength), k-1-overlapLength);
        std::transform(unmaskedNucleotides.begin(), unmaskedNucleotides.end(), unmaskedNucleotides.begin(), tolower);
        of << unmaskedNucleotides;
        // If provided, output also mask maximizing ones.
        if (maskf != nullptr) for (int j = 0; j < k - overlapLength - 1; ++j) {
            kmer_t current = simplitig_first(kmerType, simplitig, k);
            last <<= 2;
            last |= AtIndex(current, k, overlapLength + j);
            last = BitSuffix(last, k);
            if (containsKMer(kMersDict, wrapper, last, k, complements)) {
                (*maskf) << Masked(NucleotideAtIndex(last, k, 0), true);
            } else {
                (*maskf) << Masked(NucleotideAtIndex(last, k, 0), false);
            }
        }
        last = simplitig_last(kmerType, simplitig, k);
        PrintSimplitigStart(simplitig, of, maskf, k);
        start = edgeFrom[start];
    }

    // Print the trailing k-1 characters.
    last = BitSuffix(last, k-1);
    std::string unmaskedNucleotides = NumberToKMer(last, k-1);
    std::transform(unmaskedNucleotides.begin(), unmaskedNucleotides.end(), unmaskedNucleotides.begin(), tolower);
    of << unmaskedNucleotides;
    if (maskf != nullptr) {
        (*maskf) << unmaskedNucleotides << std::endl;
    }
}

/// Get the approximated shortest superstring of the given k-mers using the global greedy algorithm.
///
/// This runs in O(n k), where n is the number of k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
/// Warning: this will destroy kMers.
template <typename kmer_t, typename kh_wrapper_t>
void Global(kh_wrapper_t wrapper, kmer_t kmerType, std::vector<simplitig_t> &simplitigs, std::ostream& of, std::ostream *maskf, int k, bool complements) {
    if (simplitigs.empty()) {
        throw std::invalid_argument("input cannot be empty");
    }
    auto hamiltonianPath = OverlapHamiltonianPath(wrapper, kmerType, simplitigs, k, complements);
    SuperstringFromPath(wrapper, kmerType, hamiltonianPath, simplitigs, of, maskf, k, complements);
}

// Undefine the access macro, so it does not interfere with other files.
#undef access
