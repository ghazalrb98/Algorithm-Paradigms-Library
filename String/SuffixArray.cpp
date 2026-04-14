/*
 *  Suffix Array (Prefix-Doubling Implementation)
 *
* Authors:
*  - Mohsen Dehbashi
*  - Ghazal Rabie
 *
 *  Description:
 *  ------------
 *  This file implements a Suffix Array data structure using the
 *  prefix-doubling algorithm.
 *
 *  A suffix array stores the starting indices of all suffixes
 *  of a string, sorted in lexicographic order.
 *
 *  Example:
 *      s = "popup"
 *      suffixes:
 *          0 -> "popup"
 *          1 -> "opup"
 *          2 -> "pup"
 *          3 -> "up"
 *          4 -> "p"
 *
 *      sorted suffixes:
 *          "opup"  -> 1
 *          "p"     -> 4
 *          "popup" -> 0
 *          "pup"   -> 2
 *          "up"    -> 3
 *
 *      suffix array = [1, 4, 0, 2, 3]
 *
 *
 *  Supported Operations:
 *  ---------------------
 *      - Constructor(string s): builds suffix array
 *      - getSuffix(i): returns index of i-th smallest suffix
 *      - getArray(): returns the full suffix array
 *
 *
 *  Complexity:
 *  -----------
 *      Construction:  O(n log^2 n)
 *          - log n rounds (doubling k)
 *          - each round uses std::sort → O(n log n)
 *
 *      Query:
 *          getSuffix(i): O(1)
 *
 *      Memory:
 *          O(n)
 *
 *
 *  Notes:
 *  ------
 *  - This implementation uses comparison-based sorting for clarity.
 *  - Early stopping is used when all suffix ranks become unique.
 *
 */

#include <vector>
#include <string>
#include <algorithm>

class SuffixArray {
private:
    std::string s;
    std::vector<int> sa;      // sa[i] = starting index of i-th smallest suffix
    std::vector<int> rank;    // rank[i] = equivalence class of suffix starting at i

public:
    // Constructor: builds the suffix array
    explicit SuffixArray(const std::string& str): s(str) {
        build();
    }

    // Returns the starting index of the i-th smallest suffix
    int getSuffix(const int i) const {
        return sa[i];
    }

    // Returns the full suffix array
    const std::vector<int>& getArray() {
        return sa;
    }

private:
    // Builds the suffix array using prefix doubling
    void build() {
        const int n = static_cast<int>(s.size());
        sa.resize(n);
        rank.resize(n);

        for (int i = 0; i < n; i++) {
            sa[i] = i;
            rank[i] = static_cast<unsigned char>(s[i]);
        }

        std::vector<int> newRank(n);

        // k = current prefix length (doubles each round)
        for (int k = 1; k < n; k <<= 1) {
            // Comparator: compares suffixes by (rank[i], rank[i + k])
            auto cmp = [&](const int i, const int j) {
                if (rank[i] != rank[j])
                    return rank[i] < rank[j];

                const int ri = (i + k < n ? rank[i + k] : -1);
                const int rj = (j + k < n ? rank[j + k] : -1);
                return  ri < rj;
            };

            sort(sa.begin(), sa.end(), cmp);

            // Assign new ranks after sorting
            newRank[sa[0]] = 0;
            for (int i = 1; i < n; i++) {
                newRank[sa[i]] = newRank[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
            }
            rank = newRank;

            // Early stop: all suffixes have unique ranks
            if (rank[sa[n - 1]] == n - 1) break;
        }
    }
};
