/*
 * Suffix Array Implementation
 *
 * Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 * A suffix array stores the starting indices of all suffixes of a string, sorted in lexicographic order
 *
 * Example:
 *  s = "popup"
 *      suffixes:
 *          0 -> "popup"
 *          1 -> "opup"
 *          2 -> "pup"
 *          3 -> "up"
 *          4 -> "p"
 *      sorted suffixes:
 *          "opup"  -> 1
 *          "p"     -> 4
 *          "popup" -> 0
 *          "pup"   -> 2
 *          "up"    -> 3
 *      suffix array = [1, 4, 0, 2, 3]
 */

#include <vector>
#include <string>
#include <algorithm>

class SuffixArray {
    std::string s;
    std::vector<size_t> sa;      // sa[i] = starting index of i-th smallest suffix
    std::vector<size_t> rank;    // rank[i] = equivalence class of suffix starting at i

public:
    /**
     * Builds the suffix array using prefix doubling trick
     * Time complexity: O(n log^2 n)
     * Space complexity: O(n)
     */
    explicit SuffixArray(const std::string& str): s(str) {
        build();
    }

    /**
     * Returns the starting index of the i-th smallest suffix
     * O(1)
     */
    size_t kthSmallestSuffixStartIndex(const size_t k) const {
        return sa[k];
    }
private:
    void build() {
        const size_t n = s.size();
        sa.resize(n);
        rank.resize(n);

        for (size_t i = 0; i < n; i++) {
            sa[i] = i;
            rank[i] = static_cast<unsigned char>(s[i]);
        }

        std::vector<size_t> newRank(n);

        // k = current prefix length (doubles each round)
        for (size_t k = 1; k < n; k <<= 1) {
            // Comparator: compares suffixes by (rank[i], rank[i + k])
            auto cmp = [&](const size_t i, const size_t j) {
                if (rank[i] != rank[j])
                    return rank[i] < rank[j];

                const size_t ri = i + k < n ? rank[i + k] : 0;
                const size_t rj = j + k < n ? rank[j + k] : 0;
                return  ri < rj;
            };

            sort(sa.begin(), sa.end(), cmp);

            // Assign new ranks after sorting
            newRank[sa[0]] = 0;
            for (size_t i = 1; i < n; i++) {
                newRank[sa[i]] = newRank[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
            }
            rank = newRank;

            // Early stop: all suffixes have unique ranks
            if (rank[sa[n - 1]] == n - 1) break;
        }
    }
};