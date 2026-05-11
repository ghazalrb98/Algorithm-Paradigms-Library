/*
 * Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 * KMP (Knuth-Morris-Pratt) substring search implementation
 *
 * Provides functionality to find all occurrences of a pattern in a given text in O(|text| + |pattern|) time
 */

#include <vector>
#include <string>

class KMP {
    /*
     * Computes the LPS (Longest Prefix Suffix) array for the pattern
     *
     * lps[i] = length of the longest proper prefix of pattern[0...i] which is also a suffix of pattern[0...i]
     *
     * Time complexity: O(|pattern|)
     * Space complexity: O(|pattern|)
     */
    static std::vector<size_t> computeLPS(const std::string& pattern) {
        const size_t m = pattern.size();
        std::vector<size_t> lps(m, 0);

        size_t i = 0; // length of previous longest prefix suffix
        size_t j = 1;

        while (j < m) {
            if (pattern[i] == pattern[j]) {
                lps[j] = i + 1;
                i++;
                j++;
            } else if (i > 0) {
                i = lps[i - 1];
            } else {
                lps[j] = 0;
                j++;
            }
        }

        return lps;
    }

public:
    /*
     * Finds all occurrences of pattern in text
     *
     * @param pattern The pattern to search for
     * @param text The text to search in
     * @return A vector containing starting indices of all matches
     *
     * Time complexity: O(|text| + |pattern|)
     * Space complexity: O(|text| + |pattern|)
     */
    static std::vector<size_t> findAll(const std::string& pattern, const std::string& text) {
        std::vector<size_t> result;

        const size_t n = text.size();
        const size_t m = pattern.size();

        if (n == 0 || m == 0 || m > n)
            return result;

        const std::vector<size_t> lps = computeLPS(pattern);

        size_t i = 0; // index for text
        size_t j = 0; // index for pattern

        while (i < n) {
            if (text[i] == pattern[j]) {
                i++;
                j++;

                if (j == m) {
                    result.push_back(i - j); // match found
                    j = lps[j - 1];
                }
            } else if (j > 0) {
                j = lps[j - 1];
            } else {
                i++;
            }
        }

        return result;
    }
};