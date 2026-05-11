/*
 * Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 * Polynomial rolling hash utilities for O(1) substring hash queries
 *
 * hashes can be retrieved in constant time: hash(l, r) = s[l]*p^0 + s[l+1]*p^1 + ... + s[r]*p^(r-l)  (mod)
 */

#include <string>
#include "NumberTheory.cpp"

namespace StringHashing {

    template <typename T = long long>
    struct HashConfig {
        T p;
        T mod;
    };

    /*
     * Single polynomial rolling hash over a string
     *
     */
    template <typename T = long long>
    class SubstringHash {
        std::string str{};
        T p;
        NumberTheory::ModularArithmetic<T> ma;
        std::vector<T> prefixHash{};
        std::vector<T> pInverse{};

    public:
        /*
         * Builds prefix hash array and precomputes modular inverses of powers of p
         *
         * @param str   The input string
         * @param cfg   Hash config: {p (base), mod (modulus)}
         * @param base  Character offset
         *
         * Time complexity: O(n)
         * Space complexity: O(n)
         */
        explicit SubstringHash(const std::string &str, HashConfig<T> cfg, T base) : str(str), p(cfg.p), ma(cfg.mod) {
            size_t n = str.size();

            prefixHash.reserve(n);
            prefixHash.push_back(str[0] - base);
            T pn = p;
            for (size_t i = 1; i < n; i++) {
                prefixHash.push_back(ma.add(prefixHash.back(), ma.mul(str[i] - base, pn)));
                pn = ma.mul(pn, p);
            }

            pInverse.reserve(n);
            pInverse.push_back(1);
            T pi = ma.inverse(p);
            for (size_t i = 1; i < n; i++) {
                pInverse.push_back(ma.mul(pInverse.back(), pi));
            }
        }

        /*
         * Returns the hash of the substring s[l..r] (inclusive)
         *
         * @param l  Left index (inclusive)
         * @param r  Right index (inclusive)
         * @return   Hash value of s[l..r]
         *
         * Time complexity: O(1)
         * Space complexity: O(1)
         */
        T getHash(size_t l, size_t r) const {
            T diff = ma.sub(prefixHash[r], l ? prefixHash[l - 1] : 0);
            return ma.mul(diff, pInverse[l]);
        }
    };

    /*
     * Double polynomial rolling hash over a string
     *
     * Runs two independent SubstringHash instances with different (p, mod) pairs,
     * returning a pair of hashes per query to significantly reduce collision probability.
     */
    template <typename T = long long>
    class DoubleSubstringHash {
        SubstringHash<T> hash1;
        SubstringHash<T> hash2;

    public:
        /*
         * Builds two independent hash structures over the same string
         *
         * @param str   The input string
         * @param cfg1  Hash config for the first hash: {p, mod}
         * @param cfg2  Hash config for the second hash: {p, mod}
         * @param base  Character offset
         *
         * Time complexity: O(n)
         * Space complexity: O(n)
         */
        explicit DoubleSubstringHash(
            const std::string& str, HashConfig<T> cfg1, HashConfig<T> cfg2, T base)
            : hash1(str, cfg1, base),
              hash2(str, cfg2, base) {}

        /*
         * Returns a pair of hashes for the substring s[l..r] (inclusive)
         *
         * @param l  Left index (inclusive)
         * @param r  Right index (inclusive)
         * @return   Pair {hash1(l, r), hash2(l, r)}
         *
         * Time complexity: O(1)
         * Space complexity: O(1)
         */
        std::pair<T, T> getHash(size_t l, size_t r) const {
            return std::make_pair(hash1.getHash(l, r), hash2.getHash(l, r));
        }
    };
}
