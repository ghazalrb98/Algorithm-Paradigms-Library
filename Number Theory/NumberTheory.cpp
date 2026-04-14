/*
* Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei
 *
 *  This file provides reusable number theory utilities:
 *   - Extended Euclidean algorithm
 *   - Greatest common divisor (gcd)
 *   - Least common multiple (lcm)
 *   - Modular normalization
 *   - Modular inverse
 *   - Modular arithmetic wrapper
 *   - Chinese Remainder Theorem solver for arbitrary moduli
 *
 * Notes:
 *   - The modular arithmetic class supports addition, subtraction,
 *     multiplication, and division modulo n.
 *   - Division is performed using the modular inverse, so it only
 *     succeeds when the inverse exists.
 *   - The CRT solver supports arbitrary moduli, not only relatively
 *     prime moduli.
 *
 * Complexity summary:
 *   - normalize_mod: O(1)
 *   - add/sub/mul: O(1)
 *   - inverse_mod: O(log(min(a, mod))) via Euclidean algorithm
 *   - CRT merge: O(log(min(n, m)))
 */

#include <vector>

namespace NumberTheory {
    /*
     * Extended Euclidean Algorithm
     *
     * Computes integers x and y such that:
     *      a*x + b*y = gcd(a, b)
     *
     * Parameters:
     *  - a, b: input integers
     *  - x, y: output coefficients
     *
     * Returns:
     *  - gcd(a, b)
     *
     * Time complexity:
     *  - O(log(min(a, b)))
     */
    template <typename T>
    T egcd(T a, T b, T& x, T& y) {
        if (b == 0) {
            x = 1;
            y = 0;
            return a;
        }

        T x1, y1;
        T g = egcd(b, a % b, x1, y1);
        x = y1;
        y = x1 - (a / b) * y1;
        return g;
    }

    template <typename T>
    T gcd(T a, T b) {
        T x, y;
        return egcd(a, b, x, y);
    }

    template <typename T>
    T lcm(T a, T b) {
        return (a * b) / gcd(a, b);
    }

    /*
     * Modular normalization
     * Returns x reduced into the interval [0, mod).
     */
    template <typename T>
    T normalize_mod(T x, T mod) {
        x %= mod;
        if (x < 0) x += mod;
        return x;
    }

    /*
     * Modular inverse
     *
     * Returns y such that:
     *      a * y ≡ 1 (mod mod)
     *
     * If gcd(a, mod) != 1, the inverse does not exist
     * and the function returns -1.
     *
     * Time complexity:
     *  - O(log(min(a, mod)))
     */
    template <typename T>
    T inverse_mod(T a, T mod) {
        a = normalize_mod(a, mod);

        T x, y;
        T g = egcd(a, mod, x, y);

        if (g != 1) return -1;
        return normalize_mod(x, mod);
    }

    /*
     * ModularArithmetic
     *
     * A helper class for modular arithmetic under a fixed modulus.
     *
     * Supported operations:
     *  - normalize
     *  - add
     *  - sub
     *  - mul
     *  - div  (via modular inverse)
     *  - inverse
     *
     * Complexity:
     *  - normalize/add/sub/mul: O(1)
     *  - inverse/div: O(log mod)
     */
    template <typename T>
    class ModularArithmetic {
        T mod;

    public:
        explicit ModularArithmetic(T n) : mod(n) {}

        T normalize(T x) const {
            return normalize_mod(x, mod);
        }

        T add(T a, T b) const {
            return normalize(normalize(a) + normalize(b));
        }

        T sub(T a, T b) const {
            return normalize(normalize(a) - normalize(b));
        }

        T mul(T a, T b) const {
            return normalize(normalize(a) * normalize(b));
        }

        /*
         * Computes modular inverse of a.
         *
         * Returns:
         *  - a^{-1} mod mod, if it exists
         *  - -1 otherwise
         *
         * Time complexity:
         *  - O(log mod)
        */
        T inverse(T a) const {
            return inverse_mod(a, mod);
        }

        /*
         * Computes (a / b) mod mod.
         *
         * This is implemented as:
         *      a * b^{-1} mod mod
         *
         * Returns:
         *  - the result if inverse(b) exists
         *  - -1 otherwise
         *
         * Time complexity:
         *  - O(log mod)
         */
        T div(T a, T b) const {
            T inv = inverse(b);
            if (inv == -1) return -1;
            return mul(a, inv);
        }
    };

    /*
     * ChineseRemainder
     *
     * Solves systems of congruences using the Chinese Remainder Theorem.
     *
     * This implementation supports arbitrary moduli.
     * Therefore, it also covers the relatively prime case.
     *
     * Example:
     *      x ≡ a (mod n)
     *      x ≡ b (mod m)
     *
     * If a solution exists, the result is unique modulo lcm(n, m).
     */
    template <typename T>
    class ChineseRemainder {
    public:
        struct Congruence
        {
            T remainder;
            T modulus;
        };

        struct CRTResult
        {
            bool hasSolution;
            Congruence value;
        };

        /*
         * Merges two congruences:
         *      x ≡ a (mod n)
         *      x ≡ b (mod m)
         *
         * If a solution exists, returns:
         *      x ≡ r (mod lcm(n, m))
         *
         * Otherwise returns hasSolution = false.
         *
         * Time complexity:
         *  - O(log(min(n, m)))
         */
        static CRTResult merge(const Congruence& c1, const Congruence& c2) {
            T a = normalize_mod(c1.remainder, c1.modulus);
            T n = c1.modulus;

            T b = normalize_mod(c2.remainder, c2.modulus);
            T m = c2.modulus;

            T g = gcd(n, m);

            if ((b - a) % g != 0) {
                return {false, {0, 0}};
            }

            T n1 = n / g;
            T m1 = m / g;
            T diff = (b - a) / g;

            T inv = inverse_mod(n1, m1);
            if (inv == -1) {
                return {false, {0, 0}};
            }

            T t = normalize_mod(diff * inv, m1);

            T modulus = n * m1;
            T remainder = normalize_mod(a + n * t, modulus);

            return {true, {remainder, modulus}};
        }

        /*
         * Solves a list of congruences:
         *      x ≡ a_i (mod n_i)
         *
         * Returns:
         *  - the merged solution if one exists
         *  - hasSolution = false otherwise
         *
         * Time complexity:
         *  - O(k * log M)
         *    where k is the number of equations and M represents
         *    typical modulus size
         */
        static CRTResult solve(const std::vector<Congruence>& equations) {
            if (equations.empty()) {
                return {false, 0, 0};
            }

            CRTResult current {
            true,
    {
                normalize_mod(equations[0].remainder, equations[0].modulus),
                equations[0].modulus}
            };

            for (std::size_t i = 1; i < equations.size(); i++) {
                Congruence c = current.value;
                current = merge(c, equations[i]);
                if (!current.hasSolution) {
                    return  current;
                }
            }
            return current;
        }
    };
}
