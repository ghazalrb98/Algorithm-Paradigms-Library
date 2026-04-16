/*
* Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 *  This file provides reusable number theory utilities:
 *   - Extended Euclidean algorithm
 *   - Greatest common divisor (gcd)
 *   - Least common multiple (lcm)
 *   - Modular arithmetic (add/sub/mul/div/normalize/inverse)
 *   - Chinese Remainder Theorem solver for arbitrary moduli
 *
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
     * Time complexity: O(log(min(a, b)))
     * Space complexity: O(log(min(a, b)))
     *
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

    /*
     * Greatest Common Divisor
     *
     * Time complexity: O(log(min(a, b)))
     * Space complexity: O(log(min(a, b)))
     */
    template <typename T>
    T gcd(T a, T b) {
        T x, y;
        return egcd(a, b, x, y);
    }

    /*
     * Least Common Multiple
     *
     * Time complexity: O(log(min(a, b)))
     * Space complexity: O(log(min(a, b)))
     */
    template <typename T>
    T lcm(T a, T b) {
        return (a / gcd(a, b)) * b;
    }

    /*
     * ModularArithmetic
     *
     * A helper class for modular arithmetic under a fixed modulus.
     */
    template <typename T>
    class ModularArithmetic {
        T mod;

    public:
        explicit ModularArithmetic(T n) : mod(n) {}

        /*
         * Modular normalization
         * Returns x reduced into the interval [0, mod).
         *
         * Time complexity: O(1)
         * Space complexity: O(1)
         */
        T normalize(T x) const {
            x %= mod;
            if (x < 0) x += mod;
            return x;
        }

        /*
         * Modular addition
         *
         * Time complexity: O(1)
         * Space complexity: O(1)
         */
        T add(T a, T b) const {
            return normalize(normalize(a) + normalize(b));
        }

        /*
         * Modular subtraction
         *
         * Time complexity: O(1)
         * Space complexity: O(1)
         */
        T sub(T a, T b) const {
            return normalize(normalize(a) - normalize(b));
        }

        /*
         * Modular multiplication
         *
         * Time complexity: O(1)
         * Space complexity: O(1)
         */
        T mul(T a, T b) const {
            return normalize(normalize(a) * normalize(b));
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
         * Time complexity: O(log(min(a, mod)))
         * Space complexity: O(log(min(a, mod)))
         */
        T inverse(T a) const {
            a = normalize(a);

            T x, y;
            T g = egcd(a, mod, x, y);

            if (g != 1) return -1;
            return normalize(x);        }

        /*
         * Modular division
         *
         * This is implemented as:
         *      (a * b^{-1}) % mod
         *
         * Returns:
         *  - the result if inverse(b) exists
         *  - -1 otherwise
         *
         * Time complexity: O(log(min(b, mod))
         * Space complexity: O(log(min(b, mod))
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
         * Time complexity: O(log(min(n, m)))
         * Space complexity: O(log(min(n, m)))
         */
        static CRTResult merge(const Congruence& c1, const Congruence& c2) {
            ModularArithmetic<T> modular_arithmetic_c1(c1.modulus);
            T a = modular_arithmetic_c1.normalize(c1.remainder);
            T n = c1.modulus;

            ModularArithmetic<T> modular_arithmetic_c2(c2.modulus);
            T b = modular_arithmetic_c2.normalize(c2.remainder);
            T m = c2.modulus;

            T g = gcd(n, m);

            if ((b - a) % g != 0) {
                return {false, {0, 0}};
            }

            T n1 = n / g;
            T m1 = m / g;
            T diff = (b - a) / g;

            ModularArithmetic<T> modular_arithmetic_m1(m1);
            T inv = modular_arithmetic_m1.inverse(n1);
            if (inv == -1) {
                return {false, {0, 0}};
            }

            T t = modular_arithmetic_m1.normalize(diff * inv);

            T new_mod = n * m1;
            ModularArithmetic<T> modular_arithmetic_res(new_mod);
            T remainder = modular_arithmetic_res.normalize(a + n * t);

            return {true, {remainder, new_mod}};
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
         *    where k is the number of equations and
         *    M represents typical modulus size
         *
         *  Space complexity: O(log M)
         */
        static CRTResult solve(const std::vector<Congruence>& equations) {
            if (equations.empty()) {
                return {false,
                    {0, 0}
                };
            }

            ModularArithmetic<T> modular_arithmetic(equations[0].modulus);
            CRTResult current {
            true,
    {
                modular_arithmetic.normalize(equations[0].remainder),
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
