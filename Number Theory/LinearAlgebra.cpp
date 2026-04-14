/*
* Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei
 *
 *
 *  This file provides matrix utilities under modular arithmetic:
 *   - Matrix type alias
 *   - Identity matrix construction
 *   - Matrix multiplication modulo m
 *   - Fast matrix exponentiation modulo m
 *
 * Complexity summary:
 *   - identity(n): O(n^2)
 *   - multiply(A, B): O(n^3) for square n x n matrices
 *   - power(A, e): O(n^3 * log e) for square n x n matrices
 */

#include "NumberTheory.cpp"

namespace LinearAlgebra {
    template <typename T>
    using Matrix = std::vector<std::vector<T>>;

    /*
     * MatrixArithmetic
     *
     * Provides modular matrix operations under a fixed modulus.
     */
    template <typename T>
    class MatrixArithmetic {
        T mod;

    public:
        explicit MatrixArithmetic(T modulus) : mod(modulus) {}

        /*
         * Normalizes x into [0, mod).
         * Time complexity: O(1)
         */
        T normalize(T x) const {
            return NumberTheory::normalize_mod(x, mod);
        }

        /*
         * Constructs the n x n identity matrix.
         *
         * Time complexity:
         *  - O(n^2)
         */
        Matrix<T> identity(std::size_t n) const {
            Matrix<T> result(n, std::vector<T>(n));
            for (std::size_t i = 0; i < n; i++) {
                result[i][i] = 1;
            }
            return result;
        }

        /*
         * Multiplies matrices A and B modulo mod.
         *
         * If A is (rows x inner) and B is (inner x cols),
         * the result is (rows x cols).
         *
         * Time complexity:
         *  - O(rows * inner * cols)
         *  - O(n^3) for square n x n matrices
         */
        Matrix<T> multiply(Matrix<T>& A, Matrix<T>& B) const {
            const std::size_t rows = A.size();
            const std::size_t inner = A[0].size();
            const std::size_t cols = B[0].size();

            Matrix<T> result(rows, std::vector<T>(cols, 0));

            for (std::size_t i = 0; i < rows; i++) {
                for (std::size_t k = 0; k < inner; k++) {
                    T a = normalize(A[i][k]);
                    if (a == 0) continue;

                    for (std::size_t j = 0; j < cols; j++) {
                        T b = normalize(B[k][j]);
                        if (b == 0) continue;

                        result[i][j] = normalize(result[i][j] + a * b);
                    }
                }
            }

            return result;
        }

        /*
         * Computes base^exponent modulo mod using binary exponentiation.
         *
         * Preconditions:
         *  - base must be a square matrix
         *
         * Time complexity:
         *  - O(n^3 * log exponent) for an n x n matrix
         */
        Matrix<T> power(Matrix<T> base, long long int exponent) const {
            const std::size_t n = base.size();
            Matrix<T> result = identity(n);

            while (exponent > 0) {
                if (exponent & 1LL) {
                    result = multiply(result, base);
                }
                base = multiply(base, base);
                exponent >>= 1LL;
            }

            return result;
        }
    };
}