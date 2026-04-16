/*
* Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 *
 *  This file provides matrix utilities under modular arithmetic:
 *   - Matrix type alias
 *   - Identity matrix construction
 *   - Matrix multiplication modulo m
 *   - Fast matrix exponentiation modulo m
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
        NumberTheory::ModularArithmetic<T> modular_arithmetic;

    public:
        explicit MatrixArithmetic(T modulus) : mod(modulus), modular_arithmetic(modulus) {}

        /*
         * Constructs the n x n identity matrix.
         *
         * Time complexity: O(n^2)
         * Space complexity: O(n^2)
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
         * Preconditions:
         *  - A and B must be non-empty
         *  - A[0].size() == B.size()
         *
         * Time complexity: O(rows * inner * cols)
         * Space complexity: O(rows * cols)
         */
        Matrix<T> multiply(Matrix<T>& A, Matrix<T>& B) const {
            const std::size_t rows = A.size();
            const std::size_t inner = A[0].size();
            const std::size_t cols = B[0].size();

            Matrix<T> result(rows, std::vector<T>(cols, 0));

            for (std::size_t i = 0; i < rows; i++) {
                for (std::size_t k = 0; k < inner; k++) {
                    T a = modular_arithmetic.normalize(A[i][k]);
                    if (a == 0) continue;

                    for (std::size_t j = 0; j < cols; j++) {
                        T b = modular_arithmetic.normalize(B[k][j]);
                        if (b == 0) continue;

                        result[i][j] = modular_arithmetic.normalize(result[i][j] + a * b);
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
         *  - exponent >= 0
         *
         * Time complexity: O(n^3 * log exponent) for an n x n matrix
         * Space complextiy: O(n^2)
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