#include "biginteger.h"

template <size_t N>
class Residue {
  private:
    size_t val;

  public:
    Residue();
    size_t get_val() const;
    explicit Residue(int x);
    explicit operator int() const;
    Residue<N>& operator+=(const Residue<N>& other);
    Residue<N>& operator-=(const Residue<N>& other);
    Residue<N>& operator*=(const Residue<N>& other);
    Residue<N>& operator/=(const Residue<N>& other);
};

template <size_t N, size_t M, typename Field = Rational>
class Matrix {
  private:
    std::array<std::array<Field, M>, N> v;

  public:
    Matrix() = default;
    Matrix(std::initializer_list<std::vector<int>> a);
    Field gauss_and_det();
    void inverse_gauss();
    Field det() const;
    size_t rank() const;
    Matrix<N, M, Field> inverted() const;
    void invert();
    std::array<Field, M>& operator[](size_t i);
    const std::array<Field, M>& operator[](size_t i) const;
    Field trace() const;
    std::array<Field, M> getRow(size_t i) const;
    std::array<Field, N> getColumn(size_t i) const;
    Matrix<M, N, Field> transposed() const;
    Matrix<N, M, Field>& operator+=(const Matrix<N, M, Field>& other);
    Matrix<N, M, Field>& operator-=(const Matrix<N, M, Field>& other);
    Matrix<N, M, Field>& operator*=(const Field& a);
    Matrix<N, M, Field>& operator*=(const Matrix<N, M, Field>& other);
};

size_t mod(long long x, size_t N) {
    long long tmp = static_cast<long long>(N);
    return (x % tmp + tmp) % tmp;
}

template <size_t N, size_t M, typename Field>
bool operator==(const Matrix<N, M, Field>& a, const Matrix<N, M, Field>& b);

template <size_t N, typename Field>
void minus(std::array<Field, N>& a, const std::array<Field, N>& b,
           const Field& k);

template <size_t N, typename Field>
void multiply(std::array<Field, N>& a, const Field& k);

template <size_t N>
bool operator==(const Residue<N>& a, const Residue<N>& b);

template <size_t N>
bool operator!=(const Residue<N>& a, const Residue<N>& b);

template <size_t N>
std::ostream& operator>>(std::ostream& out, Residue<N>& a);

template <size_t N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& a);

template <size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

template <size_t N>
bool operator==(const Residue<N>& a, const Residue<N>& b) {
    return a.get_val() == b.get_val();
}

template <size_t N>
bool operator!=(const Residue<N>& a, const Residue<N>& b) {
    return !(a == b);
}

template <size_t N, size_t M, size_t K, typename Field>
Matrix<N, K, Field> operator*(const Matrix<N, M, Field>& a,
                              const Matrix<M, K, Field>& b) {
    Matrix<N, K, Field> ans;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < K; ++j) {
            for (size_t k = 0; k < M; ++k) {
                ans[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return ans;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=(const Field& a) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            v[i][j] *= a;
        }
    }
    return *this;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(const Matrix<N, M, Field>& m, const Field& a) {
    Matrix<N, M, Field> ans = m;
    ans *= a;
    return ans;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(const Field& a, const Matrix<N, M, Field>& m) {
    Matrix<N, M, Field> ans = m;
    ans *= a;
    return ans;
}

template <size_t N, size_t M, typename Field>
Matrix<M, N, Field> Matrix<N, M, Field>::transposed() const {
    Matrix<M, N, Field> ans;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            ans[j][i] = v[i][j];
        }
    }
    return ans;
}

template <size_t N, size_t M, typename Field>
Field Matrix<N, M, Field>::trace() const {
    static_assert(N == M);
    Field ans = Field(0);
    for (size_t i = 0; i < N; ++i) {
        ans += v[i][i];
    }
    return ans;
}

template <size_t N, size_t M, typename Field>
std::array<Field, M>& Matrix<N, M, Field>::operator[](size_t i) {
    return v[i];
}

template <size_t N, size_t M, typename Field>
const std::array<Field, M>& Matrix<N, M, Field>::operator[](size_t i) const {
    return v[i];
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>& plus_or_minus(Matrix<N, M, Field>& a,
                                   const Matrix<N, M, Field>& b, bool plus) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            if (plus) {
                a[i][j] += b[i][j];
            } else {
                a[i][j] -= b[i][j];
            }
        }
    }
    return a;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator+=(
    const Matrix<N, M, Field>& other) {
    return plus_or_minus(*this, other, true);
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator-=(
    const Matrix<N, M, Field>& other) {
    return plus_or_minus(*this, other, false);
}

template <size_t N, size_t M, typename Field>
std::array<Field, M> Matrix<N, M, Field>::getRow(size_t i) const {
    return v[i];
}

template <size_t N, size_t M, typename Field>
std::array<Field, N> Matrix<N, M, Field>::getColumn(size_t i) const {
    std::array<Field, N> ans;
    for (size_t j = 0; j < N; ++j) {
        ans[j] = v[j][i];
    }
    return ans;
}

template <size_t N, size_t M, typename Field>
Field Matrix<N, M, Field>::gauss_and_det() {
    int inversions_count = 0;
    size_t j = 0;
    for (size_t i = 0; i < M; ++i) {
        if (j >= N) {
            break;
        }
        if (v[j][i] == Field(0)) {
            for (size_t k = j + 1; k < N; ++k) {
                if (v[k][i] != Field(0)) {
                    std::swap(v[k], v[j]);
                    ++inversions_count;
                    break;
                }
            }
        }
        if (v[j][i] != Field(0)) {
            for (size_t k = j + 1; k < N; ++k) {
                minus(v[k], v[j], v[k][i] / v[j][i]);
            }
            ++j;
        }
    }
    if (N != M) {
        return Field(0);
    }
    Field ans = Field(1);
    for (size_t i = 0; i < N; ++i) {
        ans *= v[i][i];
    }
    if (inversions_count % 2 != 0) {
        ans *= Field(-1);
    }
    return ans;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>::Matrix(std::initializer_list<std::vector<int>> a) {
    size_t i = 0;
    for (const auto* it = a.begin(); it != a.end(); ++it) {
        for (size_t j = 0; j < M; ++j) {
            v[i][j] = Field((*it)[j]);
        }
        ++i;
    }
}

template <size_t N, size_t M, typename Field>
Field Matrix<N, M, Field>::det() const {
    Matrix<N, M, Field> tmp = *this;
    return tmp.gauss_and_det();
}

template <size_t N, size_t M, typename Field>
size_t Matrix<N, M, Field>::rank() const {
    Matrix<N, M, Field> tmp = *this;
    tmp.gauss_and_det();
    size_t ans = 0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            if (tmp[i][j] != Field(0)) {
                ++ans;
                break;
            }
        }
    }
    return ans;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> Matrix<N, M, Field>::inverted() const {
    static_assert(N == M);
    Matrix<N, 2 * N, Field> tmp;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            tmp[i][j] = v[i][j];
        }
    }
    for (size_t i = 0; i < N; ++i) {
        tmp[i][N + i] = Field(1);
    }
    tmp.gauss_and_det();
    tmp.inverse_gauss();
    Matrix<N, N, Field> ans;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            ans[i][j] = tmp[i][N + j];
        }
    }
    return ans;
}

template <size_t N, size_t M, typename Field>
void Matrix<N, M, Field>::inverse_gauss() {
    for (size_t i = 0; i < N; ++i) {
        Field x = Field(1) / v[i][i];
        multiply(v[i], x);
    }
    for (size_t i = 1; i < N; ++i) {
        for (size_t j = 0; j < i; ++j) {
            Field k = v[j][i];
            minus(v[j], v[i], k);
        }
    }
}

template <size_t N, size_t M, typename Field>
void Matrix<N, M, Field>::invert() {
    *this = inverted();
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=(
    const Matrix<N, M, Field>& other) {
    *this = *this * other;
    return *this;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator+(const Matrix<N, M, Field>& a,
                              const Matrix<N, M, Field>& b) {
    Matrix<N, M, Field> ans = a;
    ans += b;
    return ans;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator-(const Matrix<N, M, Field>& a,
                              const Matrix<N, M, Field>& b) {
    Matrix<N, M, Field> ans = a;
    ans -= b;
    return ans;
}

template <size_t N>
Residue<N>::operator int() const {
    return val;
}

template <size_t N, size_t D>
struct IsPrimeHelper;

template <bool B, size_t N, size_t D>
struct Conditional {
    static const size_t value = (N % D != 0) && IsPrimeHelper<N, D + 1>::value;
};

template <size_t N, size_t D>
struct Conditional<true, N, D> {
    static const bool value = true;
};

template <size_t N, size_t D>
struct IsPrimeHelper {
    static const bool value = Conditional<(D * D > N), N, D>::value;
};

template <size_t N>
struct IsPrime {
    static const bool value = IsPrimeHelper<N, 2>::value;
};

template <size_t N>
Residue<N>::Residue(int x) : val(mod(x, N)) {}

template <size_t N>
Residue<N>& Residue<N>::operator+=(const Residue<N>& other) {
    val += other.val;
    val %= N;
    return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator-=(const Residue<N>& other) {
    val =
        mod(static_cast<long long>(val) - static_cast<long long>(other.val), N);
    return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator*=(const Residue<N>& other) {
    val *= other.val;
    val %= N;
    return *this;
}

size_t PowMod(size_t a, size_t k, size_t N) {
    if (k == 0) {
        return 1;
    }
    if (k % 2 == 0) {
        size_t ans = PowMod(a, k / 2, N);
        return ans * ans % N;
    }
    return PowMod(a, k - 1, N) * a % N;
}

template <size_t N>
Residue<N>& Residue<N>::operator/=(const Residue<N>& other) {
    static_assert(IsPrime<N>::value);
    val *= PowMod(other.val, N - 2, N);
    val %= N;
    return *this;
}

template <size_t N>
Residue<N>::Residue() : Residue(0) {}

template <size_t N>
size_t Residue<N>::get_val() const {
    return val;
}

template <size_t N>
Residue<N> operator+(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> tmp = a;
    tmp += b;
    return tmp;
}

template <size_t N>
Residue<N> operator-(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> tmp = a;
    tmp -= b;
    return tmp;
}

template <size_t N>
Residue<N> operator*(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> tmp = a;
    tmp *= b;
    return tmp;
}

template <size_t N>
Residue<N> operator/(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> tmp = a;
    tmp /= b;
    return tmp;
}

template <size_t N>
std::istream& operator>>(std::istream& in, Residue<N>& a) {
    int x;
    in >> x;
    a = Residue<N>(x);
    return in;
}

template <size_t N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& a) {
    out << a.get_val();
    return out;
}

template <size_t N, typename Field>
void minus(std::array<Field, N>& a, const std::array<Field, N>& b,
           const Field& k) {
    for (size_t i = 0; i < N; ++i) {
        a[i] -= b[i] * k;
    }
}

template <size_t N, typename Field>
void multiply(std::array<Field, N>& a, const Field& k) {
    for (size_t i = 0; i < N; ++i) {
        a[i] *= k;
    }
}

template <size_t N, size_t M, typename Field>
bool operator==(const Matrix<N, M, Field>& a, const Matrix<N, M, Field>& b) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            if (a[i][j] != b[i][j]) {
                return false;
            }
        }
    }
    return true;
}
