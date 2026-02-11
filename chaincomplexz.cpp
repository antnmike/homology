#include <flint/fmpz_mat.h>
#include <flint/fmpz.h>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>

constexpr std::string_view PRSYM = "⊕";
constexpr std::string_view ZSYM = "Z";
constexpr std::string_view LOWERSYM[] {"₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"};
constexpr std::string_view UPPERSYM[] {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};

struct HZi {
    slong k;
    slong cnt;
};

struct Hres {
    slong torsion_value;
    std::vector<HZi> count;
};

std::string getZ(slong up, slong low) {
    std::string res;
    res += ZSYM;
    std::string tmp;
    if (low != 1) {
        tmp = std::to_string(low);
        for (char i: tmp) {
            res += LOWERSYM[i - '0'];
        }
    }
    if (up != 1) {
        tmp = std::to_string(up);
        for (char i: tmp) {
            res += UPPERSYM[i - '0'];
        }
    }
    return res;
}


std::ostream& operator<<(std::ostream &os, const Hres& a) {
    bool fl = 0;
    if (a.torsion_value == 0 && a.count.empty()) {
        return os << "0";
    } else {
        if (a.torsion_value != 0) {
            fl = 1;
            os << getZ(a.torsion_value, 1);
        }
    }
    for (HZi i: a.count) {
        if (fl) {
            os << " " << PRSYM << " ";
            fl = 1;
        }
        os << getZ(i.cnt, i.k);
    }
    return os;
}


struct FlintMatrix {
    fmpz_mat_t mat;

    FlintMatrix() {
        fmpz_mat_init(mat, 0, 0);
    }

    ~FlintMatrix() {
        fmpz_mat_clear(mat);
    }

    FlintMatrix(const FlintMatrix&) = delete;
    FlintMatrix& operator=(const FlintMatrix&) = delete;

    FlintMatrix(FlintMatrix&& other) noexcept {
        fmpz_mat_init(mat, 0, 0);
        std::swap(mat[0], other.mat[0]);
    }
};


struct ChainComplex {
    public:
        ChainComplex() {}

        ChainComplex(ChainComplex&& other) = default;

        ChainComplex(const ChainComplex& other) {
            comp.resize(other.size());
            smith.resize(smith.size());
            issmith.assign(issmith.size(), 0);
            for (size_t i = 0; i < other.size(); i++) {
                fmpz_mat_init(comp[i].mat, fmpz_mat_nrows(other[i].mat), fmpz_mat_ncols(other[i].mat));
                fmpz_mat_set(comp[i].mat, other[i].mat);
            }
        }

        ChainComplex& operator=(ChainComplex other) {
            std::swap(comp, other.comp);
            std::swap(smith, other.smith);
            std::swap(issmith, other.issmith);
            return *this;
        }

        ~ChainComplex() {
            clear();
        }

        size_t size() const{
            return comp.size();
        }

        void resize(size_t n) {
            comp.resize(n);
            smith.resize(n);
            issmith.resize(n);
        }

        FlintMatrix& operator[](const size_t& ind) {
            issmith[ind] = 0;
            return comp[ind];
        }

        const FlintMatrix& operator[](size_t ind) const {
            return comp[ind];
        }

        Hres getHomology(size_t ind) {
            assert(ind >= 1 && ind < comp.size() && "ind must be >= 1 and < n");
            assert(checkHomology(comp[ind - 1], comp[ind]) && "Homology isnt correct");
            toSmithForm(ind - 1);
            toSmithForm(ind);
            return solveHomology(smith[ind - 1], smith[ind]);
        }

        void clear() {
            comp.clear();
            issmith.clear();
            smith.clear();
        }

        bool checkComplex() {
            for (size_t i = 1; i < comp.size(); i++) {
                if (checkHomology(comp[i - 1], comp[i])) {
                    return 0;
                }
            }
            return 1;
        }

    private:
        std::vector<FlintMatrix> comp;
        std::vector<FlintMatrix> smith;
        std::vector<bool> issmith;

        void toSmithForm(size_t i) {
            if (!issmith[i]) {
                fmpz_mat_init(smith[i].mat, comp[i].mat->r, comp[i].mat->c);
                fmpz_mat_snf(smith[i].mat, comp[i].mat);
                issmith[i] = 1;
            }
        }

        bool checkAllZero(FlintMatrix& A) {
            for (slong i = 0; i < A.mat->r; i++) {
                for (slong j = 0; j < A.mat->c; j++) {
                    if (!fmpz_is_zero(fmpz_mat_entry(A.mat, i, j)))
                        return 0;
                }
            }
            return 1;
        }

        bool checkHomology(FlintMatrix& A, FlintMatrix& B) {
            if (B.mat->r == A.mat->c) {
                FlintMatrix c;
                fmpz_mat_init(c.mat, B.mat->r, A.mat->c);
                fmpz_mat_mul(c.mat, B.mat, A.mat);
                return checkAllZero(c);
            }
            return 0;
        }

        Hres solveHomology(FlintMatrix& A, FlintMatrix& B) {
            slong cnta = 0, cntb = 0;
            slong minDimA = std::min(A.mat->r, A.mat->c), minDimB = std::min(B.mat->r, B.mat->c);
            for (slong i = 0; i < minDimA; i++) {
                if (!fmpz_is_zero(fmpz_mat_entry(A.mat, i, i))) {
                    cnta++;
                }
            }

            std::vector<HZi> a;
            for (int i = 0; i < minDimB; i++) {
                if (!fmpz_is_zero(fmpz_mat_entry(B.mat, i, i))) {
                    cntb++;
                    if (!fmpz_is_one(fmpz_mat_entry(B.mat, i, i))) {
                        if (!a.empty() && a.back().k == fmpz_get_si(fmpz_mat_entry(B.mat, i, i))) {
                            a.back().cnt++;
                        } else {
                            a.push_back({fmpz_get_si(fmpz_mat_entry(B.mat, i, i)), 1});
                        }
                    }
                }
            }
            return {A.mat->c - cnta - cntb, a};
        }
};