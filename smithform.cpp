#include <flint/fmpz_mat.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <iostream>
#include <string>
#include <vector>

const std::string PRSYM = "⊕";
const std::string ZSYM = "Z";
const std::string LOWERSYM[] {"₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"};
const std::string UPPERSYM[] {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};

struct HZi {
    slong k;
    int cnt;
};

struct Hres {
    int t;
    std::vector<HZi> az;
};

std::string getZ(int up, int low) {
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
    std::cout << a.t << std::endl;
    if (a.t == 0) {
        return os << "0";
    }
    os << getZ(a.t, 1);
    for (HZi i: a.az) {
        os << " " << PRSYM << " ";
        os << getZ(i.cnt, i.k);
    }
    return os;
}

 
Hres solveHomology(fmpz_mat_t& A, fmpz_mat_t& B) {
    int cnta = 0, cntb = 0;
    for (int i = 0; i < std::min((int)A->r, (int)A->c); i++) {
        if (fmpz_get_si(fmpz_mat_entry(A, i, i)) != 0) {
            cnta++;
        }
    }
    std::vector<HZi> a;
    for (int i = 0; i < std::min((int)B->r, (int)B->c); i++) {
        if (fmpz_get_si(fmpz_mat_entry(B, i, i)) != 0) {
            cntb++;
        }
        if (fmpz_get_si(fmpz_mat_entry(B, i, i)) >= 2) {
            if (!a.empty() && a.back().k == fmpz_get_si(fmpz_mat_entry(B, i, i))) {
                a.back().cnt++;
            } else {
                a.push_back({fmpz_get_si(fmpz_mat_entry(B, i, i)), 1});
            }
        }
    }
    return {(int)A->c - cnta - cntb, a};
}

void inMatrix(fmpz_mat_t& res) {
    slong rows, cols;
    std::cin >> rows >> cols;
    fmpz_mat_init(res, rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int t;
            std::cin >> t;
            fmpz_set_si(fmpz_mat_entry(res, i, j), (slong)t);
        }
    }
 
}


int main() {
    fmpz_mat_t A, B, Sa, Sb;
    inMatrix(A);
    inMatrix(B);
    fmpz_mat_init(Sa, A->r, A->c);
    fmpz_mat_init(Sb, B->r, B->c);

    fmpz_mat_snf(Sa, A);
    fmpz_mat_snf(Sb, B);

    std::cout << solveHomology(Sa, Sb) << std::endl;

    fmpz_mat_clear(A);
    fmpz_mat_clear(Sa);
    fmpz_mat_clear(B);
    fmpz_mat_clear(Sb);
    return 0;
}
