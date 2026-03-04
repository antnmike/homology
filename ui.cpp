#include "chaincomplexz.cpp"
#include <cstddef>

int main() {
    std::cout << "Enter the size of the complex: ";
    std::cout.flush();
    size_t t;
    std::cin >> t;
    ChainComplex cmp;
    cmp.resize(t);
    for (size_t k = 0; k < t; k++) {
        std::cout << "Enter height and width: ";
        std::cout.flush();
        slong n, m;
        std::cin >> n >> m;
        fmpz_mat_init(cmp[k].mat, n, m);
        std::cout << "Enter the complex: " << std::endl;
        for (slong i = 0; i < n; i++) {
            for (slong j = 0; j < m; j++) {
                slong x;
                std::cin >> x;
                fmpz_set_si(fmpz_mat_entry(cmp[k].mat, i, j), x);
            }
        }
    }
    if (!cmp.checkComplex()) {
        std::cout << "Invalid complex" << std::endl;
    } else {
        for (size_t i = 1; i < t; i++) {
            std::cout << cmp.getHomology(i) << std::endl;
        }
    }
    return 0;
}