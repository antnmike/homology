#include "../chaincomplexz.cpp"

int main() {
    size_t t;
    std::cin >> t;
    ChainComplex a;
    a.resize(t);
    for (size_t k = 0; k < t; k++) {
        slong n, m;
        std::cin >> n >> m;
        fmpz_mat_init(a[k].mat, n, m);
        for (slong i = 0; i < n; i++) {
            for (slong j = 0; j < m; j++) {
                slong x;
                std::cin >> x;
                fmpz_set_si(fmpz_mat_entry(a[k].mat, i, j), x);
            }
        }
    }
    std::cout << a.getHomology(1) << std::endl;
    return 0;
}