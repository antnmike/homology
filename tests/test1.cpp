#include <iostream>
#include "../hclass.cpp"


int main() {
    Matrix<int> c1 ({{5, 6}, {7, -8}});
    Matrix<int> c2 ({{4, 3}, {-2, 9}});
    Matrix<int> c3 = c1 * c2;

    std::cout << "Matrixes: " << std::endl;

    for (size_t i = 0; i < c1.rows(); i++) {
        for (size_t j = 0; j < c1.cols(); j++) {
            std::cout << c1(i, j) << " ";
        }
        std::cout << std::endl;
    }

    for (size_t i = 0; i < c2.rows(); i++) {
        for (size_t j = 0; j < c2.cols(); j++) {
            std::cout << c2(i, j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Mul matrix:" << std::endl;
    for (size_t i = 0; i < c3.rows(); i++) {
        for (size_t j = 0; j < c3.cols(); j++) {
            std::cout << c3(i, j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Stairse case:" << std::endl;
    c1.toStairseCase();
    c2.toStairseCase();
    for (size_t i = 0; i < c1.rows(); i++) {
        for (size_t j = 0; j < c1.cols(); j++) {
            std::cout << c1(i, j) << " ";
        }
        std::cout << std::endl;
    }

    for (size_t i = 0; i < c2.rows(); i++) {
        for (size_t j = 0; j < c2.cols(); j++) {
            std::cout << c2(i, j) << " ";
        }
        std::cout << std::endl;
    }
}