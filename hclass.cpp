#include <cassert>
#include <algorithm>
#include <vector>


template<typename T>
class Matrix {
    public:
        Matrix() {}

        Matrix(const size_t height_, const size_t width_, const T diag_) {
            allocate(height_, width_);
            for (size_t k = 0; k < std::min(height, width); k++) {
                (*this)(k, k) = diag_;
            }
        }
    
        Matrix(const std::vector<std::vector<T>>& vec) {
            allocate(vec.size(), vec[0].size());
            for (size_t i = 0; i < vec.size(); i++) {
                for (size_t j = 0; j < vec[i].size(); j++) {
                    (*this)(i, j) = vec[i][j];
                }
            }
        }

        Matrix(const size_t height_, const size_t width_) {
            allocate(height_, width_);
        }

        Matrix(const Matrix& other) {
            allocate(other.rows(), other.cols());
            for (size_t i = 0; i < other.height; i++) {
                for (size_t j = 0; j < other.width; j++) {
                    (*this)(i, j) = other(i, j);
                }
            }
        }

        Matrix(Matrix&& other) noexcept {
            height = other.rows();
            width = other.cols();
            matrix = other.matrix;
            other.matrix = nullptr;
            other.height = 0;
            other.width = 0;
        }

        Matrix& operator=(const Matrix& other) {
            if (this != &other) {
                assign(other.height, other.width);    
                for (size_t i = 0; i < height; i++) {
                    for (size_t j = 0; j < width; j++) {
                        (*this)(i, j) = other(i, j);
                    }
                }
            }
            return *this;
        }

        Matrix& operator=(Matrix&& other) noexcept {
            if (this != &other) {
                deallocate();
                height = other.height;
                width = other.width;
                matrix = other.matrix;
                other.matrix = nullptr;
                other.height = 0;
                other.width = 0;
            }
            return *this;
        }

        ~Matrix() {
            deallocate();
        }

        T& operator()(const size_t row, const size_t col) noexcept {
            assert(row < rows() && col < cols() && "Index must be correct");
            return matrix[row * cols() + col];
        }

        const T& operator()(const size_t row, const size_t col) const noexcept {
            assert(row < rows() && col < cols() && "Index must be correct");
            return matrix[row * cols() + col];
        }


        bool operator==(const Matrix& otherMatrix) const noexcept {
            if (otherMatrix.rows() != rows() || otherMatrix.cols() != cols())
                return false;
            for (size_t i = 0; i < rows(); i++) {
                for (size_t j = 0; j < cols(); j++) {
                    if (otherMatrix(i, j) != (*this)(i, j))
                        return false;
                }
            }
            return true;
        }

        bool operator!=(const Matrix& otherMatrix) const noexcept {
            return !(*this == otherMatrix);
        }

        const size_t rows() const noexcept {
            return height;
        }

        const size_t cols() const noexcept {
            return width;
        }

        Matrix operator*(const Matrix& other) const {
            assert((rows() == other.cols()) && "Rows of the first matrix must be equal to cols of the second");
            Matrix res(rows(), other.cols(), T());
            for (size_t i = 0; i < rows(); i++) {
                for (size_t j = 0; j < other.cols(); j++) {
                    for (size_t k = 0; k < cols(); k++) {
                        res(i, j) += (*this)(i, k) * other(k, j);
                    }
                }
            }
            return res;
        }

        void toStairseCase() {
            bareissToStairseCase();
        }

        int getDimKernel() const noexcept {
            //!!check stairscase
            for (size_t i = 0; i < cols(); i++) {
                if ((*this)(rows() - 1, i) != 0) {
                    return 0;
                }
            }
            int res = 0;
            for (size_t i = rows(); i > 0; i--) {
                if ((*this)(i - 1, cols() - 1) == 0) {
                    res++;
                } else {
                    return res;
                }
            }
            return res;
        }

        int getDimImage() const noexcept {
            return rows() - getDimKernel();
        }

        void swapCols(const size_t& i, const size_t& j) {
            if (i != j) {
                for (size_t t = 0; t < (*this).rows(); t++) {
                    std::swap((*this)(i, t), (*this)(j, t));
                }
            }
        }

        void swapRows(const size_t& i, const size_t& j) {
            if (i != j) {
                for (size_t t = 0; t < (*this).cols(); t++) {
                    std::swap((*this)(t, i), (*this)(t, j));
                }
            }
        }


    private:
        size_t height = 0, width = 0;
        T * matrix = nullptr;

        void deallocate() {
            height = 0;
            width = 0;
            if (matrix) {
                delete [] matrix;
                matrix = nullptr;
            }
        }

        void allocate(const size_t height_, const size_t width_) {
            height = height_;
            width = width_;
            matrix = new T[height * width]();
        }

        void assign(const size_t height_, const size_t width_) {
            if (height * width != height_ * width_) {
                deallocate();
                allocate(height_, width_);
            } else {
                height = height_;
                width = width_;
            }
            for (size_t i = 0; i < height * width; i++) {
                matrix[i] = T();
            }
        }

        void bareissToStairseCase() {
            T prev = 1;
            size_t row = 0;
            for (size_t k = 0; k < std::min(rows(), cols()); k++) {
                size_t ind = row;
                for (; ind < rows(); ind++) {
                    if ((*this)(ind, k) != 0) {
                        for (size_t j = 0; j < width; j++) {
                            std::swap((*this)(row, j), (*this)(ind, j));
                        }
                        row = ind;
                        break;
                    }
                }

                if (ind == rows())
                    continue;
                
                for (size_t i = row + 1; i < rows(); i++) {
                    for (size_t j = k + 1; j < cols(); j++) {
                        T num = (*this)(row, k) * (*this)(i, j) - (*this)(i, k) * (*this)(row, j);
                        (*this)(i, j) = num / prev;
                    }
                    (*this)(i, k) = 0;
                }
                prev = (*this)(row, k);
                row++;
            }
        }
};


template<typename T>
class Complex {
    public:
        Complex() {}

        Complex(const size_t n) {
            allocate(n);
        }

        Complex(const size_t n, const T& elem) {
            allocate(n);
            for (size_t i = 0; i < n; i++) {
                (*this)[i] = elem;
            }
        }

        Complex(const Complex& other) {
            allocate(other.size());
            for (size_t i = 0; i < other.size(); i++) {
                (*this)[i] = other[i];
            }
        }

        Complex(Complex&& other) noexcept {
            elems = other.size();
            complex = other.complex;
            elems = 0;
            complex = nullptr;
        }

        Complex& operator=(const Complex& other) {
            if (this != &other) {
                assign(other.size());
                for (size_t i = 0; i < other.size(); i++) {
                    (*this)[i] = other[i];
                }
            }
            return *this;
        }

        Complex& operator=(Complex&& other) {
            if (this != &other) {
                deallocate();
                elems = other.size();
                complex = other.complex;
                elems = 0;
                complex = nullptr;
            }
            return *this;
        }

        ~Complex() {
            deallocate();
        }

        T& operator[](const size_t ind) noexcept{
            assert(ind < (*this).size() && "Index must be correct");
            return complex[ind];
        }

        const T& operator[](const size_t ind) const noexcept{
            assert(ind < (*this).size() && "Index must be correct");
            return complex[ind];
        }

        const size_t size() const noexcept{
            return elems;
        }

        bool checkHomology() {
            for (size_t i = 1; i < (*this).size(); i++) {
                if ((*this)[i] * (*this)[i - 1] != T((*this)[i].rows(), (*this)[i - 1].cols(), 0)) {
                    return false;
                }
            }
            return true;
        }

        void toStairseCase() {
            for (size_t i = 0; i < (*this).size(); i++) {
                complex[i].toStairseCase();
            }
        }

        int homology(const size_t ind) {
            assert(ind + 1 < (*this).size() && "Ind must be correct");
            return (*this)[ind].getDimKernel() - (*this)[ind + 1].getDimImage();
        }

    private:
        size_t elems = 0;
        T * complex;
        
        void allocate(size_t n) {
            elems = n;
            complex = new T[elems]();
        }

        void deallocate() {
            elems = 0;
            delete [] complex;
            complex = nullptr;
        }

        void assign(size_t n) {
            if (n != elems) {
                deallocate();
                allocate(n);
            } else {
                for (size_t i = 0; i < elems; i++) {
                    complex[i] = T();
                }
            }
        }
};


