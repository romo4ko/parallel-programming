#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <iomanip>

template <typename T = double>
class MatrixDiagonal : public Matrix<T> {
private:
    unsigned _dimension; // Matrix size (number of rows/columns)
    T* _elements;      // Array to store diagonal elements

public:
    MatrixDiagonal(unsigned dimension) : _dimension(dimension) {
        if (dimension == 0) {
            throw std::invalid_argument("Matrix size must be positive.");
        }
        _elements = new T[_dimension];
        for(size_t i = 0; i < _dimension; ++i)
            _elements[i] = 0;
    }

    MatrixDiagonal(const MatrixDiagonal& other) : _dimension(other._dimension) {
        _elements = new T[_dimension];
        for(size_t i = 0; i < _dimension; ++i){
            _elements[i] = other._elements[i];
        }
    }

    MatrixDiagonal& operator=(const MatrixDiagonal& other) {
        if (this == &other) return *this;
        delete[] _elements;
        _dimension = other._dimension;
        _elements = new T[_dimension];
        for(size_t i = 0; i < _dimension; ++i){
            _elements[i] = other._elements[i];
        }
        return *this;
    }

    ~MatrixDiagonal() override {
        delete[] _elements;
    }

    T& operator()(unsigned row, unsigned col) override {
        if (row >= _dimension || col >= _dimension) {
            throw std::out_of_range("Matrix index out of bounds.");
        }
        if (row == col) {
            return _elements[row];
        } else {
            static T zero = 0;
            return zero; // Return reference to static variable 0
        }
    }

    const T& operator()(unsigned row, unsigned col) const override {
        if (row >= _dimension || col >= _dimension) {
            throw std::out_of_range("Matrix index out of bounds.");
        }
        if (row == col) {
            return _elements[row];
        } else {
            static T zero = 0;
            return zero; // Return reference to static variable 0
        }
    }

    Matrix<T>* operator+(const Matrix<T>& other) const override {
        const MatrixDiagonal<T>* otherDiag = dynamic_cast<const MatrixDiagonal<T>*>(&other);
        if (!otherDiag) {
            throw std::invalid_argument("Incompatible matrix types for addition.");
        }

        if(_dimension != otherDiag->_dimension) {
            throw std::invalid_argument("Matrix dimensions must be equal for addition.");
        }

        MatrixDiagonal<T>* result = new MatrixDiagonal<T>(_dimension);
        for(size_t i = 0; i < _dimension; ++i) {
            (*result)(i, i) = (*this)(i, i) + (*otherDiag)(i, i);
        }

        return result;
    }

    Matrix<T>* operator-(const Matrix<T>& other) const override {
        const MatrixDiagonal<T>* otherDiag = dynamic_cast<const MatrixDiagonal<T>*>(&other);
        if (!otherDiag) {
            throw std::invalid_argument("Incompatible matrix types for subtraction.");
        }

        if(_dimension != otherDiag->_dimension) {
            throw std::invalid_argument("Matrix dimensions must be equal for subtraction.");
        }

        MatrixDiagonal<T>* result = new MatrixDiagonal<T>(_dimension);
        for(size_t i = 0; i < _dimension; ++i) {
            (*result)(i, i) = (*this)(i, i) - (*otherDiag)(i, i);
        }

        return result;
    }

    Matrix<T>* operator*(const Matrix<T>& other) const override {
        const MatrixDiagonal<T>* otherDiag = dynamic_cast<const MatrixDiagonal<T>*>(&other);
        if (!otherDiag) {
            throw std::invalid_argument("Incompatible matrix types for multiplication.");
        }

        if(_dimension != otherDiag->_dimension) {
            throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
        }

        MatrixDiagonal<T>* result = new MatrixDiagonal<T>(_dimension);
        for(size_t i = 0; i < _dimension; ++i) {
            (*result)(i, i) = (*this)(i, i) * (*otherDiag)(i, i);
        }

        return result;
    }

    Matrix<T>* elementWiseMultiplication(const Matrix<T>& other) const override {
        return this->operator*(other);
    }

    Matrix<T>* transpose() const override {
        return new MatrixDiagonal<T>(*this);
    }

    void importFromFile(const std::string& filename) override {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        std::string matrixType;
        file >> matrixType;
        if(matrixType != "MatrixDiagonal"){
            throw std::runtime_error("Invalid matrix type in file: " + filename);
        }

        unsigned dimension;
        file >> dimension;
        if (file.fail() || dimension == 0) {
            throw std::runtime_error("Invalid matrix dimensions in file: " + filename);
        }
        delete[] _elements;
        _dimension = dimension;
        _elements = new T[_dimension];

        for (unsigned i = 0; i < _dimension; ++i) {
            if (!(file >> _elements[i])) {
                throw std::runtime_error("Error reading matrix data from file: " + filename);
            }
        }

        file.close();
    }

    void exportToFile(const std::string& filename) const override {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        file << "MatrixDiagonal\n";
        file << _dimension << "\n";

        for (unsigned i = 0; i < _dimension; ++i) {
            file << _elements[i] << " ";
        }
        file << "\n";

        file.close();
    }

    void print() const override {
        for (unsigned i = 0; i < _dimension; ++i) {
            for (unsigned j = 0; j < _dimension; ++j) {
                if (i == j) {
                    std::cout << std::setw(10) << _elements[i] << " ";
                } else {
                    std::cout << std::setw(10) << 0 << " ";
                }
            }
            std::cout << std::endl;
        }
    }
};