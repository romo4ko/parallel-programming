#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <iomanip>

template <typename T = double>
class MatrixDense : public Matrix<T> {
private:
    unsigned _rows, _cols;
    T* _data;

public:
    MatrixDense(unsigned rows, unsigned cols) : _rows(rows), _cols(cols) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("Matrix dimensions must be positive.");
        }
        _data = new T[_rows * _cols];
        for (size_t i = 0; i < _rows * _cols; ++i)
            _data[i] = 0; // Initialize with zeros
    }

    MatrixDense(const MatrixDense& other) : _rows(other._rows), _cols(other._cols) {
        _data = new T[_rows * _cols];
        for(size_t i = 0; i < _rows * _cols; ++i){
            _data[i] = other._data[i];
        }
    }

    MatrixDense& operator=(const MatrixDense& other) {
        if (this == &other) return *this;
        delete[] _data;
        _rows = other._rows;
        _cols = other._cols;
        _data = new T[_rows * _cols];
        for(size_t i = 0; i < _rows * _cols; ++i){
            _data[i] = other._data[i];
        }
        return *this;
    }

    ~MatrixDense() override {
        delete[] _data;
    }

    T& operator()(unsigned row, unsigned col) override {
        if (row >= _rows || col >= _cols) {
            throw std::out_of_range("Matrix index out of bounds.");
        }
        return _data[row * _cols + col];
    }

    const T& operator()(unsigned row, unsigned col) const override {
        if (row >= _rows || col >= _cols) {
            throw std::out_of_range("Matrix index out of bounds.");
        }
        return _data[row * _cols + col];
    }

    Matrix<T>* operator+(const Matrix<T>& other) const override {
        const MatrixDense<T>* otherDense = dynamic_cast<const MatrixDense<T>*>(&other);
        if (!otherDense) {
            throw std::invalid_argument("Incompatible matrix types for addition.");
        }

        if(_rows != otherDense->_rows || _cols != otherDense->_cols) {
            throw std::invalid_argument("Matrix dimensions must be equal for addition.");
        }

        MatrixDense<T>* result = new MatrixDense<T>(_rows, _cols);
        for(size_t i = 0; i < _rows; ++i) {
            for(size_t j = 0; j < _cols; ++j){
                (*result)(i, j) = (*this)(i, j) + (*otherDense)(i, j);
            }
        }

        return result;
    }

    Matrix<T>* operator-(const Matrix<T>& other) const override {
        const MatrixDense<T>* otherDense = dynamic_cast<const MatrixDense<T>*>(&other);
        if (!otherDense) {
            throw std::invalid_argument("Incompatible matrix types for subtraction.");
        }

        if(_rows != otherDense->_rows || _cols != otherDense->_cols) {
            throw std::invalid_argument("Matrix dimensions must be equal for subtraction.");
        }

        MatrixDense<T>* result = new MatrixDense<T>(_rows, _cols);
        for(size_t i = 0; i < _rows; ++i) {
            for(size_t j = 0; j < _cols; ++j){
                (*result)(i, j) = (*this)(i, j) - (*otherDense)(i, j);
            }
        }

        return result;
    }

    Matrix<T>* operator*(const Matrix<T>& other) const override {
        const MatrixDense<T>* otherDense = dynamic_cast<const MatrixDense<T>*>(&other);
        if (!otherDense) {
            throw std::invalid_argument("Incompatible matrix types for multiplication.");
        }

        if(_cols != otherDense->_rows) {
            throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
        }

        MatrixDense<T>* result = new MatrixDense<T>(_rows, otherDense->_cols);
        for(size_t i = 0; i < _rows; ++i) {
            for(size_t j = 0; j < otherDense->_cols; ++j){
                T sum = 0;
                for(size_t k = 0; k < _cols; ++k){
                    sum += (*this)(i, k) * (*otherDense)(k, j);
                }
                (*result)(i, j) = sum;
            }
        }

        return result;
    }

    Matrix<T>* elementWiseMultiplication(const Matrix<T>& other) const override {
        const MatrixDense<T>* otherDense = dynamic_cast<const MatrixDense<T>*>(&other);
        if (!otherDense) {
            throw std::invalid_argument("Incompatible matrix types for element-wise multiplication.");
        }

        if(_rows != otherDense->_rows || _cols != otherDense->_cols) {
            throw std::invalid_argument("Matrix dimensions must be equal for element-wise multiplication.");
        }

        MatrixDense<T>* result = new MatrixDense<T>(_rows, _cols);
        for(size_t i = 0; i < _rows; ++i) {
            for(size_t j = 0; j < _cols; ++j){
                (*result)(i, j) = (*this)(i, j) * (*otherDense)(i, j);
            }
        }

        return result;
    }

    Matrix<T>* transpose() const override {
        MatrixDense<T>* result = new MatrixDense<T>(_cols, _rows);
        for(size_t i = 0; i < _rows; ++i){
            for(size_t j = 0; j < _cols; ++j){
                (*result)(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    void importFromFile(const std::string& filename) override {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        std::string matrixType;
        file >> matrixType;
        if(matrixType != "MatrixDense"){
            throw std::runtime_error("Invalid matrix type in file: " + filename);
        }

        unsigned rows, cols;
        file >> rows >> cols;
        if (file.fail() || rows == 0 || cols == 0) {
            throw std::runtime_error("Invalid matrix dimensions in file: " + filename);
        }
        delete[] _data;
        _rows = rows;
        _cols = cols;
        _data = new T[_rows * _cols];

        for (unsigned i = 0; i < _rows; ++i) {
            for (unsigned j = 0; j < _cols; ++j) {
                if (!(file >> _data[i * _cols + j])) {
                    throw std::runtime_error("Error reading matrix data from file: " + filename);
                }
            }
        }

        file.close();
    }

    void exportToFile(const std::string& filename) const override {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        file << "MatrixDense\n";
        file << _rows << " " << _cols << "\n";

        for (unsigned i = 0; i < _rows; ++i) {
            for (unsigned j = 0; j < _cols; ++j) {
                file << (*this)(i, j) << " ";
            }
            file << "\n";
        }

        file.close();
    }

    void print() const override {
        for (unsigned i = 0; i < _rows; ++i) {
            for (unsigned j = 0; j < _cols; ++j) {
                std::cout << std::setw(10) << (*this)(i, j) << " "; // Use setw for alignment
            }
            std::cout << std::endl;
        }
    }
};


