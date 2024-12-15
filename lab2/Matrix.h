#pragma once

#include <iostream>
#include <fstream>
#include <stdexcept>

template <typename T = double>
class Matrix {
public:

    virtual ~Matrix() = default;

    virtual T& operator()(unsigned i, unsigned j) = 0;
    virtual const T& operator()(unsigned i, unsigned j) const = 0;

    virtual Matrix<T>* operator+(const Matrix<T>& other) const = 0;
    virtual Matrix<T>* operator-(const Matrix<T>& other) const = 0;
    virtual Matrix<T>* operator*(const Matrix<T>& other) const = 0;
    virtual Matrix<T>* elementWiseMultiplication(const Matrix<T>& other) const = 0;
    virtual Matrix<T>* transpose() const = 0;

    virtual void importFromFile(const std::string& filename) = 0;
    virtual void exportToFile(const std::string& filename) const = 0;
    virtual void print() const = 0;
};