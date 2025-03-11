#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <iomanip>

// Шаблонный класс MatrixDense, наследующий от Matrix
template <typename T = double>
class MatrixDense : public Matrix<T> {
private:
    unsigned _rows, _cols; // Количество строк и столбцов
    T* _data; // Указатель на массив данных

public:
    // Конструктор, инициализирующий матрицу заданного размера
    MatrixDense(unsigned rows, unsigned cols) : _rows(rows), _cols(cols) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("Matrix dimensions must be positive."); // Проверка на положительные размеры
        }
        _data = new T[_rows * _cols]; // Выделение памяти под данные
        for (size_t i = 0; i < _rows * _cols; ++i)
            _data[i] = 0; // Инициализация нулями
    }

    // Конструктор копирования
    MatrixDense(const MatrixDense& other) : _rows(other._rows), _cols(other._cols) {
        _data = new T[_rows * _cols]; // Выделение памяти под данные
        for(size_t i = 0; i < _rows * _cols; ++i){
            _data[i] = other._data[i]; // Копирование данных
        }
    }

    // Оператор присваивания
    MatrixDense& operator=(const MatrixDense& other) {
        if (this == &other) return *this; // Проверка на самоприсваивание
        delete[] _data; // Освобождение старой памяти
        _rows = other._rows;
        _cols = other._cols;
        _data = new T[_rows * _cols]; // Выделение новой памяти
        for(size_t i = 0; i < _rows * _cols; ++i){
            _data[i] = other._data[i]; // Копирование данных
        }
        return *this;
    }

    // Деструктор
    ~MatrixDense() override {
        delete[] _data; // Освобождение памяти
    }

    // Оператор доступа к элементу (неконстантный)
    T& operator()(unsigned row, unsigned col) override {
        if (row >= _rows || col >= _cols) {
            throw std::out_of_range("Matrix index out of bounds."); // Проверка на выход за границы
        }
        return _data[row * _cols + col]; // Возврат ссылки на элемент
    }

    // Оператор доступа к элементу (константный)
    const T& operator()(unsigned row, unsigned col) const override {
        if (row >= _rows || col >= _cols) {
            throw std::out_of_range("Matrix index out of bounds."); // Проверка на выход за границы
        }
        return _data[row * _cols + col]; // Возврат константной ссылки на элемент
    }

    // Оператор сложения матриц
    Matrix<T>* operator+(const Matrix<T>& other) const override {
        const MatrixDense<T>* otherDense = dynamic_cast<const MatrixDense<T>*>(&other);
        if (!otherDense) {
            throw std::invalid_argument("Incompatible matrix types for addition."); // Проверка на совместимость типов
        }

        if(_rows != otherDense->_rows || _cols != otherDense->_cols) {
            throw std::invalid_argument("Matrix dimensions must be equal for addition."); // Проверка на совпадение размеров
        }

        MatrixDense<T>* result = new MatrixDense<T>(_rows, _cols); // Создание результирующей матрицы
        for(size_t i = 0; i < _rows; ++i) {
            for(size_t j = 0; j < _cols; ++j){
                (*result)(i, j) = (*this)(i, j) + (*otherDense)(i, j); // Сложение элементов
            }
        }

        return result; // Возврат результирующей матрицы
    }

    // Оператор вычитания матриц
    Matrix<T>* operator-(const Matrix<T>& other) const override {
        const MatrixDense<T>* otherDense = dynamic_cast<const MatrixDense<T>*>(&other);
        if (!otherDense) {
            throw std::invalid_argument("Incompatible matrix types for subtraction."); // Проверка на совместимость типов
        }

        if(_rows != otherDense->_rows || _cols != otherDense->_cols) {
            throw std::invalid_argument("Matrix dimensions must be equal for subtraction."); // Проверка на совпадение размеров
        }

        MatrixDense<T>* result = new MatrixDense<T>(_rows, _cols); // Создание результирующей матрицы
        for(size_t i = 0; i < _rows; ++i) {
            for(size_t j = 0; j < _cols; ++j){
                (*result)(i, j) = (*this)(i, j) - (*otherDense)(i, j); // Вычитание элементов
            }
        }

        return result; // Возврат результирующей матрицы
    }

    // Оператор умножения матриц
    Matrix<T>* operator*(const Matrix<T>& other) const override {
        const MatrixDense<T>* otherDense = dynamic_cast<const MatrixDense<T>*>(&other);
        if (!otherDense) {
            throw std::invalid_argument("Incompatible matrix types for multiplication."); // Проверка на совместимость типов
        }

        if(_cols != otherDense->_rows) {
            throw std::invalid_argument("Incompatible matrix dimensions for multiplication."); // Проверка на совместимость размеров
        }

        MatrixDense<T>* result = new MatrixDense<T>(_rows, otherDense->_cols); // Создание результирующей матрицы
        for(size_t i = 0; i < _rows; ++i) {
            for(size_t j = 0; j < otherDense->_cols; ++j){
                T sum = 0;
                for(size_t k = 0; k < _cols; ++k){
                    sum += (*this)(i, k) * (*otherDense)(k, j); // Умножение и суммирование элементов
                }
                (*result)(i, j) = sum; // Запись результата
            }
        }

        return result; // Возврат результирующей матрицы
    }

    // Оператор поэлементного умножения матриц
    Matrix<T>* elementWiseMultiplication(const Matrix<T>& other) const override {
        const MatrixDense<T>* otherDense = dynamic_cast<const MatrixDense<T>*>(&other);
        if (!otherDense) {
            throw std::invalid_argument("Incompatible matrix types for element-wise multiplication."); // Проверка на совместимость типов
        }

        if(_rows != otherDense->_rows || _cols != otherDense->_cols) {
            throw std::invalid_argument("Matrix dimensions must be equal for element-wise multiplication."); // Проверка на совпадение размеров
        }

        MatrixDense<T>* result = new MatrixDense<T>(_rows, _cols); // Создание результирующей матрицы
        for(size_t i = 0; i < _rows; ++i) {
            for(size_t j = 0; j < _cols; ++j){
                (*result)(i, j) = (*this)(i, j) * (*otherDense)(i, j); // Поэлементное умножение
            }
        }

        return result; // Возврат результирующей матрицы
    }

    // Метод транспонирования матрицы
    Matrix<T>* transpose() const override {
        MatrixDense<T>* result = new MatrixDense<T>(_cols, _rows); // Создание результирующей матрицы
        for(size_t i = 0; i < _rows; ++i){
            for(size_t j = 0; j < _cols; ++j){
                (*result)(j, i) = (*this)(i, j); // Перестановка элементов
            }
        }
        return result; // Возврат результирующей матрицы
    }

    // Метод импорта матрицы из файла
    void importFromFile(const std::string& filename) override {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename); // Проверка на успешное открытие файла
        }

        std::string matrixType;
        file >> matrixType;
        if(matrixType != "MatrixDense"){
            throw std::runtime_error("Invalid matrix type in file: " + filename); // Проверка на корректный тип матрицы
        }

        unsigned rows, cols;
        file >> rows >> cols;
        if (file.fail() || rows == 0 || cols == 0) {
            throw std::runtime_error("Invalid matrix dimensions in file: " + filename); // Проверка на корректные размеры матрицы
        }
        delete[] _data; // Освобождение старой памяти
        _rows = rows;
        _cols = cols;
        _data = new T[_rows * _cols]; // Выделение новой памяти

        for (unsigned i = 0; i < _rows; ++i) {
            for (unsigned j = 0; j < _cols; ++j) {
                if (!(file >> _data[i * _cols + j])) {
                    throw std::runtime_error("Error reading matrix data from file: " + filename); // Проверка на успешное чтение данных
                }
            }
        }

        file.close(); // Закрытие файла
    }

    // Метод экспорта матрицы в файл
    void exportToFile(const std::string& filename) const override {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename); // Проверка на успешное открытие файла
        }

        file << "MatrixDense\n";
        file << _rows << " " << _cols << "\n";

        for (unsigned i = 0; i < _rows; ++i) {
            for (unsigned j = 0; j < _cols; ++j) {
                file << (*this)(i, j) << " "; // Запись данных в файл
            }
            file << "\n";
        }

        file.close(); // Закрытие файла
    }

    // Метод печати матрицы на экран
    void print() const override {
        for (unsigned i = 0; i < _rows; ++i) {
            for (unsigned j = 0; j < _cols; ++j) {
                std::cout << std::setw(10) << (*this)(i, j) << " "; // Использование setw для выравнивания
            }
            std::cout << std::endl;
        }
    }
};


