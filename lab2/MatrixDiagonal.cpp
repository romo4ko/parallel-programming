#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <iomanip>

template <typename T = double>
class MatrixDiagonal : public Matrix<T> {
private:
    unsigned _dimension; // Размер матрицы (количество строк/столбцов)
    T* _elements;      // Массив для хранения диагональных элементов

public:
    // Конструктор, инициализирующий диагональную матрицу заданного размера
    MatrixDiagonal(unsigned dimension) : _dimension(dimension) {
        if (dimension == 0) {
            throw std::invalid_argument("Matrix size must be positive."); // Проверка на положительный размер
        }
        _elements = new T[_dimension]; // Выделение памяти под диагональные элементы
        for(size_t i = 0; i < _dimension; ++i)
            _elements[i] = 0; // Инициализация нулями
    }

    // Конструктор копирования
    MatrixDiagonal(const MatrixDiagonal& other) : _dimension(other._dimension) {
        _elements = new T[_dimension]; // Выделение памяти под диагональные элементы
        for(size_t i = 0; i < _dimension; ++i){
            _elements[i] = other._elements[i]; // Копирование данных
        }
    }

    // Оператор присваивания
    MatrixDiagonal& operator=(const MatrixDiagonal& other) {
        if (this == &other) return *this; // Проверка на самоприсваивание
        delete[] _elements; // Освобождение старой памяти
        _dimension = other._dimension;
        _elements = new T[_dimension]; // Выделение новой памяти
        for(size_t i = 0; i < _dimension; ++i){
            _elements[i] = other._elements[i]; // Копирование данных
        }
        return *this;
    }

    // Деструктор
    ~MatrixDiagonal() override {
        delete[] _elements; // Освобождение памяти
    }

    // Оператор доступа к элементу (неконстантный)
    T& operator()(unsigned row, unsigned col) override {
        if (row >= _dimension || col >= _dimension) {
            throw std::out_of_range("Matrix index out of bounds."); // Проверка на выход за границы
        }
        if (row == col) {
            return _elements[row]; // Возврат ссылки на диагональный элемент
        } else {
            static T zero = 0;
            return zero; // Возврат ссылки на статическую переменную 0
        }
    }

    // Оператор доступа к элементу (константный)
    const T& operator()(unsigned row, unsigned col) const override {
        if (row >= _dimension || col >= _dimension) {
            throw std::out_of_range("Matrix index out of bounds."); // Проверка на выход за границы
        }
        if (row == col) {
            return _elements[row]; // Возврат константной ссылки на диагональный элемент
        } else {
            static T zero = 0;
            return zero; // Возврат константной ссылки на статическую переменную 0
        }
    }

    // Оператор сложения матриц
    Matrix<T>* operator+(const Matrix<T>& other) const override {
        const MatrixDiagonal<T>* otherDiag = dynamic_cast<const MatrixDiagonal<T>*>(&other);
        if (!otherDiag) {
            throw std::invalid_argument("Incompatible matrix types for addition."); // Проверка на совместимость типов
        }

        if(_dimension != otherDiag->_dimension) {
            throw std::invalid_argument("Matrix dimensions must be equal for addition."); // Проверка на совпадение размеров
        }

        MatrixDiagonal<T>* result = new MatrixDiagonal<T>(_dimension); // Создание результирующей матрицы
        for(size_t i = 0; i < _dimension; ++i) {
            (*result)(i, i) = (*this)(i, i) + (*otherDiag)(i, i); // Сложение диагональных элементов
        }

        return result; // Возврат результирующей матрицы
    }

    // Оператор вычитания матриц
    Matrix<T>* operator-(const Matrix<T>& other) const override {
        const MatrixDiagonal<T>* otherDiag = dynamic_cast<const MatrixDiagonal<T>*>(&other);
        if (!otherDiag) {
            throw std::invalid_argument("Incompatible matrix types for subtraction."); // Проверка на совместимость типов
        }

        if(_dimension != otherDiag->_dimension) {
            throw std::invalid_argument("Matrix dimensions must be equal for subtraction."); // Проверка на совпадение размеров
        }

        MatrixDiagonal<T>* result = new MatrixDiagonal<T>(_dimension); // Создание результирующей матрицы
        for(size_t i = 0; i < _dimension; ++i) {
            (*result)(i, i) = (*this)(i, i) - (*otherDiag)(i, i); // Вычитание диагональных элементов
        }

        return result; // Возврат результирующей матрицы
    }

    // Оператор умножения матриц
    Matrix<T>* operator*(const Matrix<T>& other) const override {
        const MatrixDiagonal<T>* otherDiag = dynamic_cast<const MatrixDiagonal<T>*>(&other);
        if (!otherDiag) {
            throw std::invalid_argument("Incompatible matrix types for multiplication."); // Проверка на совместимость типов
        }

        if(_dimension != otherDiag->_dimension) {
            throw std::invalid_argument("Incompatible matrix dimensions for multiplication."); // Проверка на совместимость размеров
        }

        MatrixDiagonal<T>* result = new MatrixDiagonal<T>(_dimension); // Создание результирующей матрицы
        for(size_t i = 0; i < _dimension; ++i) {
            (*result)(i, i) = (*this)(i, i) * (*otherDiag)(i, i); // Умножение диагональных элементов
        }

        return result; // Возврат результирующей матрицы
    }

    // Оператор поэлементного умножения матриц
    Matrix<T>* elementWiseMultiplication(const Matrix<T>& other) const override {
        return this->operator*(other); // Поэлементное умножение совпадает с обычным умножением для диагональных матриц
    }

    // Метод транспонирования матрицы
    Matrix<T>* transpose() const override {
        return new MatrixDiagonal<T>(*this); // Транспонирование диагональной матрицы не изменяет её
    }

    // Метод импорта матрицы из файла
    void importFromFile(const std::string& filename) override {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename); // Проверка на успешное открытие файла
        }

        std::string matrixType;
        file >> matrixType;
        if(matrixType != "MatrixDiagonal"){
            throw std::runtime_error("Invalid matrix type in file: " + filename); // Проверка на корректный тип матрицы
        }

        unsigned dimension;
        file >> dimension;
        if (file.fail() || dimension == 0) {
            throw std::runtime_error("Invalid matrix dimensions in file: " + filename); // Проверка на корректные размеры матрицы
        }
        delete[] _elements; // Освобождение старой памяти
        _dimension = dimension;
        _elements = new T[_dimension]; // Выделение новой памяти

        for (unsigned i = 0; i < _dimension; ++i) {
            if (!(file >> _elements[i])) {
                throw std::runtime_error("Error reading matrix data from file: " + filename); // Проверка на успешное чтение данных
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

        file << "MatrixDiagonal\n";
        file << _dimension << "\n";

        for (unsigned i = 0; i < _dimension; ++i) {
            file << _elements[i] << " "; // Запись диагональных элементов в файл
        }
        file << "\n";

        file.close(); // Закрытие файла
    }

    // Метод печати матрицы на экран
    void print() const override {
        for (unsigned i = 0; i < _dimension; ++i) {
            for (unsigned j = 0; j < _dimension; ++j) {
                if (i == j) {
                    std::cout << std::setw(10) << _elements[i] << " "; // Печать диагонального элемента
                } else {
                    std::cout << std::setw(10) << 0 << " "; // Печать нуля для недиагональных элементов
                }
            }
            std::cout << std::endl;
        }
    }
};