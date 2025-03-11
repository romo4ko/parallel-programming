#include "Matrix.h"
#include "MatrixDense.cpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <iomanip>
#include <map>

template <typename T = double>
class MatrixBlock : public Matrix<T> {
private:
    unsigned _numBlockRows;    // Количество строк блоков
    unsigned _numBlockCols;    // Количество столбцов блоков
    unsigned _blockDim;   // Размер блока (предполагается, что блоки квадратные)
    std::vector<MatrixDense<T>*> _blockData; // Вектор указателей на блоки

public:
    // Конструктор, инициализирующий блочную матрицу заданного размера
    MatrixBlock(unsigned numBlockRows, unsigned numBlockCols, unsigned blockDim)
        : _numBlockRows(numBlockRows), _numBlockCols(numBlockCols), _blockDim(blockDim) {
            _blockData.resize(numBlockRows * numBlockCols, nullptr); // Инициализация вектора блоков
        }

    // Конструктор копирования
    MatrixBlock(const MatrixBlock& other) 
        : _numBlockRows(other._numBlockRows), _numBlockCols(other._numBlockCols), _blockDim(other._blockDim) {
        _blockData.resize(_numBlockRows * _numBlockCols, nullptr);
        std::map<MatrixDense<T>*, MatrixDense<T>*> blockMap; // Карта для отслеживания уже скопированных блоков
        for(size_t i = 0; i < _numBlockRows * _numBlockCols; ++i){
            if(other._blockData[i] != nullptr) {
                if(blockMap.find(other._blockData[i]) == blockMap.end()){
                    _blockData[i] = new MatrixDense<T>(*other._blockData[i]); // Копирование блока
                    blockMap[other._blockData[i]] = _blockData[i];
                } else {
                    _blockData[i] = blockMap[other._blockData[i]]; // Использование уже скопированного блока
                }
            }
        }
    }

    // Оператор присваивания
    MatrixBlock& operator=(const MatrixBlock& other) {
        if(this == &other) return *this; // Проверка на самоприсваивание
        for(size_t i = 0; i < _numBlockRows * _numBlockCols; ++i){
            if(_blockData[i] != nullptr){
                delete _blockData[i]; // Освобождение памяти
            }
        }
        _blockData.clear();
        _numBlockRows = other._numBlockRows;
        _numBlockCols = other._numBlockCols;
        _blockDim = other._blockDim;
        _blockData.resize(_numBlockRows * _numBlockCols);
        std::map<MatrixDense<T>*, MatrixDense<T>*> blockMap; // Карта для отслеживания уже скопированных блоков
        for(size_t i = 0; i < _numBlockRows * _numBlockCols; ++i){
            if(other._blockData[i] != nullptr) {
                if(blockMap.find(other._blockData[i]) == blockMap.end()){
                    _blockData[i] = new MatrixDense<T>(*other._blockData[i]); // Копирование блока
                    blockMap[other._blockData[i]] = _blockData[i];
                } else {
                    _blockData[i] = blockMap[other._blockData[i]]; // Использование уже скопированного блока
                }
            }
        }
        return *this;
    }

    // Деструктор
    ~MatrixBlock() override {
        for (auto block : _blockData) {
            if (block != nullptr) {
                delete block; // Освобождение памяти
            }
        }
    }

    // Оператор доступа к элементу (неконстантный)
    T& operator()(unsigned row, unsigned col) override {
        unsigned blockRow = row / _blockDim;
        unsigned blockCol = col / _blockDim;
        unsigned inBlockRow = row % _blockDim;
        unsigned inBlockCol = col % _blockDim;

        if (blockRow >= _numBlockRows || blockCol >= _numBlockCols || row >= _numBlockRows * _blockDim || col >= _numBlockCols * _blockDim) {
            throw std::out_of_range("Matrix index out of bounds."); // Проверка на выход за границы
        }
        if (_blockData[blockRow * _numBlockCols + blockCol] == nullptr){
            static T zero = 0;
            return zero; // Возврат ссылки на статическую переменную 0
        }
        return (*_blockData[blockRow * _numBlockCols + blockCol])(inBlockRow, inBlockCol); // Возврат ссылки на элемент блока
    }

    // Оператор доступа к элементу (константный)
    const T& operator()(unsigned row, unsigned col) const override {
        unsigned blockRow = row / _blockDim;
        unsigned blockCol = col / _blockDim;
        unsigned inBlockRow = row % _blockDim;
        unsigned inBlockCol = col % _blockDim;

        if (blockRow >= _numBlockRows || blockCol >= _numBlockCols || row >= _numBlockRows * _blockDim || col >= _numBlockCols * _blockDim) {
            throw std::out_of_range("Matrix index out of bounds."); // Проверка на выход за границы
        }

        if (_blockData[blockRow * _numBlockCols + blockCol] == nullptr){
            static T zero = 0;
            return zero; // Возврат ссылки на статическую переменную 0
        }
        return (*_blockData[blockRow * _numBlockCols + blockCol])(inBlockRow, inBlockCol); // Возврат ссылки на элемент блока
    }

    // Метод установки блока
    void setBlock(unsigned blockRow, unsigned blockCol, MatrixDense<T>* block) {
        if (blockRow >= _numBlockRows || blockCol >= _numBlockCols || block->getRows() != _blockDim || block->getCols() != _blockDim) {
            throw std::invalid_argument("Invalid block dimensions or position."); // Проверка на корректность размеров и позиции блока
        }
        if(_blockData[blockRow * _numBlockCols + blockCol] != nullptr){
            delete _blockData[blockRow * _numBlockCols + blockCol]; // Освобождение памяти старого блока
        }
        _blockData[blockRow * _numBlockCols + blockCol] = block; // Установка нового блока
    }

    // Метод получения количества строк
    unsigned getRows() const { return _numBlockRows * _blockDim; }

    // Метод получения количества столбцов
    unsigned getCols() const { return _numBlockCols * _blockDim; }

    // Оператор сложения матриц
    Matrix<T>* operator+(const Matrix<T>& other) const override {
        const MatrixBlock<T>* otherBlock = dynamic_cast<const MatrixBlock<T>*>(&other);
        if (!otherBlock) {
            throw std::invalid_argument("Incompatible matrix types for addition."); // Проверка на совместимость типов
        }

        if (_numBlockRows != otherBlock->_numBlockRows || _numBlockCols != otherBlock->_numBlockCols || _blockDim != otherBlock->_blockDim) {
            throw std::invalid_argument("Matrix block dimensions must be equal for addition."); // Проверка на совпадение размеров блоков
        }

        MatrixBlock<T>* result = new MatrixBlock<T>(_numBlockRows, _numBlockCols, _blockDim); // Создание результирующей матрицы

        for (size_t i = 0; i < _numBlockRows; ++i) {
            for (size_t j = 0; j < _numBlockCols; ++j) {
                if (_blockData[i * _numBlockCols + j] != nullptr && otherBlock->_blockData[i * _numBlockCols + j] != nullptr) {
                    result->setBlock(i, j, *(_blockData[i * _numBlockCols + j]->operator+(*(otherBlock->_blockData[i * _numBlockCols + j])))); // Сложение блоков
                } 
            }
        }

        return result; // Возврат результирующей матрицы
    }

    // Оператор вычитания матриц
    Matrix<T>* operator-(const Matrix<T>& other) const override {
        const MatrixBlock<T>* otherBlock = dynamic_cast<const MatrixBlock<T>*>(&other);
        if (!otherBlock) {
            throw std::invalid_argument("Incompatible matrix types for subtraction."); // Проверка на совместимость типов
        }

        if (_numBlockRows != otherBlock->_numBlockRows || _numBlockCols != otherBlock->_numBlockCols || _blockDim != otherBlock->_blockDim) {
            throw std::invalid_argument("Matrix block dimensions must be equal for subtraction."); // Проверка на совпадение размеров блоков
        }

        MatrixBlock<T>* result = new MatrixBlock<T>(_numBlockRows, _numBlockCols, _blockDim); // Создание результирующей матрицы

        for (size_t i = 0; i < _numBlockRows; ++i) {
            for (size_t j = 0; j < _numBlockCols; ++j) {
                if (_blockData[i * _numBlockCols + j] != nullptr && otherBlock->_blockData[i * _numBlockCols + j] != nullptr) {
                    result->setBlock(i, j, *(_blockData[i * _numBlockCols + j]->operator-(*(otherBlock->_blockData[i * _numBlockCols + j])))); // Вычитание блоков
                } 
            }
        }

        return result; // Возврат результирующей матрицы
    }

    // Оператор умножения матриц
    Matrix<T>* operator*(const Matrix<T>& other) const override {
        const MatrixBlock<T>* otherBlock = dynamic_cast<const MatrixBlock<T>*>(&other);
        if (!otherBlock) {
            throw std::invalid_argument("Incompatible matrix types for multiplication."); // Проверка на совместимость типов
        }

        if (_numBlockCols != otherBlock->_numBlockRows) {
            throw std::invalid_argument("Incompatible matrix dimensions for multiplication."); // Проверка на совместимость размеров
        }

        MatrixBlock<T>* result = new MatrixBlock<T>(_numBlockRows, otherBlock->_numBlockCols, _blockDim); // Создание результирующей матрицы

        for (size_t i = 0; i < _numBlockRows; ++i) {
            for (size_t j = 0; j < otherBlock->_numBlockCols; ++j) {
                MatrixDense<T>* blockSum = new MatrixDense<T>(_blockDim, _blockDim); // Создание блока для суммы
                for (size_t k = 0; k < _numBlockCols; ++k) {
                    if (_blockData[i * _numBlockCols + k] != nullptr && otherBlock->_blockData[k * otherBlock->_numBlockCols + j] != nullptr) {
                        *blockSum = *blockSum + *(_blockData[i * _numBlockCols + k]->operator*(*(otherBlock->_blockData[k * otherBlock->_numBlockCols + j]))); // Умножение и суммирование блоков
                    }
                }
                result->setBlock(i, j, blockSum); // Установка блока суммы
            }
        }

        return result; // Возврат результирующей матрицы
    }

    // Оператор поэлементного умножения матриц
    Matrix<T>* elementWiseMultiplication(const Matrix<T>& other) const override {
        const MatrixBlock<T>* otherBlock = dynamic_cast<const MatrixBlock<T>*>(&other);
        if (!otherBlock) {
            throw std::invalid_argument("Incompatible matrix types for element-wise multiplication."); // Проверка на совместимость типов
        }

        if (_numBlockRows != otherBlock->_numBlockRows || _numBlockCols != otherBlock->_numBlockCols || _blockDim != otherBlock->_blockDim) {
            throw std::invalid_argument("Matrix block dimensions must be equal for element-wise multiplication."); // Проверка на совпадение размеров блоков
        }

        MatrixBlock<T>* result = new MatrixBlock<T>(_numBlockRows, _numBlockCols, _blockDim); // Создание результирующей матрицы

        for (size_t i = 0; i < _numBlockRows; ++i) {
            for (size_t j = 0; j < _numBlockCols; ++j) {
                if (_blockData[i * _numBlockCols + j] != nullptr && otherBlock->_blockData[i * _numBlockCols + j] != nullptr) {
                    result->setBlock(i, j, *(_blockData[i * _numBlockCols + j]->elementWiseMultiplication(*(otherBlock->_blockData[i * _numBlockCols + j])))); // Поэлементное умножение блоков
                }
            }
        }

        return result; // Возврат результирующей матрицы
    }

    // Метод транспонирования матрицы
    Matrix<T>* transpose() const override {
        MatrixBlock<T>* result = new MatrixBlock<T>(_numBlockCols, _numBlockRows, _blockDim); // Создание результирующей матрицы

        for (size_t i = 0; i < _numBlockRows; ++i) {
            for (size_t j = 0; j < _numBlockCols; ++j) {
                if (_blockData[i * _numBlockCols + j] != nullptr) {
                    result->setBlock(j, i, _blockData[i * _numBlockCols + j]->transpose()); // Транспонирование блоков
                }
            }
        }

        return result; // Возврат результирующей матрицы
    }

    // Метод импорта матрицы из файла
    void importFromFile(const std::string& filename) override{
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename); // Проверка на успешное открытие файла
        }

        std::string matrixType;
        file >> matrixType;
        if(matrixType != "MatrixBlock"){
            throw std::runtime_error("Invalid matrix type in file: " + filename); // Проверка на корректный тип матрицы
        }

        unsigned numBlockRows, numBlockCols, blockDim;
        file >> numBlockRows >> numBlockCols >> blockDim;
        if (file.fail() || numBlockRows == 0 || numBlockCols == 0 || blockDim == 0) {
            throw std::runtime_error("Invalid matrix dimensions in file: " + filename); // Проверка на корректные размеры матрицы
        }
        for(size_t i = 0; i < _numBlockRows * _numBlockCols; ++i){
            if(_blockData[i] != nullptr) delete _blockData[i]; // Освобождение памяти старых блоков
        }
        _blockData.clear();
        _numBlockRows = numBlockRows;
        _numBlockCols = numBlockCols;
        _blockDim = blockDim;
        _blockData.resize(_numBlockRows * _numBlockCols);
        for(size_t i = 0; i < _numBlockRows; ++i){
            for(size_t j = 0; j < _numBlockCols; ++j){
                std::string blockType;
                file >> blockType;
                if(blockType == "nullptr"){
                    _blockData[i * _numBlockCols + j] = nullptr; // Установка пустого блока
                }else if(blockType == "MatrixDense"){
                    _blockData[i * _numBlockCols + j] = new MatrixDense<T>(blockDim, blockDim); // Создание нового блока
                    for(size_t k = 0; k < blockDim; ++k){
                        for(size_t l = 0; l < blockDim; ++l){
                            file >> (*_blockData[i * _numBlockCols + j])(k, l); // Чтение данных блока из файла
                        }
                    }
                }else{
                    throw std::runtime_error("Invalid block type in file"); // Ошибка при некорректном типе блока
                }
            }
        }
        file.close(); // Закрытие файла
    }

    // Метод экспорта матрицы в файл
    void exportToFile(const std::string& filename) const override{
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename); // Проверка на успешное открытие файла
        }
        file << "MatrixBlock\n";
        file << _numBlockRows << " " << _numBlockCols << " " << _blockDim << "\n";
        for(size_t i = 0; i < _numBlockRows; ++i){
            for(size_t j = 0; j < _numBlockCols; ++j){
                if(_blockData[i * _numBlockCols + j] == nullptr){
                    file << "nullptr\n"; // Запись пустого блока
                }else{
                    file << "MatrixDense\n";
                    for(size_t k = 0; k < _blockDim; ++k){
                        for(size_t l = 0; l < _blockDim; ++l){
                            file << (*_blockData[i * _numBlockCols + j])(k, l) << " "; // Запись данных блока в файл
                        }
                    file << "\n";
                    }
                }
            }
        }
        file.close(); // Закрытие файла
    }

    // Метод печати матрицы на экран
    void print() const override{
        for (unsigned i = 0; i < getRows(); ++i) {
            for (unsigned j = 0; j < getCols(); ++j) {
                std::cout << std::setw(10) << (*this)(i, j) << " "; // Печать элемента матрицы
            }
        std::cout << std::endl;
        }
    }
};