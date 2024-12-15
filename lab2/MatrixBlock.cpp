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
    unsigned _numBlockRows;    // Number of block rows
    unsigned _numBlockCols;    // Number of block columns
    unsigned _blockDim;   // Block dimension (assuming square blocks)
    std::vector<MatrixDense<T>*> _blockData;

public:
    MatrixBlock(unsigned numBlockRows, unsigned numBlockCols, unsigned blockDim)
        : _numBlockRows(numBlockRows), _numBlockCols(numBlockCols), _blockDim(blockDim) {
            _blockData.resize(numBlockRows * numBlockCols, nullptr);
        }

    MatrixBlock(const MatrixBlock& other) 
        : _numBlockRows(other._numBlockRows), _numBlockCols(other._numBlockCols), _blockDim(other._blockDim) {
        _blockData.resize(_numBlockRows * _numBlockCols, nullptr);
        std::map<MatrixDense<T>*, MatrixDense<T>*> blockMap;
        for(size_t i = 0; i < _numBlockRows * _numBlockCols; ++i){
            if(other._blockData[i] != nullptr) {
                if(blockMap.find(other._blockData[i]) == blockMap.end()){
                    _blockData[i] = new MatrixDense<T>(*other._blockData[i]);
                    blockMap[other._blockData[i]] = _blockData[i];
                } else {
                    _blockData[i] = blockMap[other._blockData[i]];
                }
            }
        }
    }

    MatrixBlock& operator=(const MatrixBlock& other) {
        if(this == &other) return *this;
        for(size_t i = 0; i < _numBlockRows * _numBlockCols; ++i){
            if(_blockData[i] != nullptr){
                delete _blockData[i];
            }
        }
        _blockData.clear();
        _numBlockRows = other._numBlockRows;
        _numBlockCols = other._numBlockCols;
        _blockDim = other._blockDim;
        _blockData.resize(_numBlockRows * _numBlockCols);
        std::map<MatrixDense<T>*, MatrixDense<T>*> blockMap;
        for(size_t i = 0; i < _numBlockRows * _numBlockCols; ++i){
            if(other._blockData[i] != nullptr) {
                if(blockMap.find(other._blockData[i]) == blockMap.end()){
                    _blockData[i] = new MatrixDense<T>(*other._blockData[i]);
                    blockMap[other._blockData[i]] = _blockData[i];
                } else {
                    _blockData[i] = blockMap[other._blockData[i]];
                }
            }
        }
        return *this;
    }

    ~MatrixBlock() override {
        for (auto block : _blockData) {
            if (block != nullptr) {
                delete block;
            }
        }
    }

    T& operator()(unsigned row, unsigned col) override {
        unsigned blockRow = row / _blockDim;
        unsigned blockCol = col / _blockDim;
        unsigned inBlockRow = row % _blockDim;
        unsigned inBlockCol = col % _blockDim;

        if (blockRow >= _numBlockRows || blockCol >= _numBlockCols || row >= _numBlockRows * _blockDim || col >= _numBlockCols * _blockDim) {
            throw std::out_of_range("Matrix index out of bounds.");
        }
        if (_blockData[blockRow * _numBlockCols + blockCol] == nullptr){
            static T zero = 0;
            return zero;
        }
        return (*_blockData[blockRow * _numBlockCols + blockCol])(inBlockRow, inBlockCol);
    }

    const T& operator()(unsigned row, unsigned col) const override {
        unsigned blockRow = row / _blockDim;
        unsigned blockCol = col / _blockDim;
        unsigned inBlockRow = row % _blockDim;
        unsigned inBlockCol = col % _blockDim;

        if (blockRow >= _numBlockRows || blockCol >= _numBlockCols || row >= _numBlockRows * _blockDim || col >= _numBlockCols * _blockDim) {
            throw std::out_of_range("Matrix index out of bounds.");
        }

        if (_blockData[blockRow * _numBlockCols + blockCol] == nullptr){
            static T zero = 0;
            return zero;
        }
        return (*_blockData[blockRow * _numBlockCols + blockCol])(inBlockRow, inBlockCol);
    }

    void setBlock(unsigned blockRow, unsigned blockCol, MatrixDense<T>* block) {
        if (blockRow >= _numBlockRows || blockCol >= _numBlockCols || block->getRows() != _blockDim || block->getCols() != _blockDim) {
            throw std::invalid_argument("Invalid block dimensions or position.");
        }
        if(_blockData[blockRow * _numBlockCols + blockCol] != nullptr){
            delete _blockData[blockRow * _numBlockCols + blockCol];
        }
        _blockData[blockRow * _numBlockCols + blockCol] = block;
    }

    unsigned getRows() const { return _numBlockRows * _blockDim; }
    unsigned getCols() const { return _numBlockCols * _blockDim; }

    Matrix<T>* operator+(const Matrix<T>& other) const override {
        const MatrixBlock<T>* otherBlock = dynamic_cast<const MatrixBlock<T>*>(&other);
        if (!otherBlock) {
            throw std::invalid_argument("Incompatible matrix types for addition.");
        }

        if (_numBlockRows != otherBlock->_numBlockRows || _numBlockCols != otherBlock->_numBlockCols || _blockDim != otherBlock->_blockDim) {
            throw std::invalid_argument("Matrix block dimensions must be equal for addition.");
        }

        MatrixBlock<T>* result = new MatrixBlock<T>(_numBlockRows, _numBlockCols, _blockDim);

        for (size_t i = 0; i < _numBlockRows; ++i) {
            for (size_t j = 0; j < _numBlockCols; ++j) {
                if (_blockData[i * _numBlockCols + j] != nullptr && otherBlock->_blockData[i * _numBlockCols + j] != nullptr) {
                    result->setBlock(i, j, *(_blockData[i * _numBlockCols + j]->operator+(*(otherBlock->_blockData[i * _numBlockCols + j])))); 
                } 
            }
        }

        return result;
    }

    Matrix<T>* operator-(const Matrix<T>& other) const override {
        const MatrixBlock<T>* otherBlock = dynamic_cast<const MatrixBlock<T>*>(&other);
        if (!otherBlock) {
            throw std::invalid_argument("Incompatible matrix types for subtraction.");
        }

        if (_numBlockRows != otherBlock->_numBlockRows || _numBlockCols != otherBlock->_numBlockCols || _blockDim != otherBlock->_blockDim) {
            throw std::invalid_argument("Matrix block dimensions must be equal for subtraction.");
        }

        MatrixBlock<T>* result = new MatrixBlock<T>(_numBlockRows, _numBlockCols, _blockDim);

        for (size_t i = 0; i < _numBlockRows; ++i) {
            for (size_t j = 0; j < _numBlockCols; ++j) {
                if (_blockData[i * _numBlockCols + j] != nullptr && otherBlock->_blockData[i * _numBlockCols + j] != nullptr) {
                    result->setBlock(i, j, *(_blockData[i * _numBlockCols + j]->operator-(*(otherBlock->_blockData[i * _numBlockCols + j])))); 
                } 
            }
        }

        return result;
    }

    Matrix<T>* operator*(const Matrix<T>& other) const override {
        const MatrixBlock<T>* otherBlock = dynamic_cast<const MatrixBlock<T>*>(&other);
        if (!otherBlock) {
            throw std::invalid_argument("Incompatible matrix types for multiplication.");
        }

        if (_numBlockCols != otherBlock->_numBlockRows) {
            throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
        }

        MatrixBlock<T>* result = new MatrixBlock<T>(_numBlockRows, otherBlock->_numBlockCols, _blockDim);

        for (size_t i = 0; i < _numBlockRows; ++i) {
            for (size_t j = 0; j < otherBlock->_numBlockCols; ++j) {
                MatrixDense<T>* blockSum = new MatrixDense<T>(_blockDim, _blockDim);
                for (size_t k = 0; k < _numBlockCols; ++k) {
                    if (_blockData[i * _numBlockCols + k] != nullptr && otherBlock->_blockData[k * otherBlock->_numBlockCols + j] != nullptr) {
                        *blockSum = *blockSum + *(_blockData[i * _numBlockCols + k]->operator*(*(otherBlock->_blockData[k * otherBlock->_numBlockCols + j]))); 
                    }
                }
                result->setBlock(i, j, blockSum);
            }
        }

        return result;
    }

    Matrix<T>* elementWiseMultiplication(const Matrix<T>& other) const override {
        const MatrixBlock<T>* otherBlock = dynamic_cast<const MatrixBlock<T>*>(&other);
        if (!otherBlock) {
            throw std::invalid_argument("Incompatible matrix types for element-wise multiplication.");
        }

        if (_numBlockRows != otherBlock->_numBlockRows || _numBlockCols != otherBlock->_numBlockCols || _blockDim != otherBlock->_blockDim) {
            throw std::invalid_argument("Matrix block dimensions must be equal for element-wise multiplication.");
        }

        MatrixBlock<T>* result = new MatrixBlock<T>(_numBlockRows, _numBlockCols, _blockDim);

        for (size_t i = 0; i < _numBlockRows; ++i) {
            for (size_t j = 0; j < _numBlockCols; ++j) {
                if (_blockData[i * _numBlockCols + j] != nullptr && otherBlock->_blockData[i * _numBlockCols + j] != nullptr) {
                    result->setBlock(i, j, *(_blockData[i * _numBlockCols + j]->elementWiseMultiplication(*(otherBlock->_blockData[i * _numBlockCols + j])))); 
                }
            }
        }

        return result;
    }

    Matrix<T>* transpose() const override {
        MatrixBlock<T>* result = new MatrixBlock<T>(_numBlockCols, _numBlockRows, _blockDim);

        for (size_t i = 0; i < _numBlockRows; ++i) {
            for (size_t j = 0; j < _numBlockCols; ++j) {
                if (_blockData[i * _numBlockCols + j] != nullptr) {
                    result->setBlock(j, i, _blockData[i * _numBlockCols + j]->transpose()); 
                }
            }
        }

        return result;
    }

    // Other methods: operations, importFromFile, exportToFile, print (see below)
    void importFromFile(const std::string& filename) override{
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        std::string matrixType;
        file >> matrixType;
        if(matrixType != "MatrixBlock"){
            throw std::runtime_error("Invalid matrix type in file: " + filename);
        }

        unsigned numBlockRows, numBlockCols, blockDim;
        file >> numBlockRows >> numBlockCols >> blockDim;
        if (file.fail() || numBlockRows == 0 || numBlockCols == 0 || blockDim == 0) {
            throw std::runtime_error("Invalid matrix dimensions in file: " + filename);
        }
        for(size_t i = 0; i < _numBlockRows * _numBlockCols; ++i){
            if(_blockData[i] != nullptr) delete _blockData[i];
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
                    _blockData[i * _numBlockCols + j] = nullptr;
                }else if(blockType == "MatrixDense"){
                    _blockData[i * _numBlockCols + j] = new MatrixDense<T>(blockDim, blockDim);
                    for(size_t k = 0; k < blockDim; ++k){
                        for(size_t l = 0; l < blockDim; ++l){
                            file >> (*_blockData[i * _numBlockCols + j])(k, l);
                        }
                    }
                }else{
                    throw std::runtime_error("Invalid block type in file");
                }
            }
        }
        file.close();
    }

    void exportToFile(const std::string& filename) const override{
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }
        file << "MatrixBlock\n";
        file << _numBlockRows << " " << _numBlockCols << " " << _blockDim << "\n";
        for(size_t i = 0; i < _numBlockRows; ++i){
            for(size_t j = 0; j < _numBlockCols; ++j){
                if(_blockData[i * _numBlockCols + j] == nullptr){
                    file << "nullptr\n";
                }else{
                    file << "MatrixDense\n";
                    for(size_t k = 0; k < _blockDim; ++k){
                        for(size_t l = 0; l < _blockDim; ++l){
                            file << (*_blockData[i * _numBlockCols + j])(k, l) << " ";
                        }
                    file << "\n";
                    }
                }
            }
        }
        file.close();
    }

    void print() const override{
        for (unsigned i = 0; i < getRows(); ++i) {
            for (unsigned j = 0; j < getCols(); ++j) {
                std::cout << std::setw(10) << (*this)(i, j) << " ";
            }
        std::cout << std::endl;
        }
    }
};