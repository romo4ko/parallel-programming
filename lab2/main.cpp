#include <iostream>
#include "MatrixDense.cpp"
#include "MatrixDiagonal.cpp"
#include "MatrixBlock.cpp"

int main() {
    try {
        MatrixDense<double> matrix1(3, 4);
        for (unsigned i = 0; i < 3; ++i) {
            for (unsigned j = 0; j < 4; ++j) {
                matrix1(i, j) = i * 4 + j;
            }
        }

        matrix1.exportToFile("matrix1.txt");

        MatrixDense<double> matrix2(0, 0);
        matrix2.importFromFile("matrix1.txt");

        matrix2.print();

        MatrixDense<double> matrix3 = matrix1;
        matrix3.print();

        MatrixDense<double> matrix4(2, 2);
        matrix4 = matrix1;

        MatrixDiagonal<int> diag1(3);
        diag1(0, 0) = 1;
        diag1(1, 1) = 5;
        diag1(2, 2) = 9;

        diag1.print();

        diag1.exportToFile("diag1.txt");

        MatrixDiagonal<int> diag2(0);
        diag2.importFromFile("diag1.txt");

        diag2.print();

                MatrixDiagonal<int> diag3 = diag1;
                diag3.print();
                MatrixDiagonal<int> diag4(2);
                diag4 = diag1;
                diag4.print();
        MatrixDiagonal<int> diag5(3);
        diag5(0,0) = 2;
        diag5(1,1) = 3;
        diag5(2,2) = 4;

        Matrix<int>* res = diag1 * diag5;

        res->print();

        delete res;
        
        MatrixBlock<int> blockMatrix(2, 2, 3);

        MatrixDense<int>* block1 = new MatrixDense<int>(3, 3);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                (*block1)(i, j) = i * 3 + j + 1;
            }
        }

        MatrixDense<int>* block2 = new MatrixDense<int>(3,3);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                (*block2)(i, j) = i + j + 1;
            }
        }
        blockMatrix.setBlock(0, 0, block1);
        blockMatrix.setBlock(1,1, block2);

        blockMatrix.print();
        blockMatrix.exportToFile("block_matrix.txt");

        MatrixBlock<int> blockMatrixFromFile(0,0,0);
        blockMatrixFromFile.importFromFile("block_matrix.txt");
        blockMatrixFromFile.print();

        MatrixBlock<int> blockMatrix_copy = blockMatrix;
        blockMatrix_copy.print();
        MatrixBlock<int> blockMatrix_copy2(1,1,1);
        blockMatrix_copy2 = blockMatrix;
        blockMatrix_copy2.print();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
