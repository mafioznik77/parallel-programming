#include <tbb/tbb.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>


double* CreateMatrix(int N) {
    double* matrix = new double[N * N];
    return matrix;
}

void PrintMatrix(double* matrix, int N) {
    for (int i = 0; i < N * N; i += N) {
        for (int j = 0; j < N; j++)
            std::cout << matrix[i + j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void RandMatrix(double* matrix1, double* matrix2, int N) {
    for (int i = 0; i < N * N; ++i) {
        matrix1[i] = (std::rand() % 10000) / 1000.0f;
        matrix2[i] = (std::rand() % 10000) / 1000.0f;
    }
}

void ClearMatrix(double* pCMatrix, int N) {
    for (int i = 0; i < N * N; ++i) {
        pCMatrix[i] = 0;
    }
}

void MultMatrix(double* pAMatrix, double* pBMatrix, double* pCMatrix, int blockSize, int N) {
    for (int i = 0; i < blockSize; ++i)
        for (int j = 0; j < blockSize; ++j)
            for (int k = 0; k < blockSize; ++k) {
                pCMatrix[i * N + j] += pAMatrix[i * N + k] * pBMatrix[k * N + j];
            }
}



void Canon(double* pAMatrix, double* pBMatrix, double* pCMatrix, int n, int q) {
    int blockSize = n / q;
    for (int i = 0; i < q; ++i) {
        for (int j = 0; j < q; ++j) {
            for (int k = 0; k < q; ++k) {
                MultMatrix(&pAMatrix[(i * n + (j + i + k) % q) * blockSize],
                    &pBMatrix[(((i + j + k) % q) * n + j) * blockSize], &pCMatrix[(i * n + j) * blockSize], blockSize, n);
            }
        }
    }
}

class TBB {
    int ThreadNum_;
    double* pAMatrix_;
    double* pBMatrix_;
    double* pCMatrix_;
    int RowIndex_;
    int ColIndex_;
    int Size_;
public:
    TBB(int ThreadNum, double* pAMatrix, double* pBMatrix, double* pCMatrix, int RowIndex, int
        ColIndex, int Size) : ThreadNum_(ThreadNum), pAMatrix_(pAMatrix), pBMatrix_(pBMatrix),
        pCMatrix_(pCMatrix), RowIndex_(RowIndex), ColIndex_(ColIndex), Size_(Size) {}
    void operator()() const {
        int GridSize = static_cast<int>(sqrt(static_cast<int>(ThreadNum_)));
        int BlockSize = Size_ / GridSize;
        for (int iter = 0; iter < GridSize; iter++) {
            for (int i = RowIndex_ * BlockSize; i < (RowIndex_ + 1) * BlockSize; i++)
                for (int j = ColIndex_ * BlockSize; j < (ColIndex_ + 1) * BlockSize; j++)
                    for (int k = iter * BlockSize; k < (iter + 1) * BlockSize; k++) {
                        pCMatrix_[i * Size_ + j] += pAMatrix_[i * Size_ + k] * pBMatrix_[k * Size_ + j];
                    }
        }
    }
};

int main(int argc, char** argv) {
     std::cout << "Chislo potokov (q) - polniy kvadrat" << std::endl;
    int size = 16;
    int q = 4;
    int test = 0;
    double* A, * B, * C, * SS, * S, * C1, * C3;
    if (argc > 2) {
        size = atoi(argv[1]);
        q = atoi(argv[2]);
    }
    A = CreateMatrix(size);
    SS = CreateMatrix(size);
    B = CreateMatrix(size);
    C = CreateMatrix(size);
    S = CreateMatrix(size);
    C1 = CreateMatrix(size);
    C3 = CreateMatrix(size);

    ClearMatrix(C, size);
    ClearMatrix(SS, size);
    ClearMatrix(S, size);
    ClearMatrix(C1, size);
    ClearMatrix(C3, size);

    RandMatrix(A, B, size);

    if (size < 17) {
        PrintMatrix(A, size);
        PrintMatrix(B, size);
    }

    tbb::tick_count time = tbb::tick_count::now();
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            for (int k = 0; k < size; ++k) {
                SS[i * size + j] += A[i * size + k] * B[k * size + j];
            }
    double time_easy = (tbb::tick_count::now() - time).seconds();

    time = tbb::tick_count::now();
    Canon(A, B, C, size, q);
    double time_sequence = (tbb::tick_count::now() - time).seconds();

    tbb::task_scheduler_init init(q);
    time = tbb::tick_count::now();
    tbb::task_group tg;
    for (int k = 0; k < sqrt(static_cast<int>(q)); k++)
        for (int j = 0; j < sqrt(static_cast<int>(q)); j++)
            tg.run(TBB(q, A, B, C3, k, j, size));
    tg.wait();
    double time_TBB = (tbb::tick_count::now() - time).seconds();


    if (size < 17) {
        PrintMatrix(SS, size);
        PrintMatrix(C, size);
        PrintMatrix(C3, size);
    }

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            if (fabs(SS[i * size + j] - C3[i * size + j]) < 0.1)
                test++;
            else
                test = 0;
        }

    std::cout << "Time easy version is  " << time_easy << std::endl;
    std::cout << "Time sequence version is " << time_sequence << std::endl;
    std::cout << "Time parallel_TBB version is = " << time_TBB << std::endl;
    std::cout << "Boost is " << time_easy / time_TBB << std::endl;
    if (test == 0)
        std::cout << "Result not right" << std::endl;
    else
        std::cout << "Result right!" << std::endl;

    delete[] A;
    delete[] B;
    delete[] C;
    delete[] S;
    delete[] SS;
    delete[] C3;
    delete[] C1;

    return 0;
}