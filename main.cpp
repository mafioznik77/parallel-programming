#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <chrono>

double* CreateMatrix(int N) {
    double* matrix = new double[N * N];
    return matrix;
}

bool mats_equals(double* matrix1, double* matrix2, int N)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (matrix1[i * N + j] != matrix2[i * N + j])
            {
                return false;
            }
        }
    }
    return true;
}

double* defaultMult(double* matrix1, double* matrix2, int N)
{
    double* tmp = CreateMatrix(N);
    double sum;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            sum = 0;
            for (int k = 0; k < N; k++)
            {
                sum += matrix1[i * N + k] * matrix2[k * N + j];
            }
            tmp[i * N + j] = sum;
        }
    }
    return tmp;
}

void RandMatrix(double* matrix1, double* matrix2, int N) {
    for (int i = 0; i < N * N; ++i) {
        matrix1[i] = (std::rand() % 10000) / 1000.0f;
        matrix2[i] = (std::rand() % 10000) / 1000.0f;
    }
}

void PrintMatrix(double* matrix, int N) {
    for (int i = 0; i < N * N; i += N) {
        for (int j = 0; j < N; j++)
            std::cout << matrix[i + j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void MultMatrix(double* A, double* B, double* C, int blockSize, int N) {
    for (int i = 0; i < blockSize; ++i)
        for (int j = 0; j < blockSize; ++j)
            for (int k = 0; k < blockSize; ++k) {
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
}

void ClearMatrix(double* C, int N) {
    for (int i = 0; i < N * N; ++i) {
        C[i] = 0;
    }
}

void Canon(double* A, double* B, double* C, int n, int q) {
    int blockSize = n / q;
    for (int i = 0; i < q; ++i) {
        for (int j = 0; j < q; ++j) {
            for (int k = 0; k < q; ++k) {
                MultMatrix(&A[(i * n + (j + i + k) % q) * blockSize], &B[(((i + j + k) % q) * n + j) * blockSize],
                    &C[(i * n + j) * blockSize], blockSize, n);
            }
        }
    }
}

int main(int argc, char** argv) {
    int size = 4;
    int q = 2;
    double* A, * B, * C, * DefaultRes;

    if (argc > 2) {
        size = atoi(argv[1]);
        q = atoi(argv[2]);
    }

    A = CreateMatrix(size);
    B = CreateMatrix(size);
    C = CreateMatrix(size);
    DefaultRes = CreateMatrix(size);

    ClearMatrix(C, size);

    RandMatrix(A, B, size);
    if (size < 5) {
        PrintMatrix(A, size);
        PrintMatrix(B, size);
    }

    auto start_default = std::chrono::high_resolution_clock::now();
    DefaultRes = defaultMult(A, B, size);
    auto end_default = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_def = end_default - start_default;

    auto start_Canon = std::chrono::high_resolution_clock::now();
    Canon(A, B, C, size, q);
    auto end_Canon = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_Canon = end_Canon - start_Canon;

    std::cout << "Default multiplication: " << duration_def.count() << std::endl;
    std::cout << "Canon algorithm: " << duration_Canon.count() << std::endl;

    if (mats_equals(DefaultRes, C, size)) {
        std::cout << "Matrices are equal" << std::endl;
    }
    else {
        std::cout << "Matrices are not equal" << std::endl;
    }

    if (size < 5)
        PrintMatrix(C, size);

    delete[] A;
    delete[] B;
    delete[] C;
    delete[] DefaultRes;

    return 0;
}