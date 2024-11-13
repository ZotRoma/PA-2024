#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <Eigen/Dense> 
#include <iomanip> 
#include <chrono> 

using namespace Eigen;
using namespace std;


// Операция GAXPY для обновления суммы
double gaxpy(const VectorXd& L_j, const VectorXd& L_k, int count) {
    return L_j.head(count).dot(L_k.head(count));  // Скалярное произведение (векторизация)
}

// Функция разложения Холецкого в порядке JKI с использованием GAXPY
void cholesky_jki(MatrixXd& A, MatrixXd& L, int n) {
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k <= j; ++k) {
            double sum = gaxpy(L.row(j), L.row(k), k);  // Векторизированный GAXPY для суммы произведений

            if (j == k) {
                L(j, k) = sqrt(A(j, j) - sum);
            } else {
                L(j, k) = (A(j, k) - sum) / L(k, k);
            }
        }
    }
}


void cholesky_jik(MatrixXd& A, MatrixXd& L, int n)
{
    for (int j = 0; j < n; ++j)
    {
        Eigen::VectorXd s;
        if (j > 0)
        {
            // s = A[j:n-1, j] - L[j:n-1, 0:j-1] * L[j, 0:j-1].transpose()
            s = A.block(j, j, n - j, 1);
            s.noalias() -= L.block(j, 0, n - j, j) * L.row(j).head(j).transpose();
        }
        else
        {
            // s = A[j:n-1, j]
            s = A.block(j, j, n - j, 1);
        }
        L(j, j) = sqrt(s(0));
        if (n - j - 1 > 0)
        {
            L.block(j + 1, j, n - j - 1, 1) = s.segment(1, n - j - 1) / L(j, j);
        }
    }
}

void printMatrix(const MatrixXd& M, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << setw(10) << M(i, j) << " ";
        }
        cout << endl;
    }
}

MatrixXd createRandomSimetricMatrix(int n){
    MatrixXd matrix(n,n);
    srand(time(0));

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            int value = rand() % 10;
            matrix(i,j) = value;
            matrix(j,i) = value;
        }
    }
    return matrix;
}
MatrixXd createZeroMatrix(int n){
    MatrixXd matrix(n,n);
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            matrix(i,j) = 0;
            matrix(j,i) = 0;
        }
    }
    return matrix;
}

int main() {
    int n = 10;
    // Создание двумерной матрицы для ввода данных с использованием Eigen
    MatrixXd A = createRandomSimetricMatrix(n);
    A = A * A.transpose();
    printMatrix(A, n);

    // Создание матрицы для результата разложения
    MatrixXd L_jki = MatrixXd::Zero(n, n);
    MatrixXd L_jik = MatrixXd::Zero(n, n);

    auto start_jki = std::chrono::high_resolution_clock::now();
    cholesky_jki(A, L_jki, n);
    auto end_jki = std::chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> duration_jki = end_jki - start_jki;

    auto start_jik = std::chrono::high_resolution_clock::now();
    cholesky_jik(A, L_jik, n);
    auto end_jik = std::chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> duration_jik = end_jik - start_jik;


    cout << "Lead time: " << duration_jki.count() << " ms" << std::endl;
    cout << "JKI: Result of Cholesky decomposition (matrix L):" << endl;
    printMatrix(L_jki, n);
    cout << endl << "Check (multiplying the result by the transposed matrix):" << endl;
    MatrixXd result_jki = L_jki*L_jki.transpose();
    printMatrix(result_jki, n);


    
    cout << endl << "Lead time: " << duration_jik.count() << " ms" << std::endl;
    cout << "JIK: Result of Cholesky decomposition (matrix L):" << endl;
    printMatrix(L_jik, n);
    cout << endl << "Check (multiplying the result by the transposed matrix):" << endl;
    MatrixXd result_jik = L_jik*L_jik.transpose();
    printMatrix(result_jik, n);

    return 0;
}
