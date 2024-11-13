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
    #pragma omp parallel
        {
            int j,k;
            double sum;
            #pragma omp parallel for private(j,k,sum) shared(A,L,n) schedule(dynamic, 1)
            for (int j = 0; j < n; ++j) {                
                for (k = 0; k <= j; ++k) {
                    sum = gaxpy(L.row(j), L.row(k), k);  // Векторизированный GAXPY для суммы произведений

                    if (j == k) {
                        L(j, k) = sqrt(A(j, j) - sum);
                    } else {
                        L(j, k) = (A(j, k) - sum) / L(k, k);
                    }
            }
        }
    }
}


void cholesky_jik(MatrixXd& A, MatrixXd& L, int n)
{
    #pragma omp parallel
    {
        int j;
        #pragma omp parallel for private(j) shared(A,L,n) schedule(dynamic, 1)
        for (j = 0; j < n; ++j)
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

double error(const MatrixXd& A,const MatrixXd& B){
    return (A-B).cwiseAbs().sum();
}

int main() {
    for(int n = 10; n<1001; n*=10 ){

        cout<< "n = " << n << endl;
        MatrixXd A = createRandomSimetricMatrix(n);
        A = A * A.transpose();
        
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


        cout << "Lead time JKI: " << duration_jki.count() << " ms" << std::endl;
        MatrixXd result_jki = L_jki*L_jki.transpose();
        cout  << "Check JKI (multiplying the result by the transposed matrix):" << error(A,result_jki)<< endl;
        
        cout << "Lead time JIK: " << duration_jik.count() << " ms" << std::endl;
        MatrixXd result_jik = L_jik*L_jik.transpose();
        cout  << "Check JIK (multiplying the result by the transposed matrix):" << error(A,result_jik)<< endl<<endl;
    }

    return 0;
}
