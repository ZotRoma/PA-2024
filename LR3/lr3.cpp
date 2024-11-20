#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <chrono> 
#include <iomanip> 
#include <omp.h> 

using namespace std;

void printMatrix(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

double** createRandomSimetricMatrix(int n){
    double** matrix = new double*[n];
    for(int i = 0; i<n; i++){
        matrix[i] = new double[n];
    }
    srand(time(0));

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            int value = rand() % 10;
            matrix[i][j] = value;
            matrix[j][i] = value;
        }
    }
    return matrix;
}
double** createZeroMatrix(int n){
    double** matrix = new double*[n];
    for(int i = 0; i<n; i++){
        matrix[i] = new double[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            matrix[i][j] = 0;
            matrix[j][i] = 0;
        }
    }
    return matrix;
}
double** copyMatrix(double** m, int n){
    double** matrix = new double*[n];
    for(int i = 0; i<n; i++){
        matrix[i] = new double[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            matrix[i][j] = m[i][j];
            matrix[j][i] = m[j][i];
        }
    }
    return matrix;
}


void multiplyTranspose(double** B, double** A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            A[i][j] = 0;
            for (int k = 0; k < n; k++) {
                A[i][j] += B[i][k] * B[j][k];
            }
            A[j][i] = A[i][j];  // Заполняем симметрично
        }
    }
}


void chileskyDecompositionJIK(double** A, double** L, int n){
    omp_set_num_threads(4);
    #pragma omp parallel
    {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;

            // Вычисление диагонального элемента L[j][j]
            #pragma omp single
            {
                for (int k = 0; k < j; k++) {
                    sum += L[j][k] * L[j][k];
                }
                L[j][j] = sqrt(A[j][j] - sum);
            }

            // Вычисление недиагональных элементов столбца j
            #pragma omp for schedule(dynamic)
            for (int i = j + 1; i < n; i++) {
                sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }

            // Заполняем верхний треугольник матрицы нулями (опционально)
            #pragma omp for schedule(static)
            for (int i = 0; i < j; i++) {
                L[i][j] = 0.0;
            }
        }
    }
}

void choleskyDecompositionJKI(double** A, double** L, int n){
    omp_set_num_threads(4);
    for (int j = 0; j < n; j++) {
        // Вычисление диагонального элемента L[j][j]
        double sum = 0.0;
        #pragma omp parallel for reduction(+:sum)
        for (int i = 0; i < j; i++) {
            sum += L[j][i] * L[j][i];
        }
        L[j][j] = sqrt(A[j][j] - sum);

        // Вычисление элементов ниже диагонали в столбце j
        #pragma omp parallel for
        for (int k = j + 1; k < n; k++) {
            double sum2 = 0.0;
            for (int i = 0; i < j; i++) {
                sum2 += L[k][i] * L[j][i];
            }
            L[k][j] = (A[k][j] - sum2) / L[j][j];
        }
    } 
    for(int i = 0; i <n; i++){
        for(int j = i+1; j<n; j++){
            L[i][j] = 0;
        }
    }
}


int main(void)
{
    double** pog = new double*[2];
    pog[0] = new double[10];
    pog[1] = new double[10];

    
    for(int nn = 10, t = 0; nn<1001; nn*=10, t++){
    const int n = nn;   
    

    double** A = createRandomSimetricMatrix(n);
    double** B = createZeroMatrix(n);
    multiplyTranspose(A,B,n);
    double** L = copyMatrix(B, n);
    
    double** LCopy = copyMatrix(L, n);

    auto start_jki = std::chrono::high_resolution_clock::now();
    choleskyDecompositionJKI(B,L,n);
    auto end_jki = std::chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> duration_jki = end_jki - start_jki;
    

    auto start_jik = std::chrono::high_resolution_clock::now();
    chileskyDecompositionJIK(B,LCopy,n);
    auto end_jik = std::chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> duration_jik = end_jik - start_jik;


    double** res_jik = copyMatrix(LCopy, n);
    double** res_jik_f = createZeroMatrix(n);
    double** res_jki = copyMatrix(L, n);
    double** res_jki_f = createZeroMatrix(n);
    multiplyTranspose(res_jik,res_jik_f,n);
    multiplyTranspose(res_jki,res_jki_f,n);
    cout<<endl;

    double error_jki = 0, error_jik = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            error_jki += abs(res_jki_f[i][j] - B[i][j]);
            error_jik += abs(res_jik_f[i][j] - B[i][j]);
        }
    }
    cout<<"n = "<<nn<< ", time JIK: " << duration_jik.count() << " ms, time JKI:"  << duration_jki.count() << " ms" << std::endl;
    pog[0][t] = error_jki;
    pog[1][t] = error_jik;
    for (int i = 0; i < n; i++) {
            delete[] A[i];
            delete[] LCopy[i];
            delete[] L[i];
            delete[] B[i];
            delete[] res_jik[i];
            delete[] res_jik_f[i];
            delete[] res_jki[i];
            delete[] res_jki_f[i];

        }
        delete[] A;
        delete[] L;
        delete[] LCopy;
        delete[] B;
        delete[] res_jik;
        delete[] res_jik_f;
        delete[] res_jki;
        delete[] res_jki_f;
    }
    

    cout<<"error jki:";
    for(int i = 0; i<3; i++){
        cout<<pog[0][i]<<' ';
    }
    cout<<endl<<"error jik:";
    for(int i = 0; i<3; i++){
        cout<<pog[1][i]<<' ';
    }
    delete[] pog[0];
    delete[] pog[1];
    delete[] pog;
    return 0;
}