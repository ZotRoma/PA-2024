#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <math.h>

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
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) {
            sum += L[j][k] * L[j][k];
        }
        L[j][j] = sqrt(A[j][j] - sum);

        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i][k] * L[j][k];
            }
            L[i][j] = (A[i][j] - sum) / L[j][j];
        }
    }
    for(int i = 0; i <n; i++){
        for(int j = i+1; j<n; j++){
            L[i][j] = 0;
        }
    }
}

void choleskyDecompositionJKI(double** A, double** L, int n){
    // for(int j = 0; j<n; j++){
    //     for(int k = 0; k<j; k++){
    //         for(int i = j; i<n; i++){
                
    //             L[i][j] = L[i][j]-L[i][k]*L[j][k];
    //         }
    //     }
    //     for(int i = j; i<n; i++){
    //         L[i][j] = L[i][j]/sqrt(L[j][j]);
    //     }
    //     cout<<endl;
    //     printMatrix(L, n);
    // }
    for (int J = 0; J < n; J++) {
        for (int k = 0; k <= J; k++) {
            double sum = 0;
            for (int i = 0; i < k; i++) {
                sum += L[J][i] * L[k][i];
            }
            if (J == k) {
                L[J][k] = sqrt(A[J][J] - sum);
            } else {
                L[J][k] = (A[J][k] - sum) / L[k][k];
            }
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
    double* pog = new double[10];
    
    for(int nn = 5, t = 0; nn<6; nn*=10, t++){
    const int n = nn;   
    
    // double** matrix = createRandomSimetricMatrix(n);
    
    
    // cout << "Symmetric Matrix (" << n << "x" << n << "):" << endl;
    // printMatrix(matrix, n);
    // cout << "Copy Matrix (" << n << "x" << n << "):" << endl;
    // printMatrix(L, n);
    // cout << "Cholesky Matrix (" << n << "x" << n << "):" << endl;
    // choleskyDecompositionJKI(matrix,L,n);
    // printMatrix(L, n);

    double** A = createRandomSimetricMatrix(n);
    //printMatrix(A, n);
    //cout<<endl;
    double** B = createZeroMatrix(n);
    //printMatrix(B, n);
   // cout<<endl;
    multiplyTranspose(A,B,n);
    //printMatrix(B, n);
    //cout<<endl;
    double** L = copyMatrix(B, n);
    
    double** LCopy = copyMatrix(L, n);

    //cout  << endl << "Cholesky Matrix JIK (" << n << "x" << n << "):" << endl;
    chileskyDecompositionJIK(B,L,n);
    //printMatrix(L, n);


    //cout  << endl << "Cholesky Matrix JKI (" << n << "x" << n << "):" << endl;
    choleskyDecompositionJKI(B,LCopy,n);
    //printMatrix(LCopy, n);


    double** res = copyMatrix(LCopy, n);
    double** resT = createZeroMatrix(n);
   // cout  << endl << "matrix for checker" << endl;
    multiplyTranspose(res,resT,n);
    cout<<endl;
    //printMatrix(resT, n);

    double dd = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            dd += abs(resT[i][j] - B[i][j]);
        }
    }
        cout<<"nn iter = "<<nn<<endl;
        pog[t] = dd;
        for (int i = 0; i < n; i++) {
            delete[] A[i];
            delete[] LCopy[i];
            delete[] L[i];
            delete[] B[i];
            delete[] res[i];
            delete[] resT[i];

        }
        delete[] A;
        delete[] L;
        delete[] LCopy;
        delete[] B;
        delete[] res;
        delete[] resT;
    }
    

    cout<<"discrepancy:";
    for(int i = 0; i<3; i++){
        cout<<pog[i]<<' ';
    }
    return 0;
}