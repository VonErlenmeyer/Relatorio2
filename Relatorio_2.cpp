#include <iostream>
#include <cmath>
using namespace std;

void inverseSub(double (*matrix)[][3], double (*vector_b)[3], int n);
void gaussNoPivot(double (*matrix)[][3], double (*vector_b)[3], int n);

int main(void){

double matrix[3][3]={{1,-1,-1},{0,2,-1},{0,0,-1}};
double vector_b[3]={0,-2,-7};
int n=3;

cout << "1. Inversa, R3=0" << endl;
inverseSub(&matrix, &vector_b, n);

matrix[2][0]= -2;
cout << "1. Inversa, R3=2" << endl;
inverseSub(&matrix, &vector_b, n);

cout << "Gauss R3=2" << endl;
gaussNoPivot(&matrix, &vector_b, n);

matrix[1][1]=0;
cout << "Gauss R2=0" << endl;
gaussNoPivot(&matrix, &vector_b, n);

return 0;
}


void inverseSub(double (*matrix)[][3], double (*vector_b)[3], int n){

int j, i;
double sum=0;
double vector_x[3]={0};

vector_x[n-1] = ((*vector_b)[n-1])/((*matrix)[n-1][n-1]);

for(i=n-1; i>=0; i--){
    sum=0;
    for(j=i+1; j<=n; j++){
        sum += (*matrix)[i][j]*(vector_x[j]);
    }
    vector_x[i]=((*vector_b)[i]-sum)/(*matrix)[i][i];
}

for(i=0; i<3; i++){
    cout << vector_x[i] << endl;
}
}

void gaussNoPivot(double (*matrix)[][3], double (*vector_b)[3], int n){

int i, k, j;
double m;


for(k=0; k<n-1; k++){
    for(i=k+1; i<n; i++){
    m = (*matrix)[i][k] / (*matrix)[k][k];
        for(j=k; j<n; j++){
        (*matrix)[i][j] = (*matrix)[i][j]- m*(*matrix)[k][j];
        } 
        (*vector_b)[i] = (*vector_b)[i]-m*(*vector_b)[k];
    }
}

inverseSub(matrix, vector_b, n);
}