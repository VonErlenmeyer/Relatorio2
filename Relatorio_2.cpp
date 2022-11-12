#include <iostream>
using namespace std;

void inverse_sub(double (*matrix)[][3], double (*vector_b)[3], int n);

int main(void){

double matrix[3][3]={{1,-1,-1},{0,2,-1},{0,0,-1}};
double vector_b[3]={0,-2,-7};
int n=3;
inverse_sub(&matrix, &vector_b, n);

matrix[2][0]=-2;
inverse_sub(&matrix, &vector_b, n);

return 0;
}


void inverse_sub(double (*matrix)[][3], double (*vector_b)[3], int n){

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