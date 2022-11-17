#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

void inverseSub(double (*matrix)[][3], double (*vector_b)[3], int n, bool File, int VCount);
void gauss(double (*matrix)[][3], double (*vector_b)[3], int n, bool Pivot, bool File, int VCount);
void escolhaPivot(double (*matrix)[][3], double (*vector_b)[3], int n, int k);
void function(double (*f)[], double x[]);
void jacobian(double (*jacob)[][3], double x[]);
void newton(double (&x)[], int n);

int main(void){

double matrix[3][3]={{1,-1,-1},{0,2,-1},{0,0,-1}};
double vector_b[3]={0,-2,-7};
bool comPivot = false; 
bool File = false;
int VCount = -10; //Para o grafico de I funcao de V
int n=3;

cout << "Exercício 1" << endl;
cout << "1. Inversa, R3=0" << endl;
inverseSub(&matrix, &vector_b, n, File, VCount);

matrix[2][0]= -2;
cout << "1. Inversa, R3=2" << endl;
inverseSub(&matrix, &vector_b, n, File, VCount);

matrix[2][0]=0;
cout << "1. Gauss, R3=0" << endl;
gauss(&matrix, &vector_b, n, comPivot, File, VCount);


double matrix0[3][3]={{1,-1,-1},{0,2,-1},{-2,0,-1}};
double vector_b0[3]={0,-2,-7};
cout << "1. Gauss R3=2" << endl;
gauss(&matrix0, &vector_b0, n, comPivot, File, VCount);

double matrix1[3][3]={{1,-1,-1},{0,0,-1},{-2,0,-1}};
double vector_b1[3]={0,-2,-7};
cout << "1. Gauss R2=0 e R3=2" << endl;
gauss(&matrix1, &vector_b1, n, comPivot, File, VCount);

double matrix2[3][3]={{1,-1,-1},{0,0,-1},{-2,0,-1}};
double vector_b2[3]={0,-2,-7};
comPivot = true;
cout << "1. Gauss R2=0 e R3=2 e Pivot" << endl;
gauss(&matrix2, &vector_b2, n, comPivot, File, VCount);

fstream I1;                   //Quando o codigo eh executado + que uma vez
I1.open("I1.dat", ios::out);  //Ha necessidade de apagar ficheiros antigos
if(I1){
   I1 << "";
   I1.close();
   I1.open("I2.dat", ios::out);
   I1 << "";
   I1.close();
   I1.open("I3.dat", ios::out);
   I1 << "";
   I1.close();
}

double vector_b3[3]={0,-10,-15};
for(int i=0; i<=20; i++){
    File = true;
    vector_b3[0]={0};
    vector_b3[1]=-10+i;
    vector_b3[2]=-15+i;
    double matrix3[3][3]={{1,-1,-1},{0,2,-1},{-2,0,-1}};
    gauss(&matrix3, &vector_b3, n, comPivot, File, VCount);
    VCount++;
}

double vetorX[2]={1,1};
n=2;
cout << "3. Sistema de Newton" << endl;
newton(vetorX, n);

return 0;
}


void inverseSub(double (*matrix)[][3], double (*vector_b)[3], int n, bool File, int VCount){

int j, i;
double sum=0;
double vector_x[3]={0};
fstream I1FunctionU;
fstream I2FunctionU;
fstream I3FunctionU;


I1FunctionU.open("I1.dat", ios::app);
I2FunctionU.open("I2.dat", ios::app);
I3FunctionU.open("I3.dat", ios::app);

vector_x[n-1] = ((*vector_b)[n-1])/((*matrix)[n-1][n-1]);

for(i=n-1; i>=0; i--){
    sum=0;
    for(j=i+1; j<=n; j++){
        sum += (*matrix)[i][j]*(vector_x[j]);
    }
    vector_x[i]=((*vector_b)[i]-sum)/(*matrix)[i][i];
}


if(File == true){
    I1FunctionU << VCount << "\t" << vector_x[0] << endl;
    I2FunctionU << VCount << "\t" << vector_x[1] << endl;
    I3FunctionU << VCount << "\t" << vector_x[2] << endl;
}
    else if (VCount != 100){
    cout << vector_x[0] << endl;
    cout << vector_x[1] << endl;
    cout << vector_x[2] << endl;
    }

if(VCount == 100){              //Para evitar boilerplate, em Newton, usa-se o vetor f[]  
(*vector_b)[0]=vector_x[0];     //Para armazenar os valores de x
(*vector_b)[1]=vector_x[1];
}

}

void gauss(double (*matrix)[][3], double (*vector_b)[3], int n, bool Pivot, bool File, int VCount){

int i, k, j, l, c;
double m;

for(k=0; k<n-1; k++){
    if(Pivot == 1){
        escolhaPivot(matrix, vector_b, n, k);
    }
    for(i=k+1; i<n; i++){
    m = (*matrix)[i][k] / (*matrix)[k][k];
        for(j=k; j<n; j++){
        (*matrix)[i][j] = (*matrix)[i][j]- m*(*matrix)[k][j];
        } 
        (*vector_b)[i] = (*vector_b)[i]-m*(*vector_b)[k];
    }

}

inverseSub(matrix, vector_b, n, File, VCount);
}

void escolhaPivot(double (*matrix)[][3], double (*vector_b)[3], int n, int k){

int max=k;
int i;
double temp;

for(i=k; i<n; i++){
    if(fabs((*matrix)[i][k]) > fabs((*matrix)[max][k])){
        max=i;
    }
}

if(k != max){
    for(i=k; i<n; i++){
        temp = (*matrix)[k][i];
        (*matrix)[k][i] = (*matrix)[max][i];
        (*matrix)[max][i] = temp;
    }
    temp = (*vector_b)[k];
    (*vector_b)[k] = (*vector_b)[max];
    (*vector_b)[max] = temp;
}
}

void newton(double (&x)[], int n){

double jacob[3][3]={0};
double f[3]={0};
double erro_max=10;
int k, i;


for(k=0; fabs(erro_max)>0.0001 && k<200; k++){
    jacobian(&jacob, x);
    function(&f, x);
    gauss(&jacob, &f, n, true, false, 100);
    x[0]=x[0]-f[0];
    x[1]=x[1]-f[1];
    if(f[0]/x[0] > f[1]/x[1])
        erro_max = f[0]/x[0];
    else
        erro_max = f[1]/x[1];
}

cout << x[0] << "\t" << x[1] << endl;
}

void jacobian(double (*jacob)[][3], double x[]){

(*jacob)[0][0]=2*x[0];
(*jacob)[0][1]=2*x[1];
(*jacob)[1][0]=-2*x[0];
(*jacob)[1][1]=1;

}

void function(double (*f)[], double x[]){

(*f)[0]=x[0]*x[0]+x[1]*x[1]-5;
(*f)[1]=x[1]-x[0]*x[0]+1;

}
