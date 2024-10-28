#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 512
float Mat[N][N];
float MatDD[N][N];
float V1[N];
float V2[N];
float V3[N];
float V4[N];

void InitData(){
    int i,j;
    srand(334411);
    for( i = 0; i < N; i++ )
        for( j = 0; j < N; j++ ){
            Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
            if ( (abs(i - j) <= 3) && (i != j))
                MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
            else if ( i == j )
                MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
            else MatDD[i][j] = 0.0;
        }
    for( i = 0; i < N; i++ ){
        V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
        V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
        V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
    }
}

void PrintVect( float vect[N], int from, int numel )
{
    for (int i=from; i<=numel; i++){
        printf("%f\n",vect[i]);
    }
}

void PrintRow( float mat[N][N], int row, int from, int numel ){
  for (int i=from; i<=numel; i++){
    printf("%f\n",mat[row][i]);
  }
 }
void MultEscalar( float vect[N], float vectres[N], float alfa ){
for (int i=0; i < N; i++){
vectres[i] = vect[i]*alfa;
}
}
float Scalar( float vect1[N], float vect2[N] ){
float sum = 0.0;
for (int i=0; i < N; i++){
sum += vect1[i]*vect2[i];
}
return sum;
}
float Magnitude( float vect[N] ){
 float sum = 0.0;
 for (int i=0; i < N; i++){
 sum+=pow(vect[i],2); // Elevem vect[i] al quadrat
 }
 return sqrt(sum);
 }
int Ortogonal( float vect1[N], float vect2[N] ){
float num=Scalar(vect1,vect2); // Utilitzem la funcio scalar feta previament per calcular el producte escalar
if (num==0){
return 1;
} else {
return 0;
}
}

void Projection( float vect1[N], float vect2[N], float vectres[N] ){
float num_escalar=Scalar(vect1,vect2);
float magnitude=Magnitude(vect2);
float operacio=num_escalar/magnitude;
if (magnitude!=0){
MultEscalar (vect2,vectres,operacio);
}
}
float Infininorm( float M[N][N] ){
float sum_fila_max=0.0;
    for (int i=0; i < N; i++){
        float sum_fila= 0.0;
        for (int j=0; j < N; j++){
            sum_fila+=fabs(M[i][j]);
        }
        if (sum_fila>sum_fila_max){
        sum_fila_max=sum_fila;
        }
    }
    return sum_fila_max;
}

float Onenorm( float M[N][N] ){
float sum_columna_max=0.0;
    for (int j=0; j < N; j++){
        float sum_columna= 0.0;
        for (int i=0; i < N; i++){
            sum_columna+=fabs(M[i][j]);
        }
        if (sum_columna>sum_columna_max){
        sum_columna_max=sum_columna;
        }
    }
    return sum_columna_max;
}

float NormFrobenius( float M[N][N] ){
float sum=0.0;
  for (int i=0; i < N; i++){
    for (int j=0; j < N; j++){
      sum+=M[i][j]*M[i][j];
      }
    }
  return sqrt(sum);
}
int DiagonalDom( float M[N][N] ){
float sum_element_diagonal=0.0;
for (int i=0; i < N; i++){
  for (int j=0; j < N; j++){
    if (i =! j){
    sum_element_diagonal+=fabs(M[i][j]);
    }
    }

  if (M[i][i]<sum_element_diagonal){
    return 0;
    }
  }
  return 1;
}

void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] ){
  for (int i=0; i < N; i++){
    for (int j=0; j < N; j++){
        vectres[i] += M[i][j]*vect[j];

        }
    }
  }
int main(){
    InitData();
    PrintRow(MatDD,0,0,9);
    PrintRow(MatDD,100,95,104);
}
