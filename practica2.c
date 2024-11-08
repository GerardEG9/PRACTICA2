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
        printf("%f ",vect[i]);
    }
}

void PrintRow( float mat[N][N], int row, int from, int numel ){
  for (int i=from; i<=numel; i++){
    printf("%f ",mat[row][i]);
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
for (int i=0; i < N; i++){
    float sum_element_diagonal=0.0;
  for (int j=0; j < N; j++){
    if (i != j){
    sum_element_diagonal+=fabs(M[i][j]);
    }
    }
  if (fabs(M[i][i])<sum_element_diagonal){
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

int Jacobi(float M[N][N], float vect[N], float vectres[N], unsigned iter) {
for(int i=0; i < N; i++){
vectres[i]=0.0;
}
    if (DiagonalDom(M)==1) {
        for (unsigned k = 0; k < iter; k++) {

            for (int i = 0; i < N; i++) {
                float sum = 0.0;

                for (int j = 0; j < N; j++) {
                    if (i != j) {
                        sum -= M[i][j] * vectres[j];
                    }
                }


                vectres[i] = (vect[i] - sum) / M[i][i];
            }


        }
        return 1;
    } else {
        return 0;
    }
}

int main(){
InitData();
float vectres[N];
printf("D.\n\n");
printf("Infininorma de Mat = %f\n",Infininorm(Mat));
printf("Norma u de Mat = %f\n",Onenorm(Mat));
printf("Norma de Forbenius = %f\n",NormFrobenius(Mat));
    if (DiagonalDom(Mat)==1){
    printf("La matriu Mat és diagonal dominant.\n");
    } else {
    printf("La matriu Mat no és diagonal dominant.\n");
}
printf("\n");
    printf("Infininorma de MatDD = %f\n",Infininorm(MatDD));
    printf("Norma u de MatDD = %f\n",Onenorm(MatDD));
    printf("Norma de Forbenius = %f\n",NormFrobenius(MatDD));
    if (DiagonalDom(MatDD)==1){
        printf("La matriu MatDD és diagonal dominant.\n");
    } else {
        printf("La matriu MatDD no és diagonal dominant.\n\n");
    }
    printf("E.\n\n");
printf("Escalar <V1,V2> = %f Escalar <V1,V3> = %f Escalar <V2,V3> = %f\n\n",Scalar(V1,V2),Scalar(V1,V3),Scalar(V2,V3));
printf("F.\n\n");
printf("Magnitud V1,V2 i V3 = %f %f %f\n\n",Magnitude(V1),Magnitude(V2),Magnitude(V3));
printf("G.\n\n");
if (Ortogonal(V1,V2)==1){
printf("V1 i V2 són ortogonals\n");
}
if (Ortogonal(V1,V3)==1){
printf("V1 i V3 són ortogonals\n");
}
if (Ortogonal(V2,V3)==1){
printf("V2 i V3 són ortogonals\n");
}
printf("\n");
printf("H.\n\n");
MultEscalar(V3,vectres,2.0);
printf("Els elements 0 al 9 i 256 al 265 del resultat de multiplicar V3x2.0 són:\n");
PrintVect( vectres, 0, 9 );
printf("\n");
PrintVect( vectres, 256, 265 );
printf("\n\n");
printf("I.\n\n");
Projection( V2, V3,vectres );
printf("Els elements 0 a 9 del resultat de la projecció de V2 sobre V3 són:\n");
PrintVect( vectres, 0, 9 );
printf("\n");
Projection(V1,V2,vectres);
printf("Els elements 0 a 9 del resultat de la projecció de V1 sobre V2 són:\n");
PrintVect( vectres, 0, 9 );
printf("\n\n");
printf("J.\n\n");
Matriu_x_Vector( Mat,V2,vectres );
printf("Els elements 0 a 9 del resultat de la multiplicació de Mat per v2 són:\n");
PrintVect( vectres, 0, 9 );
printf("\n\n");

printf("K.\n\n");
Jacobi(MatDD,V3,vectres,1);
printf("Els elements 0 a 9 de la solució (1 iter) del sistema d'equacions són:\n");
PrintVect( vectres, 0, 9 );
printf("\n");
Jacobi(MatDD,V3,vectres,1000);
printf("Els elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
PrintVect( vectres, 0, 9 );
printf("\n");




}
