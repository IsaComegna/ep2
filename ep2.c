#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "fftpack4.h"


void imprimeVetor (double *vetor, int n){ //imprime vetor de tamanho n;
    int i;
    for (i=0;i<n;i++){
        printf ("%d: %lf \n", i, vetor[i]);
    }
}

void imprimeVetorComplx(double complex *vetor, int n){ //imprime um vetor de complexos de tamanho n;
    int i;
    for (i=0;i<n;i++){
        printf ("%d: %lf + j%lf \n", i, creal(vetor[i]), cimag(vetor[i]) );
    }
}

void fourier(double complex *F, double complex *x, double complex *c, int n){ //recebe os vetores F[x] e x, calcula a transformada e coloca no vetor c

    double complex soma;
    int k, j;
    for (k=0; k<2*n; k++){
        soma = 0;
        for (j=0; j<2*n; j++){
            soma += F[j]*(cpow(M_E , (-I*k*x[j]) ));
        }
        c[k] = soma*((double)1/(double)(2*n));
    }
}

void antiFourier (double complex *F, double complex *x, double complex *c, int n){ //recebe os vetores C e X, calcula a antitransformada e coloca no vetor F

    int j,k;
    double complex soma;

    for (j=0;j<2*n;j++){
        soma = 0;
        for (k=0;k<2*n;k++){
            soma += c[k]*(cpow(M_E, I*k*x[j]));
        }
        F[j] = soma;
    }

}

void fftrec(double complex *c, double complex *f, int n, bool dir){

    double complex *even = (double complex *)malloc(n * sizeof(double complex));
    double complex *odd = (double complex *)malloc(n * sizeof(double complex));
    double complex *fe = (double complex *)malloc(n * sizeof(double complex));
    double complex *fo = (double complex *)malloc(n * sizeof(double complex));
    double complex eij_int;
    double complex eij;
    int j;
    float nn = (float)n;

    if (n==1){
       c[0] = f[0] + f[1];
       c[1] = f[0] - f[1];
   }
   else {
       for (j=0;j<n;j++){
           fe[j]=f[2*j];
           fo[j]=f[2*j+1];
       }
       fftrec(even,fe, n/2, dir);
       fftrec(odd, fo, n/2, dir);
       for (j=0;j<n;j++){
           if (dir){
               eij = cpow(M_E, (-I*j*M_PI)/nn);
               //printf("%le",((I*j*M_PI)/nn));
           }
           else{
               eij = cpow(M_E, (I*j*M_PI)/nn);
                //printf("\n %le",  M_PI);
               //printf("%le",(I*j));
           }
           c[j] = even[j]+eij*odd[j];
           c[j+n] = even[j]-eij*odd[j];
       }
   }
}


int main (){ //teste inicial

    int i, n,k,j = 0;
    double pi = 3.1415926;
    n = 2;
    double complex *x = (double complex *)malloc(2*n * sizeof(double complex)); //x
    double complex *F = (double complex *)malloc(2*n * sizeof(double complex)); //F(x)
    double complex *c = (double complex *)malloc(2*n * sizeof(double complex)); //ck

    x[0] = 0;
    x[1] = pi/2;
    x[2] = pi;
    x[3] = (3*pi)/2;

    F[0] = 5;
    F[1] = -1;
    F[2] = 3;
    F[3] = 1;

    c[0] = 0;
    c[1] = 0;
    c[2] = 0;
    c[3] = 0;

    fourier(F, x, c, n); //calcula a transformada de fourier e coloca os ck's no vetor c

    printf ("Teste a1) Transformada de Fourier de F(x)\n");
    printf ("F(x): \n");
    imprimeVetorComplx(F, 2*n);
    printf ("ck's: \n");
    imprimeVetorComplx(c, 2*n);

    printf ("Teste a2) Antitransformada utilizando os coeficientes ck's: \n");
    double complex *F2 = (double complex *)malloc(2*n * sizeof(double complex));
    antiFourier(F2, x, c, n); //calcula a antitransformada de fourir e coloca os F(X) no vetor F2
    printf ("F(x): \n");
    imprimeVetorComplx(F2, 2*n);


    printf ("Teste a3) Transformada utilizando FFT recursiva \n");
    printf ("Ck's: \n");
    double complex *c2 = (double complex *)malloc(2*n * sizeof(double complex));
    fftrec(c2, F, n, 1);
    for (j=0;j<2*n;j++){   //dividir por 2n pra normalizar, faz caso for transformação direta
        c2[j] = c2[j]/(2*n);
    }
    imprimeVetorComplx(c2, 2*n);


    printf ("Teste a4) Antitransformada utilizando FFT recursiva \n");
    printf ("F(x): \n");
    double complex *F3 = (double complex *)malloc(2*n * sizeof(double complex));
    fftrec(F, c2, n, 0);
    imprimeVetorComplx(F, 2*n);


    n = 8;
    double complex *xb = (double complex *)malloc(2*n * sizeof(double complex)); //x
    double complex *Fb = (double complex *)malloc(2*n * sizeof(double complex)); //F(x)
    double complex *cb = (double complex *)malloc(2*n * sizeof(double complex)); //ck

    F[0] = 6;
    F[1] = 2;
    F[2] = 5;
    F[3] = 2;
    F[4] = 11;
    F[5] = 2;
    F[6] = 8;
    F[7] = 8;

    printf ("Teste b) Comparação com rotinas FFTPACK4 \n");
    printf ("F(x): \n");


    //test04();
    double *a;//array de numeros reais contendo as partes reais dos coef. complexos de fourier,
                //se n é impar, b tem tamanhao n/2, senao tamanho (n-1)/2

    double azero; //constante de fourier A0
    double *b; //array de numeros reais contendo as partes imaginarias dos coef. complexos de fourier,
                //se n é impar, b tem tamanhao n/2, senao tamanho (n-1)/2
    int *ifac;
    n = 8; //tamanho da sequencia
    int nh; //tamanho de a e b, dependendo se n é par ou nao
    //int seed;
    double *wsave; //work array inicializado pelo ezffti
    double *xx = (double *)malloc(2*n * sizeof(double));

    xx[0] = 6;
    xx[1] = 2;
    xx[2] = 5;
    xx[3] = 2;
    xx[4] = 11;
    xx[5] = 2;
    xx[6] = 8;
    xx[7] = 8;

    wsave = ( double * ) malloc ( ( 3 * n + 15 ) * sizeof ( double ) );
    ifac = ( int * ) malloc ( 8 * sizeof ( int ) );
    ezffti ( &n, wsave, ifac );
    imprimeVetorComplx(xx, 8);
    //r8vec_print_part ( n, x, 10, "  The original data:" );
    nh = n / 2;
    a = ( double * ) malloc ( nh * sizeof ( double ) );
    b = ( double * ) malloc ( nh * sizeof ( double ) );
    ezfftf ( &n, xx, &azero, a, b, wsave, ifac );
    printf ( "\n" );
    printf ( "  The A0 coefficient:\n" );
    printf ( "\n" );
    printf ( "  %g\n", azero );
    r8vec_print_part ( n/2, a, n, "  The A coefficients:" );
    r8vec_print_part ( n/2, b, n, "  The B coefficients:" );
    printf ( "\n" );
    printf ( "  Retrieve data from FFT coeficients.\n" );
    ezfftb ( &n, x, &azero, a, b, wsave, ifac );
    r8vec_print_part ( n, x, n, " The retrieved data:" );

    timestamp();



    return 0;
}
