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


void pegarDados(double complex *c1, double complex *c2, double complex *x, int tam, char *nome, int *ncanais, int *sampleRate){ //pega os dados do arquivo de dados e preenche x, c1, c2 e armazena ncanais e sampleRate

    FILE *arq = fopen (nome, "r");
    char c[10];
    int i=0;

    fscanf(arq, "%s", &c); //ponto e virgula
    fscanf(arq, "%s", &c); //Sample
    fscanf(arq, "%s", &c); //Rate
    fscanf(arq, "%d", sampleRate);
    fscanf(arq, "%s", &c); //ponto e virgula
    fscanf(arq, "%s", &c); //Channels
    fscanf(arq, "%d", ncanais);

    for (i=0;i<tam;i++){
        fscanf(arq, "%lf", &x[i]);
        fscanf(arq, "%lf", &c1[i]);
        if (*ncanais == 2){
            fscanf(arq, "%lf", &c2[i]);
        }
    }

    fclose(arq);

}

int calcQtdDados(char *nome){    //calcula a quantidade de dados no arquivo com nome "nome"

    FILE *arq = fopen (nome, "r");
    int tam = 0;
    char c[10];
    double aux =0;


    fscanf(arq, "%s", &c); //ponto e virgula
    fscanf(arq, "%s", &c); //Sample
    fscanf(arq, "%s", &c); //Rate
    fscanf(arq, "%lf", &aux);
    fscanf(arq, "%s", &c); //ponto e virgula
    fscanf(arq, "%s", &c); //Channels
    fscanf(arq, "%lf", &aux);

    while (!feof(arq)){
        fscanf(arq, "%lf", &aux);
        fscanf(arq, "%lf", &aux);
        fscanf(arq, "%lf", &aux);
        tam++;
    }
    fclose(arq);
    tam--;

    return tam;
}

void imprimirDados(double complex *c1, double complex *c2, double complex *x, int tam, char *nome, int ncanais, int sampleRate){

    FILE *arq = fopen (nome, "w");

    fprintf(arq, "%s ", ";");
    fprintf(arq, "%s", "Sample Rate ");
    fprintf(arq, "%d", sampleRate);
    fprintf (arq, "\n%s ", ";");
    fprintf (arq, "%s", "Channels ");
    fprintf (arq, "%d", ncanais);
    fprintf (arq, "\n");
    int i;
    for (i=0;i<tam;i++){
        fprintf(arq, "%lf ", creal(x[i]));
        fprintf(arq, "%lf ", creal(c1[i]));
        if (ncanais == 1){
            fprintf(arq, "\n");
        }
        else {
            fprintf(arq, "%lf \n", creal(c2[i]));
        }
    }

}

void fPassaBaixas(double complex *c, int k, int tam){
    int i;

    for (i=k+1;i<tam;i++){
        c[i] = 0;
    }
}

void fPassaAltas(double complex *c, int k, int tam){
    int i;
    for (i=k-1;i>0;i--){
        c[i] = 0;
    }
}

void fPassaBandas(double complex *c, int k1, int k2, int tam){

    int i;
    for (i=k1-1;i>0;i--){
        c[i] = 0;
    }
    for (i=k2+1;i<tam;i++){
        c[i] = 0;
    }
}

void compressao(double complex *c, double rate, int tam){

    int i;
    double amp = 0;

    for (i=0;i<tam;i++){
        amp = 2*cabs(c[i]);
        if (amp < rate){
            c[i] = 0;
        }
    }

}


void testeA (){
    int i, n,k,j = 0;
    n = 2;
    double complex *x = (double complex *)malloc(2*n * sizeof(double complex)); //x
    double complex *F = (double complex *)malloc(2*n * sizeof(double complex)); //F(x)
    double complex *c = (double complex *)malloc(2*n * sizeof(double complex)); //ck

    x[0] = 0;
    x[1] = M_PI/2; //pi na biblioteca math.h
    x[2] = M_PI;
    x[3] = (3*M_PI)/2;

    F[0] = 5;
    F[1] = -1;
    F[2] = 3;
    F[3] = 1;

    c[0] = 0;
    c[1] = 0;
    c[2] = 0;
    c[3] = 0;

    fourier(F, x, c, n); //calcula a transformada de fourier e coloca os ck's no vetor c

    printf ("Teste a1) Transformada de Fourier de F(x)\n\n");
    printf ("F(x): \n");
    imprimeVetorComplx(F, 2*n);
    printf ("ck's: \n");
    imprimeVetorComplx(c, 2*n);

    printf ("\nTeste a2) Antitransformada utilizando os coeficientes ck's: \n\n");
    double complex *F2 = (double complex *)malloc(2*n * sizeof(double complex));
    antiFourier(F2, x, c, n); //calcula a antitransformada de fourir e coloca os F(X) no vetor F2
    printf ("F(x): \n");
    imprimeVetorComplx(F2, 2*n);


    printf ("\nTeste a3) Transformada utilizando FFT recursiva \n\n");
    printf ("Ck's: \n");
    double complex *c2 = (double complex *)malloc(2*n * sizeof(double complex));
    fftrec(c2, F, n, 1);
    for (j=0;j<2*n;j++){   //dividir por 2n pra normalizar, faz caso for transformação direta
        c2[j] = c2[j]/(2*n);
    }
    imprimeVetorComplx(c2, 2*n);


    printf ("\nTeste a4) Antitransformada utilizando FFT recursiva \n\n");
    printf ("F(x): \n");
    double complex *F3 = (double complex *)malloc(2*n * sizeof(double complex));
    fftrec(F3, c2, n, 0);
    imprimeVetorComplx(F3, 2*n);

    free(F);
    free(F2);
    free(F3);
    free(c);
    free(c2);
}

void calcularX (double complex *x, int tam){ //calcula o vetor x a partir do tamanho do vetor

    double complex delta_x = (2*M_PI)/tam;
    int i;

    x[0] = 0;
    for (i=1;i<tam;i++){
        x[i] = x[i-1] + delta_x;
    }

}

void fttpack4Cks(double complex *c, double A0, double *A, double *B, int n){
    /* ak = Ak + i*Bk */

    double complex *ak = (double complex *)malloc(n * sizeof(double complex));

    c[0] = A0;
/*
    for (int i=1; i<n;i++){

    }*/

/*    printf("to aqui moms \n");
    imprimeVetorComplx(c, n);
    printf("\n\n");*/

}

void testeB() {

    int i, n,k,j = 0;
    n = 4;
    double complex *x = (double complex *)malloc(2*n * sizeof(double complex)); //x
    double complex *F = (double complex *)malloc(2*n * sizeof(double complex)); //F(x)
    double complex *c = (double complex *)malloc(2*n * sizeof(double complex)); //ck

    calcularX(x, 8);

    F[0] = 6;
    F[1] = 2;
    F[2] = 5;
    F[3] = 2;
    F[4] = 11;
    F[5] = 2;
    F[6] = 8;
    F[7] = 8;

    fourier(F, x, c, n); //calcula a transformada de fourier e coloca os ck's no vetor c

    printf ("Teste b1) Transformada de Fourier de F(x)\n");
    printf ("F(x): \n");
    imprimeVetorComplx(F, 2*n);
    printf ("ck's: \n");
    imprimeVetorComplx(c, 2*n);

    printf ("Teste b2) Antitransformada utilizando os coeficientes ck's: \n");
    double complex *F2 = (double complex *)malloc(2*n * sizeof(double complex));
    antiFourier(F2, x, c, n); //calcula a antitransformada de fourir e coloca os F(X) no vetor F2
    printf ("F(x): \n");
    imprimeVetorComplx(F2, 2*n);

    printf ("Teste b3) Transformada utilizando FFT recursiva \n");
    printf ("Ck's: \n");
    double complex *c2 = (double complex *)malloc(2*n * sizeof(double complex));
    fftrec(c2, F, n, 1);
    for (j=0;j<2*n;j++){   //dividir por 2n pra normalizar, faz caso for transformação direta
        c2[j] = c2[j]/(2*n);
    }
    imprimeVetorComplx(c2, 2*n);

    printf ("Teste b4) Antitransformada utilizando FFT recursiva \n");
    printf ("F(x): \n");
    double complex *F3 = (double complex *)malloc(2*n * sizeof(double complex));
    fftrec(F3, c2, n, 0);
    imprimeVetorComplx(F3, 2*n);

    free(F);
    free(F2);
    free(F3);
    free(c);
    free(c2);


    //double *xb = (double *)malloc(2*n * sizeof(double)); //x
    double *Fb = (double *)malloc(2*n * sizeof(double)); //F(x)
    //double complex *cs = (double complex *)malloc(2*n * sizeof(double complex)); //ck


    Fb[0] = 6;
    Fb[1] = 2;
    Fb[2] = 5;
    Fb[3] = 2;
    Fb[4] = 11;
    Fb[5] = 2;
    Fb[6] = 8;
    Fb[7] = 8;

    printf ("Teste b5) Transformada de Fourier utilizando a rotina do FFTPACK4 \n");
    double *a;//array de numeros reais contendo as partes reais dos coef. complexos de fourier,
                //se n é impar, b tem tamanhao n/2, senao tamanho (n-1)/2
    double azero; //constante de fourier A0
    double *b; //array de numeros reais contendo as partes imaginarias dos coef. complexos de fourier,
                //se n é impar, b tem tamanhao n/2, senao tamanho (n-1)/2
    int *ifac;
    double *wsave; //work array inicializado pelo ezffti
    wsave = ( double * ) malloc ( ( 3 * n + 15 ) * sizeof ( double ) );
    ifac = ( int * ) malloc ( 8 * sizeof ( int ) );
    a = ( double * ) malloc ( (n/2) * sizeof ( double ) );
    b = ( double * ) malloc ( (n/2) * sizeof ( double ) );
    n=8;

    r8vec_print_part ( n, Fb, n, "F(x)" );

    ezffti ( &n, wsave, ifac ); //inicialização para funções a serem utilizadas

    ezfftf ( &n, Fb, &azero, a, b, wsave, ifac ); //calculo da transformada

    printf ( "\n" );
    printf ( "  Coeficiente A0: %g \n", azero);
    r8vec_print_part ( n/2, a, 10, "  Coeficientes A:" );
    r8vec_print_part ( n/2, b, 10, "  Coeficientes B:" );
    printf ( "\n" );

    printf ("Teste b6) Antitransformada de Fourier utilizando a rotina do FFTPACK4 \n");
    ezfftb ( &n, Fb, &azero, a, b, wsave, ifac );
    r8vec_print_part ( n, Fb, n, " The retrieved data:" );

 //  fttpack4Cks(c, azero,a, b, n);
//    imprimeVetorComplx(c, n*2);



}

void testeC() {

    int i, n,k,j = 0;
    n = 1024/2;

    double complex *x = (double complex *)malloc(2*n * sizeof(double complex));
    double complex *F = (double complex *)malloc(2*n * sizeof(double complex));
    double complex *c = (double complex *)malloc(2*n * sizeof(double complex));

    calcularX(x, 2*n);

    for (i=0;i<2*n;i++){
        F[i] = 10*sin(x[i]) + 7*cos(30*x[i]) + 11*sin(352*x[i]) - 8*cos(711*x[i]);
    }


    fourier(F, x, c, n); //calcula a transformada de fourier e coloca os ck's no vetor c

    printf ("Teste c1) Transformada de Fourier de F(x)\n");
    printf ("F(x): \n");
    imprimeVetorComplx(F, 2*n);
    printf ("ck's: \n");
    imprimeVetorComplx(c, 2*n);

    printf ("Teste c2) Antitransformada utilizando os coeficientes ck's: \n");
    double complex *F2 = (double complex *)malloc(2*n * sizeof(double complex));
    antiFourier(F2, x, c, n); //calcula a antitransformada de fourir e coloca os F(X) no vetor F2
    printf ("F(x): \n");
    imprimeVetorComplx(F2, 2*n);

    free(F);
    free(c);
    free(F2);

}


int main (){ //teste inicial

    printf ("---- TESTE A ----  \n \n");
    testeA();
    system("pause");

    printf ("---- TESTE B ----  \n \n");
    testeB();
    system("pause");

    printf ("---- TESTE C ----  \n \n");
    testeC();
    system("pause");


    return 0;
}
