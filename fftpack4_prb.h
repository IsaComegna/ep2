#ifndef FFTPACK4_PRB
#define FFTPACK4_PRB


# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "fftpack4.h"
# include "fftpack4_precision.h"

#ifdef __cplusplus
extern "C" {
#endif
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void r8vec_print_part ( int n, double a[], int max_print, char *title );
double *r8vec_uniform_01_new ( int n, int *seed );
void rr8vec_print_part ( int n, double a[], int max_print, char *title );
double *rr8vec_uniform_01_new ( int n, int *seed );
void timestamp ( );

#ifdef __cplusplus
}
#endif

#endif
