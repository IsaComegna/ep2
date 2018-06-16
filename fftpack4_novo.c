# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include "fftpack4.h"
# include "fftpack4_precision.h"
/*
int main ( );

void test04 ( );
void r8vec_print_part ( int n, double a[], int max_print, char *title );
double *r8vec_uniform_01_new ( int n, int *seed );
void rr8vec_print_part ( int n, double a[], int max_print, char *title );
double *rr8vec_uniform_01_new ( int n, int *seed );
void timestamp ( );
*/


void r8vec_print_part ( int n, double a[], int max_print, char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT_PART prints "part" of an R8VEC.

  Discussion:

    The user specifies MAX_PRINT, the maximum number of lines to print.

    If N, the size of the vector, is no more than MAX_PRINT, then
    the entire vector is printed, one entry per line.

    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    followed by a line of periods suggesting an omission,
    and the last entry.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 February 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double A[N], the vector to be printed.

    Input, int MAX_PRINT, the maximum number of lines
    to print.

    Input, char *TITLE, a title.
*/
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
    }
    fprintf ( stdout, "  ......  ..............\n" );
    i = n - 1;
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
    }
    i = max_print - 1;
    fprintf ( stdout, "  %8d: %14f  ...more entries...\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

double *r8vec_uniform_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01_NEW returns a unit pseudorandom R8VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

void rr8vec_print_part ( int n, double a[], int max_print, char *title )

/******************************************************************************/
/*
  Purpose:

    RR8VEC_PRINT_PART prints "part" of an RR8VEC.

  Discussion:

    The user specifies MAX_PRINT, the maximum number of lines to print.

    If N, the size of the vector, is no more than MAX_PRINT, then
    the entire vector is printed, one entry per line.

    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    followed by a line of periods suggesting an omission,
    and the last entry.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 May 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double A[N*2], the vector to be printed.

    Input, int MAX_PRINT, the maximum number of lines
    to print.

    Input, char *TITLE, a title.
*/
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      fprintf ( stdout, "  %8d: %14f  %14f\n", i, a[i*2+0], a[i*2+1] );
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      fprintf ( stdout, "  %8d: %14f  %14f\n", i, a[i*2+0], a[i*2+1] );
    }
    fprintf ( stdout, "  ......  ..............\n" );
    i = n - 1;
    fprintf ( stdout, "  %8d: %14f  %14f\n", i, a[i*2+0], a[i*2+1] );
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      fprintf ( stdout, "  %8d: %14f  %14f\n", i, a[i*2+0], a[i*2+1] );
    }
    i = max_print - 1;
    fprintf ( stdout, "  %8d: %14f  %14f  ...more entries...\n", i, a[i*2+0], a[i*2+1] );
  }
  return;
}
/******************************************************************************/

double *rr8vec_uniform_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    RR8VEC_UNIFORM_01_NEW returns a unit pseudorandom RR8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 May 2013

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double RR8VEC_UNIFORM_01_NEW[N*2], the vector of pseudorandom values.
*/
{
  double *c;
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;
  double r;
  double pi = 3.1415926;
  double theta;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "RR8VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  c = ( double * ) malloc ( n * 2 * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    theta = 2.0 * pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

    c[i*2+0] = r * cos ( theta );
    c[i+2+1] = r * sin ( theta );
  }
  return c;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}


void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests EZFFTB, EZFFTF, EZFFTI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double azero;
  double *b;
  int i;
  int *ifac;
  int n = 4096;
  int nh;
  int seed;
  double *wsave;
  double *x;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  For real fast Fourier transform,\n" );
  printf ( "  EZFFTI initializes the transform.\n" );
  printf ( "  EZFFTF does a forward transform;\n" );
  printf ( "  EZFFTB does a backward transform.\n" );
  printf ( "\n" );
  printf ( "  The number of data items is N = %d\n", n );
/*
  Set the data values.
*/
  seed = 1973;

  x = r8vec_uniform_01_new ( n, &seed );

  r8vec_print_part ( n, x, 10, "  The original data:" );
/*
  Initialize the WSAVE array.
*/
  wsave = ( double * ) malloc ( ( 3 * n + 15 ) * sizeof ( double ) );
  ifac = ( int * ) malloc ( 8 * sizeof ( int ) );

  ezffti ( &n, wsave, ifac );
/*
  Compute FFT
*/
  nh = n / 2;
  a = ( double * ) malloc ( nh * sizeof ( double ) );
  b = ( double * ) malloc ( nh * sizeof ( double ) );

  ezfftf ( &n, x, &azero, a, b, wsave, ifac );

  printf ( "\n" );
  printf ( "  The A0 coefficient:\n" );
  printf ( "\n" );
  printf ( "  %g\n", azero );

  r8vec_print_part ( n/2, a, 10, "  The A coefficients:" );

  r8vec_print_part ( n/2, b, 10, "  The B coefficients:" );
/*
  Now compute inverse FFT of coefficients.  Should get back the
  original data.  First destroy original data so we're sure
  that the routine had to recreate them!
*/
  printf ( "\n" );
  printf ( "  Retrieve data from FFT coeficients.\n" );

  ezfftb ( &n, x, &azero, a, b, wsave, ifac );

  r8vec_print_part ( n, x, 10, "  The retrieved data:" );

  free ( a );
  free ( b );
  free ( ifac );
  free ( wsave );
  free ( x );

  return;
}
/*************************************************************************/


int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FFTPACK4_PRB.

  Discussion:

    FFTPACK4_PRB tests the FFTPACK4 library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FFTPACK4_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FFTPACK4 library.\n" );


  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FFTPACK4_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
