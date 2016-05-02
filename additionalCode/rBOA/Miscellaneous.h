/*#############################################################
  ##                                                         ##
  ##                       Miscellaneous.h                   ##
  ##           -------------------------------------         ##
  ##             Copyright (c) 2004 Chang Wook Ahn           ##
  ##         (Original Copyright (c) 2001 Peter Bosman)      ##
  ##       (Some addition and correction have been made.)    ##
  ##     Gwangju Institute of Science & Technology (GIST)    ##
  ##                                                         ##
  ##                Main Work: Macro functions               ##
  ##                                                         ##
  ##  This library contains utilities for memory allocation, ##
  ##    time measurement, and some mathematical functions.   ##
  ##                                                         ##
  #############################################################*/


#include <time.h>
#include <sys/timeb.h>


/************************************************************************
 * Global variables.
 ***********************************************************************/

/* pacific standard & daylight savings */
char *tzstr = (char*)"TZ=PST8PDT";





/************************************************************************
 * General macro definitions.
 ***********************************************************************/

#define SAFELOG(VALUE)       (VALUE <= 0 ? 0.0 : log( VALUE ))
#define SAFESQRT(VALUE)      (VALUE <= 0 ? 0.0 : sqrt( VALUE ))
#define LOGBASE(VALUE,BASE)  return( log( VALUE ) / log( BASE ) );
#define PI                   3.1415926536
#define E                    2.7182818284
#define OVERFLOWVALUE        2.247116970028124842000000000000000000000e+100
#define UNDERFLOWVALUE       5.0e-100
#define SMALLVALUE           5.0e-25



/************************************************************************
 * Function prototypes.
 ***********************************************************************/

void   *Malloc( long );
long    timeMe( int );
int     powint( int, int );
double  moduloDouble( double, double );
double  sgn( double );
int     justsgn( double );
int     evensgn( double );
double  doublefac( int );
int     twopow( int );


/*
 * Allocate memory, abort program in case of allocation failure.
 */
void *Malloc( long size )
{
  void *memory;

  memory = (void *) malloc( size );
  if( !memory )
  {
    printf("\n");
    printf("Error: in allocating memory in Malloc( %ld ), aborting program.", size);
    printf("\n");

    exit(1);
  }

  return( memory );
}


/*
 * A function to measure elapsed time.
 * A first call results in -1 as a return value.
 * A second call results in the time since the first call as a return value.
 */
short msufcset1 = 0, msufcset2 = 0;
long  milli_seconds_upon_first_call, milli_seconds_upon_first_call2;
long timeMe( int timer )
{
  short        msufcset;
  long         result, msufc;
  struct timeb t;

  if( timer == 1 )
  {
    msufc    = milli_seconds_upon_first_call;
    msufcset = msufcset1;
  }
  else
  {
    msufc    = milli_seconds_upon_first_call2;
    msufcset = msufcset2;
  }

  if( !msufcset )
  {
    putenv(tzstr);
    tzset();
    ftime(&t);
    msufc = t.time*1000 + t.millitm;

    if( timer == 1 )
    {
      milli_seconds_upon_first_call = msufc;
      msufcset1                     = 1;
    }
    else
    {
      milli_seconds_upon_first_call2 = msufc;
      msufcset2                     = 1;
    }

    return( -1 );
  }

  putenv(tzstr);
  tzset();
  ftime(&t);
  result = t.time*1000 + t.millitm - msufc;

  if( timer == 1 )
  {
    milli_seconds_upon_first_call = -1;
    msufcset1                     = 0;
  }
  else
  {
    milli_seconds_upon_first_call2 = -1;
    msufcset2                     = 0;
  }

  return( result );
}

/*
 * Computes and returns base^power for integers using recursion
 */
int powint( int base, int power )
{
  int value;

  if( power == 0 )
    return( 1 );

  if( power == 1 )
    return( base );

  if( power == 2 )
    return( base*base );

  if( !(power % 2) )
  {
    value = powint( base, power / 2 );

    return( value * value );
  }

  return( base*powint( base, power - 1 ) );
}

/*
 * Returns x modulo m for floats: x - m*floor(x/m)
 */
double moduloDouble( double x, double m )
{
  return( x - m*floor(x/m) );
}

/*
 * Returns +1 if x > 0, 0 if x = 0, -1 if x < 0
 */
double sgn( double x )
{
  if( x > 0 )
    return( 1 );
  if( x < 0 )
    return( -1 );

  return( 0 );
}

/*
 * Returns +1 if x >= 0, -1 if x < 0
 */
int justsgn( double x )
{
  if( x >= 0 )
    return( 1 );

  return( -1 );
}

/*
 * Returns +1 if floor(x) is even, -1 otherwise
 */
int evensgn( double x )
{
  return( int( -sgn( moduloDouble( x, 2.0 ) - 0.5 ) ) );
}

/*
 * Returns factorial as a double value.
 */
double doublefac( int val )
{
  static double facs[13] = {1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600};

  if( val > 12 )
    return( 0.0 );

  return( facs[val] );
}

/*
 * Returns 2^{val} for integers
 */
int twopow( int val )
{
  static int twopows[30] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576,2097152,4194304,8388608,16777216,33554432,67108864,134217728,268435456,536870912};

  return( twopows[val] );
}
