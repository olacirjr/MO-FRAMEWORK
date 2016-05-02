/*#############################################################
  ##                                                         ##
  ##                         Matrix.h                        ##
  ##           -------------------------------------         ##
  ##             Copyright (c) 2004 Chang Wook Ahn           ##
  ##         (Original Copyright (c) 2001 Peter Bosman)      ##
  ##  (Some correction has been made in the orginal codes.)  ##
  ##     Gwangju Institute of Science & Technology (GIST)    ##
  ##                                                         ##
  ##               Main Work: Matrix operations              ##
  ##                                                         ##
  ##  Matrix operations such as inverse and determinant can  ##
  ##    be performed by solving a linear equation.           ##
  ##  The library contains some functions required for       ##
  ##    speedily solving linear equations.                   ##
  ##                                                         ##
  #############################################################*/
 

/************************************************************************
 * Type definitions
 ***********************************************************************/

typedef double **Matrix;
typedef double  *Vector;


/************************************************************************
 * Function prototypes.
 ***********************************************************************/

Matrix Matrix_new( int, int );
Matrix Matrix_clone( Matrix, int, int );
Matrix Matrix_inverse( Matrix, int );
Matrix Matrix_symmetric_inverse( Matrix, int );
Matrix Matrix_inverse_waste_input( Matrix, int );
Matrix Matrix_symmetric_inverse_waste_input( Matrix, int );
void   Matrix_symmetric_inverse_waste_input_place_in( Matrix, int, Matrix );
Vector Vector_new( int );
Vector Vector_clone( Vector, int );
Vector Matrix_Square_Axb( Matrix, int, Vector, short * );
void   Matrix_Square_Axb_waste_input_place_in( Matrix, int, Vector, Vector, short * );
void   Matrix_print( Matrix, int, int );
void   Vector_print( Vector, int );


/*
 * Returns a new Matrix of size n x m (with undefined contents).
 */
Matrix Matrix_new( int n, int m )
{
  int    i;
  Matrix result;

  result = (double **) Malloc( n*( sizeof( double * ) ) );
  for( i = 0; i < n; i++ )
    result[i] = (double *) Malloc( m*( sizeof( double ) ) );

  return( result );
}

/*
 * Returns an indentical clone of a Matrix of size n x m.
 */
Matrix Matrix_clone( Matrix matrix, int n, int m )
{
  int    i, j;
  Matrix result;

  result = Matrix_new( n, m );
  for( i = 0; i < n; i++ )
    for( j = 0; j < m; j++ )
      result[i][j] = matrix[i][j];

  return( result );
}

/*
 * Returns the inverse of a Matrix of size n x n in a new Matrix.
 */
Matrix Matrix_inverse( Matrix matrix, int n )
{
  int    i;
  Matrix clone, inverse;

  clone   = Matrix_clone( matrix, n, n );
  inverse = Matrix_inverse_waste_input( clone, n );
  for( i = 0; i < n; i++ )
    free( clone[i] );
  free( clone );

  return( inverse );
}

/*
 * Returns the inverse of a symmetric Matrix of size n x n in a new Matrix.
 * Because of the symmetry, the inverse can be computed more efficiently.
 */
Matrix Matrix_symmetric_inverse( Matrix matrix, int n )
{
  int    i;
  Matrix clone, inverse;

  clone   = Matrix_clone( matrix, n, n );
  inverse = Matrix_symmetric_inverse_waste_input( clone, n );
  for( i = 0; i < n; i++ )
    free( clone[i] );
  free( clone );

  return( inverse );
}

/*
 * Returns the inverse of a Matrix of size n x n in a new Matrix,
 * altering the given Matrix for computations. The input
 * has therefore most likely changed upon completion.
 */
Matrix Matrix_inverse_waste_input( Matrix matrix, int n )
{
  int    i, j, row, column, row2;
  double factor;
  Matrix inverse;

  /* Initialize: inverse <- identity matrix */
  inverse = Matrix_new( n, n );
  for( i = 0; i < n; i++ )
    for( j = 0; j < n; j++ )
    {
      if( i == j )
        inverse[i][j] = 1.0;
      else
        inverse[i][j] = 0.0;
    }

  /* Downward: perform row operations to make upper triangular */
  for( row = 0; row < n-1; row++ )
  {
    for( row2 = row+1; row2 < n; row2++ )
    {
      factor = matrix[row2][row]/matrix[row][row];
      /* Only do upper trangle for input, lower triangle for inverse */
      for(column = 0; column < n; column++ ) {
        inverse[row2][column] -= factor * inverse[row][column];
        matrix[row2][column] -= factor * matrix[row][column]; 
       }
    }
  }

  /* Upward  : perform row operations to make diagonal with off-diagonal == 0 */
  for( row = n-1; row >= 0; row-- )
  {
    for( row2 = row-1; row2 >= 0; row2-- )
    {
      factor = matrix[row2][row]/matrix[row][row];
      for( column = n-1; column >= 0; column-- )
        inverse[row2][column] -= factor * inverse[row][column];
    }
    /* Divide row by pivot values */
    for( column = n-1; column >= 0; column-- )
      inverse[row][column] /= matrix[row][row];
  }

  return( inverse );
}

/*
 * Returns the inverse of a symmetric Matrix of size n x n in a new Matrix,
 * altering the given Matrix for computations. 
 * The input has therefore most likely changed upon completion.
 */
Matrix Matrix_symmetric_inverse_waste_input( Matrix matrix, int n )
{
  int    i, j, row, column, row2;
  double factor;
  Matrix inverse;

  /* Initialize: inverse <- identity matrix */
  inverse = Matrix_new( n, n );
  for( i = 0; i < n; i++ )
    for( j = 0; j < n; j++ )
    {
      if( i == j )
        inverse[i][j] = 1.0;
      else
        inverse[i][j] = 0.0;
    }

  /* Downward: perform row operations to make upper triangular */
  for( row = 0; row < n-1; row++ )
  {
    for( row2 = row+1; row2 < n; row2++ )
    {
      factor = matrix[row2][row]/matrix[row][row];
      /* Only do upper trangle for input, lower triangle for inverse */
      
      for( column = 0; column < n; column++ ) {
        inverse[row2][column] -= factor * inverse[row][column];
        matrix[row2][column] -= factor * matrix[row][column];
      }
    }
  }

  /* Upward  : perform row operations to make diagonal with off-diagonal == 0 */
  for( row = n-1; row >= 0; row-- )
  {
    for( row2 = row-1; row2 >= 0; row2-- )
    {
      factor = matrix[row2][row]/matrix[row][row];
      /* Only do lower trangle */
      for( column = n-1; column >= 0; column-- )
        inverse[row2][column] -= factor * inverse[row][column];
    }
    /* Divide row by pivot values */
    for( column = n-1; column >= 0; column-- )
      inverse[row][column] /= matrix[row][row];
  }

  /* Copy values to upper triangle */
  for( row = 0; row < n; row++ )
    for( column = row+1; column < n; column++ )
      inverse[row][column] = inverse[column][row];

  return( inverse );
}

/*
 * Computes the inverse of a symmetric Matrix of size n x n
 * and places the result in a given matrix of the same size.
 * The given Matrix is altered for computations. 
 * The input has therefore most likely changed upon completion.
 */
void Matrix_symmetric_inverse_waste_input_place_in( Matrix matrix, int n, Matrix inverse )
{
  int    i, j, row, column, row2;
  double factor, value;

  /* Initialize: inverse <- identity matrix */
  for( i = 0; i < n; i++ )
    for( j = 0; j < n; j++ )
    {
      if( i == j )
        inverse[i][j] = 1.0;
      else
        inverse[i][j] = 0.0;
    }

  /* Downward: perform row operations to make upper triangular */
  for( row = 0; row < n-1; row++ )
  {
    for( row2 = row+1; row2 < n; row2++ )
    {
      factor = (fabs(matrix[row][row]) < UNDERFLOWVALUE) ? 0 : (matrix[row2][row]/matrix[row][row]);
      /* Only do upper trangle for input, lower triangle for inverse */
      for( column = 0; column < n; column++ )
      {
        value = factor * inverse[row][column];
        if( fabs(value) < OVERFLOWVALUE )
          inverse[row2][column] -= value;
      
		  value = factor * matrix[row][column];

        if( fabs(value) < OVERFLOWVALUE )
          matrix[row2][column] -= value;
      }
    }
  }

  /* Upward  : perform row operations to make diagonal with off-diagonal == 0 */
  for( row = n-1; row >= 0; row-- )
  {
    for( row2 = row-1; row2 >= 0; row2-- )
    {
      factor = (fabs(matrix[row][row]) < UNDERFLOWVALUE) ? 0 : (matrix[row2][row]/matrix[row][row]);
      /* Only do lower trangle */
      for( column = n-1; column >= 0; column-- )
      {
        value = factor * inverse[row][column];
        if( fabs(value) < OVERFLOWVALUE )
          inverse[row2][column] -= value;
      }
    }
    /* Divide row by pivot values */
    for( column = n-1; column >= 0; column-- )
      if( fabs(matrix[row][row]) < UNDERFLOWVALUE )
        inverse[row][column] = 0;
      else
        inverse[row][column] /= matrix[row][row];
  }

  /* Copy values to upper triangle */
  for( row = 0; row < n; row++ )
    for( column = row+1; column < n; column++ )
      inverse[row][column] = inverse[column][row];
}

/*
 * Returns a new Vector of size size.
 */
Vector Vector_new( int size )
{
  Vector result;

  result = (double *) Malloc( size*sizeof( double ) );

  return( result );
}

/*
 * Returns an indentical clone of a Vector of size size.
 */
Vector Vector_clone( Vector input, int size )
{
  int    i;
  Vector result;

  result = Vector_new( size );
  for( i = 0; i < size; i++ )
    result[i] = input[i];

  return( result );
}

/*
 * Computes Vector x from Ax = b for a Matrix A of size n x n and a
 * vector b of size n. It returns x, leaving the input unharmed.
 * The short variable singular is set if the Matrix is not invertible.
 */
Vector Matrix_Square_Axb( Matrix matrix, int n, Vector input, short *singular )
{
  int    i;
  Vector result, clone;
  Matrix mclone;

  clone  = Vector_clone( input, n );
  mclone = Matrix_clone( matrix, n, n );
  result = Vector_new( n );
  Matrix_Square_Axb_waste_input_place_in( mclone, n, clone, result, singular );
  free( clone );
  for( i = 0; i < n; i++ )
    free( mclone[i] );
  free( mclone );

  return( result );
}

/*
 * Computes Vector x from Ax = b for a Matrix A of size n x n and a
 * vector b of size n. It places the values of x in a specified
 * cector and also most likely alters the values of both the Matrix
 * and the input vector. The short variable singular is set 
 * if the matrix is not invertible.
 */
void Matrix_Square_Axb_waste_input_place_in( Matrix matrix, int n, Vector input, Vector result, short *singular )
{
  int    row, column, row2;
  double factor;

  (*singular) = 0;
  /* Downward: perform row operations to make upper triangle */
  for( row = 0; row < n-1; row++ )
  {
    for( row2 = row+1; row2 < n; row2++ )
    {
      if( matrix[row][row] == 0 )
        (*singular) = 1;

      factor = matrix[row][row] == 0 ? 0 : matrix[row2][row]/matrix[row][row];
      /* Only do upper trangle for matrix */
      input[row2] -= factor*input[row];
      for( column = 0; column < n; column++ )
        matrix[row2][column] -= factor * matrix[row][column];
    }
  }

  /* Upward: compute values for x by substitution */
  for( row = n-1; row >= 0; row-- )
  {
    for( column = n-1; column > row; column-- )
      input[row] -= result[column]*matrix[row][column];
    result[row] = matrix[row][row] == 0 ? 0 : input[row]/matrix[row][row];
  }
}

/*
 * Prints a Matrix of size n x m on screen row by row.
 */
void Matrix_print( Matrix matrix, int n, int m )
{
  int i, j;

  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < m; j++ )
      printf("%lf ", matrix[i][j]);
    printf("\n");
  }
}

/*
 * Prints a Vector of size size (in a column).
 */
void Vector_print( Vector vector, int size )
{
  int i;

  for( i = 0; i < size; i++ )
    printf("%lf\n", vector[i]);
}
