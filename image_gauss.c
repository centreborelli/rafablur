/*----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*----------------------------------------------------------------------------*/
/*
   Fatal error, print a message to standard-error output and exit.
 */
static void error(char * msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/** Memory allocation, print an error and exit if fail.
 */
static void * xmalloc(size_t size)
{
  void * p;
  if( size == 0 ) error("xmalloc: zero size");
  p = malloc(size);
  if( p == NULL ) error("xmalloc: out of memory");
  return p;
}

/*----------------------------------------------------------------------------*/
/** Compute a Gaussian kernel of length 'n', standard deviation 'sigma',
    and centered at value 'mean'.

    For example, if mean=0.5, the Gaussian will be centered
    in the middle point between values 'kernel[0]'
    and 'kernel[1]'.
 */
static void gaussian_kernel(double * kernel, int n, double sigma, double mean)
{
  double sum = 0.0;
  double val;
  int i;

  /* check parameters */
  if( kernel == NULL || n < 1 )
    error("gaussian_kernel: invalid 'kernel'.");
  if( sigma <= 0.0 ) error("gaussian_kernel: 'sigma' must be positive.");

  /* compute Gaussian kernel */
  for(i=0;i<n;i++)
    {
      val = ( (double) i - mean ) / sigma;
      kernel[i] = exp( -0.5 * val * val );
      sum += kernel[i];
    }

  /* normalization */
  if( sum >= 0.0 ) for(i=0;i<n;i++) kernel[i] /= sum;
}

/*----------------------------------------------------------------------------*/
void gaussian_filter(double * image, int X, int Y, double sigma)
{
  int x,y,offset,i,j,nx2,ny2,n;
  double * kernel;
  double * tmp;
  double val,prec;

  if( sigma <= 0.0 ) error("gaussian_filter: 'sigma' must be positive.");
  if( image == NULL || X < 1 || Y < 1 )
    error("gaussian_filter: invalid image.");

  /* create temporary image */
  tmp = (double *) xmalloc( X * Y * sizeof(double) );

  /* compute gaussian kernel */
  /*
     The size of the kernel is selected to guarantee that the
     the first discarted term is at least 10^prec times smaller
     than the central value. For that, h should be larger than x, with
       e^(-x^2/2sigma^2) = 1/10^prec.
     Then,
       x = sigma * sqrt( 2 * prec * ln(10) ).
   */
  prec = 3.0;
  offset = ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
  n = 1 + 2 * offset; /* kernel size */
  kernel = (double *) xmalloc( n * sizeof(double) );
  gaussian_kernel(kernel, n, sigma, (double) offset);

  /* auxiliary double image size variables */
  nx2 = 2*X;
  ny2 = 2*Y;

  /* x axis convolution */
  for(x=0; x<X; x++)
    for(y=0; y<Y; y++)
      {
        val = 0.0;
        for(i=0; i<n; i++)
          {
            j = x - offset + i;

            /* symmetry boundary condition */
            while(j<0) j += nx2;
            while(j>=nx2) j -= nx2;
            if( j >= X ) j = nx2-1-j;

            val += image[j+y*X] * kernel[i];
          }
        tmp[x+y*X] = val;
      }

  /* y axis convolution */
  for(x=0; x<X; x++)
    for(y=0; y<Y; y++)
      {
        val = 0.0;
        for(i=0; i<n; i++)
          {
            j = y - offset + i;

            /* symmetry boundary condition */
            while(j<0) j += ny2;
            while(j>=ny2) j -= ny2;
            if( j >= Y ) j = ny2-1-j;

            val += tmp[x+j*X] * kernel[i];
          }
        image[x+y*X] = val;
      }

  /* free memory */
  free( (void *) kernel );
  free( (void *) tmp );
}

#ifndef IMAGE_GAUSS_OMIT_MAIN

/*----------------------------------------------------------------------------*/
/** Open file, print an error and exit if fail.
 */
static FILE * xfopen(const char * path, const char * mode)
{
  FILE * f = fopen(path,mode);
  if( f == NULL ) error("xfopen: unable to open file");
  return f;
}

/*----------------------------------------------------------------------------*/
/** Close file, print an error and exit if fail.
 */
static int xfclose(FILE * f)
{
  if( fclose(f) == EOF ) error("xfclose: unable to close file");
  return 0;
}

/*----------------------------------------------------------------------------*/
static double * read_asc(char * name, int * X, int * Y, int * Z, int * C)
{
  FILE * f;
  int i,n;
  double val;
  double * image;

  /* open file */
  f = xfopen(name,"r");

  /* read header */
  n = fscanf(f,"%u%*c%u%*c%u%*c%u",X,Y,Z,C);
  if( n!=4 || *X<=0 || *Y<=0 || *Z<=0 || *C<=0 )
    error("read_asc: invalid asc file A");

  /* get memory */
  image = (double *) xmalloc( *X * *Y * *Z * *C * sizeof(double) );

  /* read data */
  for(i=0; i<(*X * *Y * *Z * *C); i++)
    {
      n = fscanf(f,"%lf%*[^0-9.eE+-]",&val);
      if( n!=1 ) error("read_asc: invalid asc file");
      image[i] = val;
    }

  /* close file */
  xfclose(f);

  return image;
}

/*----------------------------------------------------------------------------*/
static void write_asc(double * image, int X, int Y, int Z, int C, char * name)
{
  FILE * f;
  int i;

  /* check input */
  if( image == NULL || X < 1 || Y < 1 || Z < 1 || X < 1 )
    error("write_asc: invalid image");

  f = xfopen(name,"w");                                  /* open file */
  fprintf(f,"%u %u %u %u\n",X,Y,Z,C);                    /* write header */
  for(i=0; i<X*Y*Z*C; i++) fprintf(f,"%.16g ",image[i]); /* write data */
  xfclose(f);                                            /* close file */
}


/*----------------------------------------------------------------------------*/
/*                                    Main                                    */
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv)
{
  double * image;
  int X,Y,Z,C;
  double sigma;

  /* get input */
  if( argc < 4 ) error("use: image_gauss <sigma> <input.asc> <output.asc>");
  sigma = atof(argv[1]);
  image = read_asc(argv[2],&X,&Y,&Z,&C);
  if( Z!=1 || C!=1 ) error("Z!=1 or C!=1");

  /* process */
  gaussian_filter(image,X,Y,sigma);

  /* write output */
  write_asc(image,X,Y,Z,C,argv[3]);

  return 0;
}
/*----------------------------------------------------------------------------*/

#endif//IMAGE_GAUSS_OMIT_MAIN
