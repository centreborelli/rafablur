/*----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

/*----------------------------------------------------------------------------*/
/*
   Fatal error, print a message to standard-error output and exit.
 */
void error(char * msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/** Memory allocation, print an error and exit if fail.
 */
void * xmalloc(size_t size)
{
  void * p;
  if( size == 0 ) error("xmalloc: zero size");
  p = malloc(size);
  if( p == NULL ) error("xmalloc: out of memory");
  return p;
}

/*----------------------------------------------------------------------------*/
/** Open file, print an error and exit if fail.
 */
FILE * xfopen(const char * path, const char * mode)
{
  FILE * f = fopen(path,mode);
  if( f == NULL ) error("xfopen: unable to open file");
  return f;
}

/*----------------------------------------------------------------------------*/
/** Close file, print an error and exit if fail.
 */
int xfclose(FILE * f)
{
  if( fclose(f) == EOF ) error("xfclose: unable to close file");
  return 0;
}

/*----------------------------------------------------------------------------*/
double * read_asc(char * name, int * X, int * Y, int * Z, int * C)
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
void write_asc(double * image, int X, int Y, int Z, int C, char * name)
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
/*----------------------------- Gaussian filter ------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Compute a Gaussian kernel of length 'n',
    standard deviation 'sigma', and centered at value 'mean'.

    For example, if mean=0.5, the Gaussian will be centered
    in the middle point between values 'kernel[0]' and 'kernel[1]'.
 */
static void gaussian_kernel(double * kernel, int n, double sigma, double mean)
{
  double sum = 0.0;
  double val;
  int i;

  /* check parameters */
  if( kernel == NULL ) error("gaussian_kernel: null 'kernel' pointer.");
  if( sigma <= 0.0 ) error("gaussian_kernel: 'sigma' must be positive.");

  /* compute Gaussian kernel */
  for(i=0; i<n; i++)
    {
      val = ( (double) i - mean ) / sigma;
      kernel[i] = exp( -0.5 * val * val );
      sum += kernel[i];
    }

  /* normalization */
  if( sum >= 0.0 )
    for(i=0; i<n; i++)
      kernel[i] /= sum;
}

/*----------------------------------------------------------------------------*/
/** Scale the input image 'in' by a factor 'scale' by Gaussian sub-sampling.

    For example, scale=0.8 will give a result at 80% of the original size.

    The image is convolved with a Gaussian kernel
    @f[
        G(x,y) = \frac{1}{2\pi\sigma^2} e^{-\frac{x^2+y^2}{2\sigma^2}}
    @f]
    before the sub-sampling to prevent aliasing.

    The standard deviation sigma given by:
    -  sigma = sigma_scale / scale,   if scale <  1.0
    -  sigma = sigma_scale,           if scale >= 1.0

    To be able to sub-sample at non-integer steps, some interpolation
    is needed. In this implementation, the interpolation is done by
    the Gaussian kernel, so both operations (filtering and sampling)
    are done at the same time. The Gaussian kernel is computed
    centered on the coordinates of the required sample. In this way,
    when applied, it gives directly the result of convolving the image
    with the kernel and interpolated to that particular position.

    A fast algorithm is done using the separability of the Gaussian
    kernel. Applying the 2D Gaussian kernel is equivalent to applying
    first a horizontal 1D Gaussian kernel and then a vertical 1D
    Gaussian kernel (or the other way round). The reason is that
    @f[
        G(x,y) = G(x) * G(y)
    @f]
    where
    @f[
        G(x) = \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{x^2}{2\sigma^2}}.
    @f]
    The algorithm first applies a combined Gaussian kernel and sampling
    in the x axis, and then the combined Gaussian kernel and sampling
    in the y axis.
 */
double * gaussian_sampler( double * in, int X, int Y, double scale,
                           double sigma_scale, int * XX, int * YY )
{
  double * aux;
  double * out;
  double * kernel;
  int N,M,h,n,x,y,i,xc,yc,j,double_x_size,double_y_size;
  double sigma,xx,yy,sum,prec;

  /* check parameters */
  if( in == NULL || X == 0 || Y == 0 )
    error("gaussian_sampler: invalid image.");
  if( scale <= 0.0 ) error("gaussian_sampler: 'scale' must be positive.");
  if( sigma_scale <= 0.0 )
    error("gaussian_sampler: 'sigma_scale' must be positive.");

  /* compute new image size and get memory for images */
  if( X * scale > (double) UINT_MAX || Y * scale > (double) UINT_MAX )
    error("gaussian_sampler: the output image size exceeds the handled size.");
  N = (int) ceil( X * scale );
  M = (int) ceil( Y * scale );
  aux = (double *) xmalloc( N * Y * sizeof(double) );
  out = (double *) xmalloc( N * M * sizeof(double) );

  /* sigma, kernel size and memory for the kernel */
  sigma = scale < 1.0 ? sigma_scale / scale : sigma_scale;
  /*
     The size of the kernel is selected to guarantee that the
     the first discarded term is at least 10^prec times smaller
     than the central value. For that, h should be larger than x, with
       e^(-x^2/2sigma^2) = 1/10^prec.
     Then,
       x = sigma * sqrt( 2 * prec * ln(10) ).
   */
  prec = 3.0;
  h = (int) ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
  n = 1+2*h; /* kernel size */
  kernel = (double *) xmalloc( n * sizeof(double) );

  /* auxiliary double image size variables */
  double_x_size = 2 * X;
  double_y_size = 2 * Y;

  /* First subsampling: x axis */
  for(x=0; x<N; x++)
    {
      /*
         x   is the coordinate in the new image.
         xx  is the corresponding x-value in the original size image.
         xc  is the integer value, the pixel coordinate of xx.
       */
      xx = (double) x / scale;
      /* coordinate (0.0,0.0) is in the center of pixel (0,0),
         so the pixel with xc=0 get the values of xx from -0.5 to 0.5 */
      xc = (int) floor( xx + 0.5 );
      gaussian_kernel( kernel, n, sigma, (double) h + xx - (double) xc );
      /* the kernel must be computed for each x because the fine
         offset xx-xc is different in each case */

      for(y=0; y<Y; y++)
        {
          sum = 0.0;
          for(i=0; i<n; i++)
            {
              j = xc - h + i;

              /* symmetry boundary condition */
              while( j < 0 ) j += double_x_size;
              while( j >= double_x_size ) j -= double_x_size;
              if( j >= X ) j = double_x_size-1-j;

              sum += in[j+y*X] * kernel[i];
            }
          aux[x+y*N] = sum;
        }
    }

  /* Second subsampling: y axis */
  for(y=0; y<M; y++)
    {
      /*
         y   is the coordinate in the new image.
         yy  is the corresponding x-value in the original size image.
         yc  is the integer value, the pixel coordinate of xx.
       */
      yy = (double) y / scale;
      /* coordinate (0.0,0.0) is in the center of pixel (0,0),
         so the pixel with yc=0 get the values of yy from -0.5 to 0.5 */
      yc = (int) floor( yy + 0.5 );
      gaussian_kernel( kernel, n, sigma, (double) h + yy - (double) yc );
      /* the kernel must be computed for each y because the fine
         offset yy-yc is different in each case */

      for(x=0; x<N; x++)
        {
          sum = 0.0;
          for(i=0; i<n; i++)
            {
              j = yc - h + i;

              /* symmetry boundary condition */
              while( j < 0 ) j += double_y_size;
              while( j >= double_y_size ) j -= double_y_size;
              if( j >= Y ) j = double_y_size-1-j;

              sum += aux[x+j*N] * kernel[i];
            }
          out[x+y*N] = sum;
        }
    }

  /* free memory */
  free( (void *) kernel );
  free( (void *) aux );

  *XX = N;
  *YY = M;

  return out;
}

/*----------------------------------------------------------------------------*/
/*                                    Main                                    */
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv)
{
  double * image;
  double * output;
  int X,Y,Z,C,XX,YY;
  double scale;
  double sigma_scale = 0.6;

  /* read input */
  if( argc != 4 && argc != 5 )
    error("use: image_gaussian_sampler <scale> <input.asc> <output.asc>"
          " [sigma_scale]");
  scale = atof(argv[1]);
  image = read_asc(argv[2],&X,&Y,&Z,&C);
  if( Z!=1 || C!=1 ) error("Z!=1 or C!=1");
  if( argc == 5 ) sigma_scale = atof(argv[4]);

  /* sub-sample image */
  output = gaussian_sampler( image, X, Y, scale, sigma_scale, &XX, &YY );

  /* write output */
  write_asc(output,XX,YY,1,1,argv[3]);

  /* free memory */
  free( (void *) image );
  free( (void *) output );

  return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
