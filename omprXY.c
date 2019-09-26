/*omprXY.c is for R, generates x and y data for n random points, 
 * and outputs data.frame xydf and file xy<n>.csv.
 * or use ./ompCD -n <n> -d 1 to get same files with numPt=<n>.
 */

#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>
#include <math.h>
#include <mathimf.h>
#include <omp.h>

#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define MILLION 1000000L

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MinRows 4
#define MaxColumns 2
#define Data(i,j) data[(j)*(*row)+(i)]
/*R uses column-major order in converting matrix to array
 */ 

uint32_t numThread, verify, debug, readfile;
uint64_t numPt, seed;
uint64_t *sdnc, *sd, *sdbf, *td;
double *x, *y, *ang, *sdSec, *tdSec, *cdSec;

typedef struct {
  double angle;
  int charge;
} Node;
Node **node;

double drand()   /*uniform distribution, (0..1]*/
{
  return (rand()+1.0)/(RAND_MAX+1.0);
}

double random_normal()
/*normal distribution, centered on 0, std dev 1 
 *mod to use drand48() not drand() above
 */
{
   return sqrt(-2*log(drand48())) * cos(2*M_PI*drand48());
}

static int cmpNode(const void *a, const void *b) {
  Node *n1, *n2;

  n1 = (Node*) a;
  n2 = (Node*) b;

  /*sort by increasing angle*/ 
  if (n1->angle < n2->angle) return -1;
  if (n1->angle > n2->angle) return 1;

  /*break tie by placing data points (charge == -1) ahead of
   *antipodal points (charge == 1), important for sd correctness
   */
  return n1->charge - n2->charge;
}


/*int main(int argc, char *argv[]) {*/
void omprXY(int *row,double *xvec,double *yvec){
  struct timeval start, stop;
  uint64_t i,j,k;
  float sec;
  int c;
  FILE *fptr, *fptr2;
  /*numPt = 100100; use -n 100100 for this*/
  seed = 1;/*opt -s*/
  numThread = omp_get_max_threads();/*opt -t*//*min for omp is 2*/


  if ( (*row) < MinRows ) {
    fprintf(stderr, "error: minimum %u rows\n", MinRows);
    exit(1);
  }
  
  numPt = (*row);
  numThread = omp_get_max_threads();
  omp_set_num_threads(numThread);

 
  srand48(seed);

  x = (double*) malloc(numPt * sizeof(double));
  y = (double*) malloc(numPt * sizeof(double));

 
  //else{
    printf("\nUse random points; numPt = %lu.",numPt);

    /*generate random points*/
    for (i=0; i < numPt; i++) {
      x[i] = random_normal();
      y[i] = random_normal();
    }
  //}
  ///*end else Use random points*/

  for (i = 0; i < numPt; i++){
     xvec[i] = (double)x[i];
     yvec[i] = (double)y[i];
  }

  printf("\nResults:\n");
 

  /*fptr xy.csv used for testing readfile*/
  //if (debug){
    fptr = fopen("xy.csv","w");
    if (fptr != NULL){
      printf("File xy.csv created.\n");
      fprintf(fptr,"x, y\n");
      for (i = 0; i < numPt; i++){
        fprintf(fptr, "%5.2f, %5.2f\n",x[i],y[i]);
      }
    }
    else{
      printf("Failure on xy.csv creation.\n");
    }
    fclose(fptr);
  //}


  printf("\n");/*to newline after the run*/

  free(y);
  free(x);
  return;
}
