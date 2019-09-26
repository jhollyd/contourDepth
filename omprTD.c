/*omprTD.c is for R, computes parallelized Tukey depth,
 * via method omprTD(), outputs to R vector tdepth, 
 * takes input from data frame named xydf in R 
 * that has 2 D data for x, y.
 * Note: In R use xydf <- read.csv(file="xy.csv") to make xydf. 
 *   A file of random Gaussian data of length numPt can be generated.
 *   Outside R, type:
 *   ./ompCD -d 1 to get xy.csv of default size numPt=20100, or 
 *   ./ompCD -n xxxxx -d 1 for numPt = xxxxx.
 */

#include <R.h>/*for R*/
#include <stdint.h>
#include <omp.h>

#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define MILLION 1000000L

/*for R: 3 defines*/
#define MinRows 4 
#define MaxColumns 2 
#define Data(i,j) data[(j) * numPt + (i)]
/*R uses column-major order in converting matrix to array*/

uint32_t numThread;
uint64_t numPt;
uint64_t *td, *sdnc;
double *x, *y;

typedef struct {
  double angle;
  int charge;
} Node;
static Node **node;

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

int cmpDoub2 (const void * a, const void * b) {
  const double *p, *q;
  p = (const double*)a;
  q = (const double*)b;
  if (*p > *q) return 1;
  else if (*p < *q) return -1;
  else return 0;
  /*return ( *(double*)a - *(double*)b );*/
}

/*Rousseeuw and Ruts, Bivariate location depth, 1996
 *moving each point to the origin
 *counting the number of points in the halfplane 
 */
void Tukey(void) {
  uint64_t i, j, k, h, first, tdsp, q;
  double angle, antipodal;
  uint32_t tid;

tid=0;
#pragma omp parallel for private(tid,j,k,h,first,angle,antipodal,tdsp,q)
  for (i = 0; i < numPt; i++) { /*L1 - find Tukey depth of point i*/
    tid = omp_get_thread_num();
    for (j = k = 0; j < numPt; j++) {/*L2*/ 
      if (j == i) /*skip self*/ 
	continue;

      /*move point i to the origin, calculate angles of other points*/
      angle = atan2d(y[j] - y[i], x[j] - x[i]);
      /*-180 <= angle <= 180*/
      angle = angle >= 0 ? angle : angle + 360;
      /*0 <= angle <= 360*/

      /*data point*/
      node[tid][k].angle = angle;
      node[tid][k].charge = -1;
      k++;
    
      /*antipodal point*/
      node[tid][k].angle = angle < 180 ? angle + 180 : angle - 180;
      node[tid][k].charge = 1;
      k++;
    }
   
    qsort((void *)node[tid], (numPt - 1) * 2, sizeof(Node), cmpNode);

    /*L3 - find first*/
    for (first = 0; node[tid][first].charge == 1; first++)
      ;

    if (node[tid][first].angle > 180) {/*highest y valued point */
      td[i] = 1; /*first halfplane initial td is final td for this point*/
      continue;
    }

    /*compute angle of the first antipodal point*/
    antipodal = node[tid][first].angle < 180 ? node[tid][first].angle + 180:      node[tid][first].angle - 180;

    /*L4 - find h: the number of data points in the first half plane (HP1),
     * including the first data point
     */
    for (j = first, h = 0;
	 j < (numPt - 1) * 2 && node[tid][j].angle < antipodal;
	 j++) /*L4*/
      if (node[tid][j].charge == -1) /*a data point*/
	h++;

    /*initialize the Tukey depth, high value:
     *RR used: td[i]=min(h,(numPt-1)-h);
     *use numPt-1.
     *Pseudocode: calc tdsp, min(h, n-h) so spokedepth-f(n,h)
     *td[i]=min(n-1,tdsp) so initial value = n-1
     *min(min(td[i],h),numPt-h)
     */
    td[i] = numPt - 1;/*numPt or (numPt-1) ok here*/

    /*Rotate the halfplane in increasing angles - counterclockwise.
     *The line of the halfplane does not contain any data points or
     *any antipodal points.  The line is conceptually positioned slightly
     *clockwise from node[i] -- between node[i-1] and node[i].
     *In one iteration of the for loop, the line rotates past node[i],
     *and node[i+1], ... node[i+k], if these nodes have the same
     *angle in degrees as node[i].  The new line is positioned between
     *node[i+k] and node[i+k+1].
     */

    /*R&R wheel goes through 2*(numPt-1) positions after HP1, each just 
     *clockwise of the line to a specified angle, but we go through
     *only 180 degrees, covering HPs in other halfplanes by deltas.
     */

    /*L5 again from first, not first + 1, and stopping at 180 degrees*/
    for (j = first; j < (numPt - 1) * 2  && node[tid][j].angle <=180; j++) {
      /*update the number of data points in HP1*/
      h += node[tid][j].charge;
     
      /*L6 for angles that are equal to the one selected*/
      q = 0;
      //while(j+1<(numPt-1)*2 && node[tid][j+1].angle==node[tid][j].angle){
        //h += node[tid][j].charge;
        //q++;
      //}/*end L6, was hang-prone in R and unlikely*/

      /*calculate tdsp of angle j for this i out of 2*(numPt-1) nodes*/
      if (node[tid][j].charge == -1){/*a data point angle*/
        if ( h >= (numPt-1)/2 )
          tdsp = numPt - 1 - h + q;
        else
          tdsp = h + 1 + q;
      }
      else {/*if antipodal angle*/
        if ( h >= (numPt-1)/2 ) 
          tdsp = numPt - h + q;
        else
          tdsp = h + q;
      }

      /*use tdsp to update the running Tukey depth*/
      td[i] = min( td[i], tdsp );

    }/*end of loop on j*/
  }/*end of for loop on i*/
  /*RR action for a data point: nf++*/
  /*RR action for an antipodal, beta_i: h = nf-i, i++, goes to < antipodal*/
}

void omprTD( int *row, int *col, double *data, double *depth){
  struct timeval start, stop;
  uint64_t i,j;
  uint64_t maxtd;
  float sec;
  FILE *fptr, *fptr2;
  numThread = omp_get_max_threads();/*min for omp is 2*/

  if ( (*row) < MinRows ) {
    fprintf(stderr, "error: minimum %u rows\n", MinRows);
    exit(1);
  }
  if ( (*col) > MaxColumns ) {
    fprintf(stderr, "error: maximum %u columns\n", MaxColumns);
    exit(1);
  }
  numPt = (*row);
  numThread = omp_get_max_threads();
  omp_set_num_threads(numThread);
 
  x = (double*) malloc(numPt * sizeof(double));
  y = (double*) malloc(numPt * sizeof(double));
  sdnc  = (uint64_t*) malloc(numPt * sizeof(uint64_t));
  //sd = (uint64_t*) malloc(numPt * sizeof(uint64_t));
  td = (uint64_t*) malloc(numPt * sizeof(uint64_t));

  /*pull in the x,y data*/
  for (j = 0; j < (*col); j++){
    for (i = 0; i < numPt; i++){
      x[i] = Data(i,0);
      y[i] = Data(i,1);
    }
  }

  node = (Node**) malloc(sizeof(Node*) * numThread);
  /*for each thread, malloc for 2*(numPt-1) nodes and store in node*/
  for (i = 0; i < numThread; i++)
    node[i] = (Node*) malloc(sizeof(Node) * (numPt - 1) * 2);

  printf("\nResults:\n");
 
  /*run Tukey() and time it*/
  gettimeofday(&start, NULL);
  Tukey();
  gettimeofday(&stop, NULL);
  sec = (stop.tv_sec - start.tv_sec) +
    (stop.tv_usec - start.tv_usec) / (float) MILLION;
  printf("Tukey depth, %lu points, %u threads, time %.2f.\n",numPt, numThread, sec);

  for (i = 0; i < numPt; i++){
     depth[i] = (double)td[i];
   }

  /*find max TD */
  qsort(td, numPt, sizeof(double), cmpDoub2);
  maxtd = td[numPt - 1]; 
  fprintf(stderr, "  Found Tukey depth max=%lu.\n",maxtd);

  printf("\n");/*newline after the run*/

  for (i = 0; i < numThread; i++)
    free(node[i]);
  free(node);
  free(td);
  free(sdnc);
  free(y);
  free(x);
  return;
}
