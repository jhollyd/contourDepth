/*omprSD.c is for R, computes parallelized simplicial depth,
 * via method omprSD(), outputs to R vector sdepth, 
 * takes input from data frame named xydf in R 
 * that has 2 D data for x, y.
 * Note: In R use xydf <- read.csv(file="xy.csv") to make df. 
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
uint64_t *sdnc, *sd;
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

int cmpDoub (const void * a, const void * b) {
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
 *counting the number of triangles that do not contain the origin
 */
void sdRR(void) {
  uint64_t i, j, k, h, first;
  double angle, antipodal;
  uint32_t tid;

tid=0;
#pragma omp parallel for private(tid,j,k,h,first,angle,antipodal)/*unc*/
  for (i = 0; i < numPt; i++) { /*L1 - find simplicial depth of point i*/
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

    /*L3 - find the first data point*/
    for (first = 0; node[tid][first].charge == 1; first++)
      ;

    if (node[tid][first].angle > 180) {/*highest y valued point sdRR*/
      sdnc[i] = (numPt-1)*(numPt-2)*(numPt-3)/6; /*number of non-containing triangles, was 0 when point was origin*/
      sd[i] = (numPt-1)*(numPt-2)/2; /*adjusted (n*n-1*n-2)/6-sd*/
      continue;
    }

    /*angle of the first antipodal point*/
    antipodal = node[tid][first].angle < 180 ? node[tid][first].angle + 180 :
      node[tid][first].angle - 180;

    /*L4 - find h: the number of data points in the first half plane (HP1)*/
    for (j = first + 1, h = 0;
	 j < (numPt - 1) * 2 && node[tid][j].angle < antipodal;
	 j++)
      if (node[tid][j].charge == -1) /*a data point*/
	h++;

    /*number of non-containing triangles in HP1*/
    sdnc[i] = h > 1 ? h * (h - 1) / 2 : 0;

    /*L5 - rotating the half plane in increasing degrees, counterclockwise*/
    for (j = first + 1; j < (numPt - 1) * 2; j++) {
      /*update the number of data points in HP1*/
      h += node[tid][j].charge;
      
      /*if the incoming point is a data point
       *calculate the number of non-containing triangles
       */
      if (node[tid][j].charge == -1)
	sdnc[i] += h > 1 ? h * (h - 1) / 2 : 0;
    }

    /*simplicial depth: (number of all triangles) -
     *(number of non-containing triangles), but must adjust to incl origin
     * (no non-containing triangles, since we add back the origin)
     * use the #Tri(n) to match bruteForce
     */
    sd[i] = (numPt) * (numPt - 1) * (numPt - 2) / 6  -  sdnc[i];
  }
}

void omprSD( int *row, int *col, double *data, double *depth){
  struct timeval start, stop;
  uint64_t i,j;
  uint64_t maxsd;
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
  sd = (uint64_t*) malloc(numPt * sizeof(uint64_t));

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

  printf("Results:\n");
 
  /*run sdRR() and time it*/
  gettimeofday(&start, NULL);
  sdRR();
  gettimeofday(&stop, NULL);
  sec = (stop.tv_sec - start.tv_sec) +
    (stop.tv_usec - start.tv_usec) / (float) MILLION;
  printf("Rousseeuw and Ruts simplicial depth, %lu points, %u threads, time %.2f.\n",numPt, numThread, sec);

  for (i = 0; i < numPt; i++){
     depth[i] = (double)sd[i];
   }

  /*find max SD */
  qsort(sd, numPt, sizeof(double), cmpDoub);
  maxsd = sd[numPt - 1]; 
  fprintf(stderr, "  Found simplicial depth max=%lu.\n",maxsd);

  printf("\n");/*newline after the run*/

  for (i = 0; i < numThread; i++)
    free(node[i]);
  free(node);
  free(sdnc);
  free(sd);
  free(y);
  free(x);
  return;
}
