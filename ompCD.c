/*ompCD.c computes parallelized simplicial and Tukey depths, 
 * and has main, outputs to files, input from file.
 * Use options to see all results.
 */

/*for R: #include <R.h>
 */ 
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

/*for R: #define MinRows 4, MaxColumns 2, Data(i,j) data[(j)*(*row)+(i)]
 *R uses column-major order in converting matrix to array
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

/*point p is inside the closed triangle of (i, j, k)*/
uint32_t insideTriangle(uint32_t p, uint32_t i, uint32_t j, uint32_t k) {
  uint32_t count = 0;
  double sign; 

  /*crossProduct(P_i, P_j)*/
  sign = (x[i] - x[p]) * (y[j] - y[p]) - (y[i] - y[p]) * (x[j] - x[p]);
  if (sign == 0.0 && (x[i] - x[p]) * (x[j] - x[p]) <= 0.0 &&
      (y[i] - y[p]) * (y[j] - y[p]) <= 0.0)
    return 1; /*the origin is on the line segment [P_i, P_j]*/
  count += sign < 0 ? -1 : 1;

  /*crossProduct(P_j, P_k)*/
  sign = (x[j] - x[p]) * (y[k] - y[p]) - (y[j] - y[p]) * (x[k] - x[p]);
  if (sign == 0.0 && (x[j] - x[p]) * (x[k] - x[p]) <= 0.0 &&
      (y[j] - y[p]) * (y[k] - y[p]) <= 0.0)
    return 1; //the origin is on the line segment [P_j, P_k]
  count += sign < 0 ? -1 : 1;

  /*crossProduct(P_k, P_i)*/
  sign = (x[k] - x[p]) * (y[i] - y[p]) - (y[k] - y[p]) * (x[i] - x[p]);
  if (sign == 0.0 && (x[i] - x[p]) * (x[k] - x[p]) <= 0.0 &&
      (y[i] - y[p]) * (y[k] - y[p]) <= 0.0)
    return 1; /*the origin is on the line segment [P_i, P_k]*/
  count += sign < 0 ? -1 : 1;

  /*if all three cross products have the same sign
   *the origin is inside the triangle
   */
  return count == 3 || count == -3;
}

/*brute force the simplicial depth*/
/*form each i-j-k triangle; is origin p in or on or out?*/
void sdBruteForce(void) {
  uint32_t p, i, j, k;

#pragma omp parallel for private(i,j,k)
  for (p = 0; p < numPt; p++) {
    /*sdbf[p] = 0; not = 0 with origin as a point*/
    sdbf[p] = (numPt - 1) *  (numPt - 2) / 2;/*add this for origin included*/
    for (i = 0; i < numPt - 2; i++) {
      if (i == p)
	continue;
      for (j = i + 1; j < numPt - 1; j++) {
	if (j == p)
	  continue;
	for (k = j + 1; k < numPt; k++) {
	  if (k == p)
	    continue;
	  sdbf[p] += insideTriangle(p, i, j, k);
	}
      }
    }
  }
}

void Tukey(void) {
  uint64_t i, j, k, h, first;
  uint64_t tdsp, q;
  double angle, antipodal;
  uint32_t tid;

tid=0;
#pragma omp parallel for private(tid,j,k,h,first,angle,antipodal,tdsp,q)
  for (i = 0; i < numPt; i++) { //L1 - find Tukey depth of point i
    tid = omp_get_thread_num();
    for (j = k = 0; j < numPt; j++) {/*L2*/ 
      if (j == i){/*skip self*/
        continue;
      }

      /*move point i to the origin, calculate angles of other points*/
      angle = atan2d(y[j] - y[i], x[j] - x[i]);
      /*-180 <= angle <= 180*/
      angle = angle >= 0 ? angle : angle + 360;
      /*0 <= angle <= 360 CCW*/ 

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

    if (node[tid][first].angle > 180) {/*affects point at top of graph*/
      td[i] = 1; /*first halfplane initial td is final td*/
      continue;  
    }

    /*compute angle of the first antipodal point*/
    antipodal = node[tid][first].angle < 180 ? node[tid][first].angle + 180 :  node[tid][first].angle - 180;

    /*L4 - count data points in the first half plane (HP1), including
     *the first data point and up to just less than its antipodal
     */ 
    for (j = first, h = 0;
         j < (numPt - 1) * 2 && node[tid][j].angle < antipodal;
         j++){/*L4*/ 
      if (node[tid][j].charge == -1) {/*a data point*/
        h++;
      }
    }

    /*initialize the Tukey depth:
     *RR used: td2[i]=min(h,(numPt-1)-h);
     *use numPt - 1.
     *pseudocode: calc tdsp, min(h,n-h)
     *so spokedepth=f(n,h}
     *td[i]=min(n-1,tdsp)
     *so initial value=n-1
     *min(min(td[i],h),numPt-h) 
     */ 
     td[i] = numPt - 1; /**** numPt or (numPt - 1) ok here ****/

    /*Rotate the half plane in increasing degrees -- counterclockwise.
     *The line of the half plane does not contain any data points or
     *antipodal points.
     *This line is conceptually positioned slightly clockwise from
     *node[i] -- between node[i - 1] and node[i].
     *In one iteration of the for loop, the line rotates past node[i],
     *and node[i + 1], ... node[i + k], if these nodes have the same
     *degree as node[i].
     *The new line is positioned between node[i + k] and node[i + k + 1].
     */    
    
    /*Wheel goes through 2*(numPt-1) positions after HP1
     *each is just CW of the line to a specified angle
     *but we only go to 180 deg, covering HPs in other half planes by deltas
     */

    /*L5*/
    for (j = first; j <(numPt - 1) * 2 && node[tid][j].angle <= 180; j++) {
       h += node[tid][j].charge;

       /*L6 for angles that are equal to the one selected*/
       q = 0;
       while (j + 1 < (numPt - 1) * 2 &&
              node[tid][j + 1].angle == node[tid][j].angle){
         h += node[tid][j].charge;
         q++; 
       }/*end L6*/

       /*calc tdsp of angle j for this i out of 2*(numPt-1) nodes*/
       if (node[tid][j].charge == -1){/*data point angle*/
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

int main(int argc, char *argv[]) {
  struct timeval start, stop;
  uint64_t i,j,k;
  float sec;
  int c;
  FILE *fptr, *fptr2;
  /*numPt = 100100; use -n 100100 for this*/
  numPt = 20100;/*opt -n*/ 
  seed = 1;/*opt -s*/
  numThread = omp_get_max_threads();/*opt -t*//*min for omp is 2*/

  verify = 0;/*opt -v*/
  debug = 0;/*opt -d*/
  readfile = 0;/*opt -r <numPt>*/
  const int numCase = 11;
  int numm = 99;

  while ((c = getopt(argc, argv, "n:s:t:v:d:r:")) != -1) {
    switch (c) {
    case 'n': sscanf(optarg, "%lu", &numPt);    break; /*number of data pts*/
    case 's': sscanf(optarg, "%lu", &seed);     break; /*seed srand48*/
    case 't': sscanf(optarg, "%u", &numThread); break; /*number of threads*/
    case 'v': sscanf(optarg, "%u", &verify);    break; /*show bruteForceSD*/
    case 'd': sscanf(optarg, "%u", &debug);     break; /*outfile xy.csv*/
    case 'r': sscanf(optarg, "%u", &readfile);  break; /*number of readfile data points*/
    default: break;
    }
  }
 
  if (readfile > 0){
    numPt = readfile;
  }
 
  /* adjust numThread; not doing perfpar serial run here*/
  if (numThread > omp_get_max_threads()){
    printf("Requested numThread, %u > max available, %u, so using max available.\n", numThread,omp_get_max_threads()); 
     numThread = omp_get_max_threads();
  } 
  else if (numThread == 1) {
     printf("Requested numThread, %u < min allowable for omp, =2, so using 2.  Comment out omp to make serial run. \n",numThread);
     numThread = 2;
     omp_set_num_threads(numThread);
  }
  else {
     //printf("Running with requested numThread, %u.\n",numThread);
     omp_set_num_threads(numThread);
  }
 
  srand48(seed);

  x = (double*) malloc(numPt * sizeof(double));
  y = (double*) malloc(numPt * sizeof(double));
  sdnc  = (uint64_t*) malloc(numPt * sizeof(uint64_t));
  sdbf = (uint64_t*) malloc(numPt * sizeof(uint64_t));/*bruteForce sd*/
  sd = (uint64_t*) malloc(numPt * sizeof(uint64_t));
  /* malloc for node is below */
  td  = (uint64_t*) malloc(numPt * sizeof(uint64_t));
  ang = (double*) malloc(numPt * sizeof(double));/*angle fr (0,0)*/

  /*set the performance graph values for numPt; numCase is set above*/  
  uint64_t numPoint[numCase] = {100,10100,20100,30100,40100,50100,60100,70100,80100,90100,100100}; 
  sdSec = (double*) malloc(numCase * sizeof(double));
  tdSec = (double*) malloc(numCase * sizeof(double));
  cdSec = (double*) malloc(numCase * sizeof(double));

  char chunk[64];
  size_t len = sizeof(chunk);
  chunk[len-1] = '\0';
  int row = 0;
  int col = 2;
  int num = 0;
  float v1, v2;
  int count = 0;
  int count2 = 0;
 
  if (readfile){
    printf("\nUse read from input file xy.csv. ");/*this is readfile -r <numPt>*/
    fptr2 = fopen("xy.csv","r");
    if (fptr2 != NULL) {
      printf("File found. ");

      while(fgets(chunk, len, fptr2)!=NULL){
        num = sscanf(chunk,"%f,%f",&v1,&v2);
        if (num == 2){
          if (count2 == 0){
            x[row] = v1;
            y[row] = v2; 
          } 
          row++;
          if (row > numPt){
            printf("  Read interrupted at numPt = %lu row = %d rows of values.\n",numPt,row);
            count2++; 
            continue;/*this repeats because fgets can't be stopped*/
          }
        }
        else{
          count++;
          /*printf("A non 2 (num) value read occurred, num = %d; count = %d.\n",num,count);*/
          if (count > 1){
            printf("  Too few in xy.csv, count = %d; stop reading.\n",count);
            continue;
          } 
          else{ 
            printf("  Header line skipped in xy.csv, count = %d; read on.\n",count);
          }
        }
      }/*end while*/

      fclose(fptr2);

      printf("  File xy.csv: has %d rows; numPt=%lu.",row,numPt);
      /* if numPt > row reset numPt*/
      if (numPt > row){
        numPt = row;
        printf(" Reset numPt = row so numPt = %lu.",numPt);
      }
      else{
        printf(" No reset of numPt required.");
      }
 
      for (i=0; i < numPt; i++){
        ang[i] = atan2d( y[i],x[i] );
        /*-180 <= angle <= 180*/
        ang[i] = ang[i] >= 0 ? ang[i] : ang[i] + 360;
        /*0 <= angle <= 360*/
      } 
    }/*fileptr not NULL*/

    else {
      printf("Failed to find file xy.csv.\n");
      fclose(fptr2);/*?*/
    }
  }/*end if Use readfile*/

  else{
    printf("\nUse random points; numPt = %lu.",numPt);

    /*generate random points*/
    for (i=0; i < numPt; i++) {
      x[i] = random_normal();
      y[i] = random_normal();
      ang[i] = atan2d( y[i], x[i] );
      /*-180 <= angle <= 180*/
      ang[i] = ang[i] >= 0 ? ang[i] : ang[i] + 360;
      /*0 <= angle <= 360*/
    }
  }/*end else Use random points*/

  /* for R */
  /*qsort(ang, numPt, sizeof(double), cmpDoub);*/
  
  node = (Node**) malloc(sizeof(Node*) * numThread);
  /*for each thread, malloc for 2*(numPt-1) nodes and store in node*/
  for (i = 0; i < numThread; i++)
    node[i] = (Node*) malloc(sizeof(Node) * (numPt - 1) * 2);

  printf("\nResults:\n");
 
  /* to store times */
  if (numPt % 10000 == 100) {
    numm = (int)( ( numPt - 100) / 10000);
    /*printf("Time will be stored at numm=%d.\n",numm);*/
  }

  gettimeofday(&start, NULL);
  sdRR();
  gettimeofday(&stop, NULL);
  sec = (stop.tv_sec - start.tv_sec) +
    (stop.tv_usec - start.tv_usec) / (float) MILLION;
  printf("Rousseeuw and Ruts simplicial depth, %lu points, %u threads, time %.2f\n",numPt, numThread, sec);
  if(!(numm == 99)) 
    sdSec[numm] = sec;
 
  gettimeofday(&start, NULL);
  Tukey();
  gettimeofday(&stop, NULL);
  sec = (stop.tv_sec - start.tv_sec) +
    (stop.tv_usec - start.tv_usec) / (float) MILLION;
  printf("Tukey depth, %lu points, %u threads, time %.2f\n", numPt, numThread, sec);
  if(!(numm == 99)) 
    tdSec[numm] = sec;

  if (verify) {
    gettimeofday(&start, NULL);
    sdBruteForce();
    gettimeofday(&stop, NULL);
    sec = (stop.tv_sec - start.tv_sec) +
      (stop.tv_usec - start.tv_usec) / (float) MILLION;
    printf("  and brute force sd time %.2f\n", sec);

    for (i = 0; i < numPt; i++)
      if (sd[i] != sdbf[i])
	printf("for i=%lu, mismatch sd by sdRR vs sdBF %lu %lu\n", i, sd[i], sdbf[i]);
  }

  if (verify) {
	  printf("\nResults2: Verify is on, numPt is %lu, here are \n   id:  td:   x:    y:  ang:  Tri(n-1): Tri(n): sdNonConT: sd=T(n)-sNCT: sdBFConT:\n",numPt);
    for (i = 0; i < numPt; i++){
      printf("%5lu, %5lu, % 6.2f, % 6.2f, % 6.2f, %5lu, %5lu, %5lu, %5lu, %5lu\n",i,td[i],x[i],y[i],ang[i],(numPt-1)*(numPt-2)*(numPt-3)/6,(numPt)*(numPt-1)*(numPt-2)/6,sdnc[i],sd[i],sdbf[i]);
    }
  }

  /*start CD timer: sdRR, Tukey and file outputs (repeats sdRR(), Tukey())*/
  gettimeofday(&start, NULL);
  sdRR();
  Tukey();

  /* output to file as i, x, y, ang, sd, (sdBF,) td */
  fptr = fopen("xyst.csv","w");
  if (fptr != NULL){
    printf("File xyst.csv created. \n");
    if (verify){
      fprintf(fptr, "i,x,y,ang,sd,sdBF,td\n");
    }
    else{
      fprintf(fptr, "i,x,y,ang,sd,td\n");
    }
    for (i = 0; i < numPt; i++) {/*for each of randomly chosen numPt pts */
      if (verify) {
        fprintf(fptr, "%2lu, %5.2f, %5.2f, %5.2f, %2lu, %2lu, %2lu\n",i,x[i],y[i],ang[i],sd[i],sdbf[i],td[i]);
      }
      else {
        fprintf(fptr, "%2lu, %5.2f, %5.2f, %5.2f, %3lu, %3lu\n",i,x[i],y[i],ang[i],sd[i],td[i]);
      }
    }   
  }
  else {
    printf("Failure on xyst.csv creation.\n");
  }
  fclose(fptr);

  /*fptr xy.csv used for testing readfile*/
  if (debug){
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
  }

  /*fptr ixy.csv used for testing readfile*/
  if (debug){
    fptr = fopen("ixy.csv","w");
    if (fptr != NULL){
       printf("File ixy.csv created.\n");
       fprintf(fptr,"i, x, y\n");
       for (i = 0; i < numPt; i++){
         fprintf(fptr, "%3lu, %5.2f, %5.2f\n",i,x[i],y[i]);
       }
     }
     else{
       printf("Failure on ixy.csv creation.\n");
     }
     fclose(fptr);
   }

  /*stop timer for CD*/
  gettimeofday(&stop, NULL);
  sec = (stop.tv_sec - start.tv_sec) +
     (stop.tv_usec - start.tv_usec) / (float) MILLION;
  printf("RR simplicial depth, Tukey depth, file outputs, %lu points, %u threads, time %.2f\n", numPt, numThread, sec);
  if(!(numm == 99)) 
    cdSec[numm] = sec;

  /* output to file nstc.csv, third timer*/
  fptr = fopen("nstc.csv","w");
  if (fptr != NULL){
    printf("File nstc.csv created.\n");
    fprintf(fptr, "numPoint,sdSec,tdSec,cdSec\n");
    for (i = 0; i < numCase; i++){
      fprintf(fptr, "%6lu, %5.2f, %5.2f, %5.2f\n",numPoint[i],sdSec[i],tdSec[i],cdSec[i]);
    }
  }
  else {
    printf("Failure on nstc.csv creation.\n");
  }
  fclose(fptr);

  printf("\n");/*to newline after the run*/

  free(cdSec);
  free(tdSec);
  free(sdSec);
  for (i = 0; i < numThread; i++)
    free(node[i]);
  free(node);
  free(td);
  free(sdnc);
  free(sdbf);
  free(sd);
  free(y);
  free(x);
  free(ang);
  return 0;
}
