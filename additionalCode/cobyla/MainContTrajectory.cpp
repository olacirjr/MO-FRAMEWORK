#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h> 
#include <time.h> 
#include <stdio.h> 
#include <string.h> 
#include <stdlib.h> 

#include "mga_dsm.h"
#include "auxfunc.h"  
#include "ContPopul.h"  
#include "EDA.h" 
#include "AbstractTree.h"  
#include "FDA.h"  
#include "MixtureTrees.h" 
#include "GaussianModel.h" 
#include "cobyla.h"



cobyla_function calcfc;
typedef struct 
{
  int nprob;
} example_state;



#include "cobyla.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define abs(x) ((x) >= 0 ? (x) : -(x))

/*
 * Return code strings
 */
char *cobyla_rc_string[6] =
{
  "N<0 or M<0",
  "Memory allocation failed",
  "Normal return from cobyla",
  "Maximum number of function evaluations reached",
  "Rounding errors are becoming damaging",
  "User requested end of minimization"
};

static int cobylb(int *n, int *m, int *mpp, double *x, double *rhobeg,
  double *rhoend, int *iprint, int *maxfun, double *con, double *sim,
  double *simi, double *datmat, double *a, double *vsig, double *veta,
  double *sigbar, double *dx, double *w, int *iact, cobyla_function *calcfc,
  void *state);
static int trstlp(int *n, int *m, double *a, double *b, double *rho,
  double *dx, int *ifull, int *iact, double *z__, double *zdota, double *vmultc,
  double *sdirn, double *dxnew, double *vmultd);

/* ------------------------------------------------------------------------ */

int cobyla(int n, int m, double *x, double rhobeg, double rhoend, int iprint,
  int *maxfun, cobyla_function *calcfc, void *state)
{
  int icon, isim, isigb, idatm, iveta, isimi, ivsig, iwork, ia, idx, mpp, rc;
  int *iact;
  double *w;

/*
 * This subroutine minimizes an objective function F(X) subject to M
 * inequality constraints on X, where X is a vector of variables that has 
 * N components. The algorithm employs linear approximations to the 
 * objective and constraint functions, the approximations being formed by 
 * linear interpolation at N+1 points in the space of the variables. 
 * We regard these interpolation points as vertices of a simplex. The 
 * parameter RHO controls the size of the simplex and it is reduced 
 * automatically from RHOBEG to RHOEND. For each RHO the subroutine tries 
 * to achieve a good vector of variables for the current size, and then 
 * RHO is reduced until the value RHOEND is reached. Therefore RHOBEG and 
 * RHOEND should be set to reasonable initial changes to and the required 
 * accuracy in the variables respectively, but this accuracy should be 
 * viewed as a subject for experimentation because it is not guaranteed. 
 * The subroutine has an advantage over many of its competitors, however, 
 * which is that it treats each constraint individually when calculating 
 * a change to the variables, instead of lumping the constraints together 
 * into a single penalty function. The name of the subroutine is derived 
 * from the phrase Constrained Optimization BY Linear Approximations. 
 *
 * The user must set the values of N, M, RHOBEG and RHOEND, and must 
 * provide an initial vector of variables in X. Further, the value of 
 * IPRINT should be set to 0, 1, 2 or 3, which controls the amount of 
 * printing during the calculation. Specifically, there is no output if 
 * IPRINT=0 and there is output only at the end of the calculation if 
 * IPRINT=1. Otherwise each new value of RHO and SIGMA is printed. 
 * Further, the vector of variables and some function information are 
 * given either when RHO is reduced or when each new value of F(X) is 
 * computed in the cases IPRINT=2 or IPRINT=3 respectively. Here SIGMA 
 * is a penalty parameter, it being assumed that a change to X is an 
 * improvement if it reduces the merit function 
 *      F(X)+SIGMA*MAX(0.0,-C1(X),-C2(X),...,-CM(X)), 
 * where C1,C2,...,CM denote the constraint functions that should become 
 * nonnegative eventually, at least to the precision of RHOEND. In the 
 * printed output the displayed term that is multiplied by SIGMA is 
 * called MAXCV, which stands for 'MAXimum Constraint Violation'. The 
 * argument MAXFUN is an int variable that must be set by the user to a 
 * limit on the number of calls of CALCFC, the purpose of this routine being 
 * given below. The value of MAXFUN will be altered to the number of calls 
 * of CALCFC that are made. The arguments W and IACT provide real and 
 * int arrays that are used as working space. Their lengths must be at 
 * least N*(3*N+2*M+11)+4*M+6 and M+1 respectively. 
 *
 * In order to define the objective and constraint functions, we require 
 * a subroutine that has the name and arguments 
 *      SUBROUTINE CALCFC (N,M,X,F,CON) 
 *      DIMENSION X(*),CON(*)  . 
 * The values of N and M are fixed and have been defined already, while 
 * X is now the current vector of variables. The subroutine should return 
 * the objective and constraint functions at X in F and CON(1),CON(2), 
 * ...,CON(M). Note that we are trying to adjust X so that F(X) is as 
 * small as possible subject to the constraint functions being nonnegative. 
 *
 * Partition the working space array W to provide the storage that is needed 
 * for the main calculation.
 */

  if (n == 0)
  {
    if (iprint>=1) fprintf(stderr, "cobyla: N==0.\n");
    *maxfun = 0;
    return 0;
  }

  if (n < 0 || m < 0)
  {
    if (iprint>=1) fprintf(stderr, "cobyla: N<0 or M<0.\n");
    *maxfun = 0;
    return -2;
  }

  /* workspace allocation */
  w = (double*) malloc((n*(3*n+2*m+11)+4*m+6)*sizeof(*w));
  if (w == NULL)
  {
    if (iprint>=1) fprintf(stderr, "cobyla: memory allocation error.\n");
    *maxfun = 0;
    return -1;
  }
  iact = (int*) malloc((m+1)*sizeof(*iact));
  if (iact == NULL)
  {
    if (iprint>=1) fprintf(stderr, "cobyla: memory allocation error.\n");
    free(w);
    *maxfun = 0;
    return -1;
  }
  
  /* Parameter adjustments */
  --iact;
  --w;
  --x;

  /* Function Body */
  mpp = m + 2;
  icon = 1;
  isim = icon + mpp;
  isimi = isim + n * n + n;
  idatm = isimi + n * n;
  ia = idatm + n * mpp + mpp;
  ivsig = ia + m * n + n;
  iveta = ivsig + n;
  isigb = iveta + n;
  idx = isigb + n;
  iwork = idx + n;
  rc = cobylb(&n, &m, &mpp, &x[1], &rhobeg, &rhoend, &iprint, maxfun,
      &w[icon], &w[isim], &w[isimi], &w[idatm], &w[ia], &w[ivsig], &w[iveta],
      &w[isigb], &w[idx], &w[iwork], &iact[1], calcfc, state);

  /* Parameter adjustments (reverse) */
  ++iact;
  ++w;

  free(w);
  free(iact);
  
  return rc;
} /* cobyla */

/* ------------------------------------------------------------------------- */
int cobylb(int *n, int *m, int *mpp, double 
    *x, double *rhobeg, double *rhoend, int *iprint, int *
    maxfun, double *con, double *sim, double *simi, 
    double *datmat, double *a, double *vsig, double *veta,
     double *sigbar, double *dx, double *w, int *iact, cobyla_function *calcfc,
     void *state)
{
  /* System generated locals */
  int sim_dim1, sim_offset, simi_dim1, simi_offset, datmat_dim1, 
      datmat_offset, a_dim1, a_offset, i__1, i__2, i__3;
  double d__1, d__2;

  /* Local variables */
  double alpha, delta, denom, tempa, barmu;
  double beta, cmin = 0.0, cmax = 0.0;
  double cvmaxm, dxsign, prerem = 0.0;
  double edgmax, pareta, prerec = 0.0, phimin, parsig = 0.0;
  double gamma;
  double phi, rho, sum = 0.0;
  double ratio, vmold, parmu, error, vmnew;
  double resmax, cvmaxp;
  double resnew, trured;
  double temp, wsig, f;
  double weta;
  int i__, j, k, l;
  int idxnew;
  int iflag = 0;
  int iptemp;
  int isdirn, nfvals, izdota;
  int ivmc;
  int ivmd;
  int mp, np, iz, ibrnch;
  int nbest, ifull, iptem, jdrop;
  int rc = 0;

/* Set the initial values of some parameters. The last column of SIM holds */
/* the optimal vertex of the current simplex, and the preceding N columns */
/* hold the displacements from the optimal vertex to the other vertices. */
/* Further, SIMI holds the inverse of the matrix that is contained in the */
/* first N columns of SIM. */

  /* Parameter adjustments */
  a_dim1 = *n;
  a_offset = 1 + a_dim1 * 1;
  a -= a_offset;
  simi_dim1 = *n;
  simi_offset = 1 + simi_dim1 * 1;
  simi -= simi_offset;
  sim_dim1 = *n;
  sim_offset = 1 + sim_dim1 * 1;
  sim -= sim_offset;
  datmat_dim1 = *mpp;
  datmat_offset = 1 + datmat_dim1 * 1;
  datmat -= datmat_offset;
  --x;
  --con;
  --vsig;
  --veta;
  --sigbar;
  --dx;
  --w;
  --iact;

  /* Function Body */
  iptem = min(*n,4);
  iptemp = iptem + 1;
  np = *n + 1;
  mp = *m + 1;
  alpha = .25;
  beta = 2.1;
  gamma = .5;
  delta = 1.1;
  rho = *rhobeg;
  parmu = 0.;
  if (*iprint >= 2) {
    fprintf(stderr,
      "cobyla: the initial value of RHO is %12.6E and PARMU is set to zero.\n",
      rho);
  }
  nfvals = 0;
  temp = 1. / rho;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sim[i__ + np * sim_dim1] = x[i__];
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
      sim[i__ + j * sim_dim1] = 0.;
      simi[i__ + j * simi_dim1] = 0.;
    }
    sim[i__ + i__ * sim_dim1] = rho;
    simi[i__ + i__ * simi_dim1] = temp;
  }
  jdrop = np;
  ibrnch = 0;

/* Make the next call of the user-supplied subroutine CALCFC. These */
/* instructions are also used for calling CALCFC during the iterations of */
/* the algorithm. */

L40:
  if (nfvals >= *maxfun && nfvals > 0) {
    if (*iprint >= 1) {
      fprintf(stderr,
        "cobyla: maximum number of function evaluations reach.\n");
    }
    rc = 1;
    goto L600;
  }
  ++nfvals;
  if (calcfc(*n, *m, &x[1], &f, &con[1], state))
  {
    if (*iprint >= 1) {
      fprintf(stderr, "cobyla: user requested end of minimization.\n");
    }
    rc = 3;
    goto L600;
  }
  resmax = 0.;
  if (*m > 0) {
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
      d__1 = resmax, d__2 = -con[k];
      resmax = max(d__1,d__2);
    }
  }
  if (nfvals == *iprint - 1 || *iprint == 3) {
    fprintf(stderr, "cobyla: NFVALS = %4d, F =%13.6E, MAXCV =%13.6E\n",
      nfvals, f, resmax);
    i__1 = iptem;
    fprintf(stderr, "cobyla: X =");
    for (i__ = 1; i__ <= i__1; ++i__) {
      if (i__>1) fprintf(stderr, "  ");
      fprintf(stderr, "%13.6E", x[i__]);
    }
    if (iptem < *n) {
      i__1 = *n;
      for (i__ = iptemp; i__ <= i__1; ++i__) {
        if (!((i__-1) % 4)) fprintf(stderr, "\ncobyla:  ");
        fprintf(stderr, "%15.6E", x[i__]);
      }
    }
    fprintf(stderr, "\n");
  }
  con[mp] = f;
  con[*mpp] = resmax;
  if (ibrnch == 1) {
    goto L440;
  }

/* Set the recently calculated function values in a column of DATMAT. This */
/* array has a column for each vertex of the current simplex, the entries of */
/* each column being the values of the constraint functions (if any) */
/* followed by the objective function and the greatest constraint violation */
/* at the vertex. */

  i__1 = *mpp;
  for (k = 1; k <= i__1; ++k) {
    datmat[k + jdrop * datmat_dim1] = con[k];
  }
  if (nfvals > np) {
    goto L130;
  }

/* Exchange the new vertex of the initial simplex with the optimal vertex if */
/* necessary. Then, if the initial simplex is not complete, pick its next */
/* vertex and calculate the function values there. */

  if (jdrop <= *n) {
    if (datmat[mp + np * datmat_dim1] <= f) {
      x[jdrop] = sim[jdrop + np * sim_dim1];
    } else {
      sim[jdrop + np * sim_dim1] = x[jdrop];
      i__1 = *mpp;
      for (k = 1; k <= i__1; ++k) {
        datmat[k + jdrop * datmat_dim1] = datmat[k + np * datmat_dim1]
            ;
        datmat[k + np * datmat_dim1] = con[k];
      }
      i__1 = jdrop;
      for (k = 1; k <= i__1; ++k) {
        sim[jdrop + k * sim_dim1] = -rho;
        temp = 0.f;
        i__2 = jdrop;
        for (i__ = k; i__ <= i__2; ++i__) {
          temp -= simi[i__ + k * simi_dim1];
        }
        simi[jdrop + k * simi_dim1] = temp;
      }
    }
  }
  if (nfvals <= *n) {
    jdrop = nfvals;
    x[jdrop] += rho;
    goto L40;
  }
L130:
  ibrnch = 1;

/* Identify the optimal vertex of the current simplex. */

L140:
  phimin = datmat[mp + np * datmat_dim1] + parmu * datmat[*mpp + np * 
      datmat_dim1];
  nbest = np;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    temp = datmat[mp + j * datmat_dim1] + parmu * datmat[*mpp + j * 
        datmat_dim1];
    if (temp < phimin) {
      nbest = j;
      phimin = temp;
    } else if (temp == phimin && parmu == 0.) {
      if (datmat[*mpp + j * datmat_dim1] < datmat[*mpp + nbest * 
          datmat_dim1]) {
        nbest = j;
      }
    }
  }

/* Switch the best vertex into pole position if it is not there already, */
/* and also update SIM, SIMI and DATMAT. */

  if (nbest <= *n) {
    i__1 = *mpp;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = datmat[i__ + np * datmat_dim1];
      datmat[i__ + np * datmat_dim1] = datmat[i__ + nbest * datmat_dim1]
          ;
      datmat[i__ + nbest * datmat_dim1] = temp;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = sim[i__ + nbest * sim_dim1];
      sim[i__ + nbest * sim_dim1] = 0.;
      sim[i__ + np * sim_dim1] += temp;
      tempa = 0.;
      i__2 = *n;
      for (k = 1; k <= i__2; ++k) {
        sim[i__ + k * sim_dim1] -= temp;
        tempa -= simi[k + i__ * simi_dim1];
      }
      simi[nbest + i__ * simi_dim1] = tempa;
    }
  }

/* Make an error return if SIGI is a poor approximation to the inverse of */
/* the leading N by N submatrix of SIG. */

  error = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
      temp = 0.;
      if (i__ == j) {
        temp += -1.;
      }
      i__3 = *n;
      for (k = 1; k <= i__3; ++k) {
        temp += simi[i__ + k * simi_dim1] * sim[k + j * sim_dim1];
      }
      d__1 = error, d__2 = abs(temp);
      error = max(d__1,d__2);
    }
  }
  if (error > .1) {
    if (*iprint >= 1) {
      fprintf(stderr, "cobyla: rounding errors are becoming damaging.\n");
    }
    rc = 2;
    goto L600;
  }

/* Calculate the coefficients of the linear approximations to the objective */
/* and constraint functions, placing minus the objective function gradient */
/* after the constraint gradients in the array A. The vector W is used for */
/* working space. */

  i__2 = mp;
  for (k = 1; k <= i__2; ++k) {
    con[k] = -datmat[k + np * datmat_dim1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      w[j] = datmat[k + j * datmat_dim1] + con[k];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = 0.;
      i__3 = *n;
      for (j = 1; j <= i__3; ++j) {
        temp += w[j] * simi[j + i__ * simi_dim1];
      }
      if (k == mp) {
        temp = -temp;
      }
      a[i__ + k * a_dim1] = temp;
    }
  }

/* Calculate the values of sigma and eta, and set IFLAG=0 if the current */
/* simplex is not acceptable. */

  iflag = 1;
  parsig = alpha * rho;
  pareta = beta * rho;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    wsig = 0.;
    weta = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      d__1 = simi[j + i__ * simi_dim1];
      wsig += d__1 * d__1;
      d__1 = sim[i__ + j * sim_dim1];
      weta += d__1 * d__1;
    }
    vsig[j] = 1. / sqrt(wsig);
    veta[j] = sqrt(weta);
    if (vsig[j] < parsig || veta[j] > pareta) {
      iflag = 0;
    }
  }

/* If a new vertex is needed to improve acceptability, then decide which */
/* vertex to drop from the simplex. */

  if (ibrnch == 1 || iflag == 1) {
    goto L370;
  }
  jdrop = 0;
  temp = pareta;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    if (veta[j] > temp) {
      jdrop = j;
      temp = veta[j];
    }
  }
  if (jdrop == 0) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      if (vsig[j] < temp) {
        jdrop = j;
        temp = vsig[j];
      }
    }
  }

/* Calculate the step to the new vertex and its sign. */

  temp = gamma * rho * vsig[jdrop];
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dx[i__] = temp * simi[jdrop + i__ * simi_dim1];
  }
  cvmaxp = 0.;
  cvmaxm = 0.;
  i__1 = mp;
  for (k = 1; k <= i__1; ++k) {
    sum = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      sum += a[i__ + k * a_dim1] * dx[i__];
    }
    if (k < mp) {
      temp = datmat[k + np * datmat_dim1];
      d__1 = cvmaxp, d__2 = -sum - temp;
      cvmaxp = max(d__1,d__2);
      d__1 = cvmaxm, d__2 = sum - temp;
      cvmaxm = max(d__1,d__2);
    }
  }
  dxsign = 1.;
  if (parmu * (cvmaxp - cvmaxm) > sum + sum) {
    dxsign = -1.;
  }

/* Update the elements of SIM and SIMI, and set the next X. */

  temp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dx[i__] = dxsign * dx[i__];
    sim[i__ + jdrop * sim_dim1] = dx[i__];
    temp += simi[jdrop + i__ * simi_dim1] * dx[i__];
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    simi[jdrop + i__ * simi_dim1] /= temp;
  }
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    if (j != jdrop) {
      temp = 0.;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        temp += simi[j + i__ * simi_dim1] * dx[i__];
      }
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        simi[j + i__ * simi_dim1] -= temp * simi[jdrop + i__ * 
            simi_dim1];
      }
    }
    x[j] = sim[j + np * sim_dim1] + dx[j];
  }
  goto L40;

/* Calculate DX=x(*)-x(0). Branch if the length of DX is less than 0.5*RHO. */

L370:
  iz = 1;
  izdota = iz + *n * *n;
  ivmc = izdota + *n;
  isdirn = ivmc + mp;
  idxnew = isdirn + *n;
  ivmd = idxnew + *n;
  trstlp(n, m, &a[a_offset], &con[1], &rho, &dx[1], &ifull, &iact[1], &w[
      iz], &w[izdota], &w[ivmc], &w[isdirn], &w[idxnew], &w[ivmd]);
  if (ifull == 0) {
    temp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      d__1 = dx[i__];
      temp += d__1 * d__1;
    }
    if (temp < rho * .25 * rho) {
      ibrnch = 1;
      goto L550;
    }
  }

/* Predict the change to F and the new maximum constraint violation if the */
/* variables are altered from x(0) to x(0)+DX. */

  resnew = 0.;
  con[mp] = 0.;
  i__1 = mp;
  for (k = 1; k <= i__1; ++k) {
    sum = con[k];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      sum -= a[i__ + k * a_dim1] * dx[i__];
    }
    if (k < mp) {
      resnew = max(resnew,sum);
    }
  }

/* Increase PARMU if necessary and branch back if this change alters the */
/* optimal vertex. Otherwise PREREM and PREREC will be set to the predicted */
/* reductions in the merit function and the maximum constraint violation */
/* respectively. */

  barmu = 0.;
  prerec = datmat[*mpp + np * datmat_dim1] - resnew;
  if (prerec > 0.) {
    barmu = sum / prerec;
  }
  if (parmu < barmu * 1.5) {
    parmu = barmu * 2.;
    if (*iprint >= 2) {
      fprintf(stderr, "cobyla: increase in PARMU to %12.6E\n", parmu);
    }
    phi = datmat[mp + np * datmat_dim1] + parmu * datmat[*mpp + np * 
        datmat_dim1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      temp = datmat[mp + j * datmat_dim1] + parmu * datmat[*mpp + j * 
          datmat_dim1];
      if (temp < phi) {
        goto L140;
      }
      if (temp == phi && parmu == 0.f) {
        if (datmat[*mpp + j * datmat_dim1] < datmat[*mpp + np * 
            datmat_dim1]) {
          goto L140;
        }
      }
    }
  }
  prerem = parmu * prerec - sum;

/* Calculate the constraint and objective functions at x(*). Then find the */
/* actual reduction in the merit function. */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    x[i__] = sim[i__ + np * sim_dim1] + dx[i__];
  }
  ibrnch = 1;
  goto L40;
L440:
  vmold = datmat[mp + np * datmat_dim1] + parmu * datmat[*mpp + np * 
      datmat_dim1];
  vmnew = f + parmu * resmax;
  trured = vmold - vmnew;
  if (parmu == 0. && f == datmat[mp + np * datmat_dim1]) {
    prerem = prerec;
    trured = datmat[*mpp + np * datmat_dim1] - resmax;
  }

/* Begin the operations that decide whether x(*) should replace one of the */
/* vertices of the current simplex, the change being mandatory if TRURED is */
/* positive. Firstly, JDROP is set to the index of the vertex that is to be */
/* replaced. */

  ratio = 0.;
  if (trured <= 0.f) {
    ratio = 1.f;
  }
  jdrop = 0;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    temp = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      temp += simi[j + i__ * simi_dim1] * dx[i__];
    }
    temp = abs(temp);
    if (temp > ratio) {
      jdrop = j;
      ratio = temp;
    }
    sigbar[j] = temp * vsig[j];
  }

/* Calculate the value of ell. */

  edgmax = delta * rho;
  l = 0;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    if (sigbar[j] >= parsig || sigbar[j] >= vsig[j]) {
      temp = veta[j];
      if (trured > 0.) {
        temp = 0.;
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          d__1 = dx[i__] - sim[i__ + j * sim_dim1];
          temp += d__1 * d__1;
        }
        temp = sqrt(temp);
      }
      if (temp > edgmax) {
        l = j;
        edgmax = temp;
      }
    }
  }
  if (l > 0) {
    jdrop = l;
  }
  if (jdrop == 0) {
    goto L550;
  }

/* Revise the simplex by updating the elements of SIM, SIMI and DATMAT. */

  temp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sim[i__ + jdrop * sim_dim1] = dx[i__];
    temp += simi[jdrop + i__ * simi_dim1] * dx[i__];
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    simi[jdrop + i__ * simi_dim1] /= temp;
  }
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    if (j != jdrop) {
      temp = 0.;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        temp += simi[j + i__ * simi_dim1] * dx[i__];
      }
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        simi[j + i__ * simi_dim1] -= temp * simi[jdrop + i__ * 
            simi_dim1];
      }
    }
  }
  i__1 = *mpp;
  for (k = 1; k <= i__1; ++k) {
    datmat[k + jdrop * datmat_dim1] = con[k];
  }

/* Branch back for further iterations with the current RHO. */

  if (trured > 0. && trured >= prerem * .1) {
    goto L140;
  }
L550:
  if (iflag == 0) {
    ibrnch = 0;
    goto L140;
  }

/* Otherwise reduce RHO if it is not at its least value and reset PARMU. */

  if (rho > *rhoend) {
    rho *= .5;
    if (rho <= *rhoend * 1.5) {
      rho = *rhoend;
    }
    if (parmu > 0.) {
      denom = 0.;
      i__1 = mp;
      for (k = 1; k <= i__1; ++k) {
        cmin = datmat[k + np * datmat_dim1];
        cmax = cmin;
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          d__1 = cmin, d__2 = datmat[k + i__ * datmat_dim1];
          cmin = min(d__1,d__2);
          d__1 = cmax, d__2 = datmat[k + i__ * datmat_dim1];
          cmax = max(d__1,d__2);
        }
        if (k <= *m && cmin < cmax * .5) {
          temp = max(cmax,0.) - cmin;
          if (denom <= 0.) {
            denom = temp;
          } else {
            denom = min(denom,temp);
          }
        }
      }
      if (denom == 0.) {
        parmu = 0.;
      } else if (cmax - cmin < parmu * denom) {
        parmu = (cmax - cmin) / denom;
      }
    }
    if (*iprint >= 2) {
      fprintf(stderr, "cobyla: reduction in RHO to %12.6E and PARMU =%13.6E\n",
        rho, parmu);
    }
    if (*iprint == 2) {
      fprintf(stderr, "cobyla: NFVALS = %4d, F =%13.6E, MAXCV =%13.6E\n",
        nfvals, datmat[mp + np * datmat_dim1], datmat[*mpp + np * datmat_dim1]);

      fprintf(stderr, "cobyla: X =");
      i__1 = iptem;
      for (i__ = 1; i__ <= i__1; ++i__) {
        if (i__>1) fprintf(stderr, "  ");
        fprintf(stderr, "%13.6E", sim[i__ + np * sim_dim1]);
      }
      if (iptem < *n) {
        i__1 = *n;
        for (i__ = iptemp; i__ <= i__1; ++i__) {
          if (!((i__-1) % 4)) fprintf(stderr, "\ncobyla:  ");
          fprintf(stderr, "%15.6E", x[i__]);
        }
      }
      fprintf(stderr, "\n");
    }
    goto L140;
  }

/* Return the best calculated values of the variables. */

  if (*iprint >= 1) {
    fprintf(stderr, "cobyla: normal return.\n");
  }
  if (ifull == 1) {
    goto L620;
  }
L600:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    x[i__] = sim[i__ + np * sim_dim1];
  }
  f = datmat[mp + np * datmat_dim1];
  resmax = datmat[*mpp + np * datmat_dim1];
L620:
  if (*iprint >= 1) {
    fprintf(stderr, "cobyla: NFVALS = %4d, F =%13.6E, MAXCV =%13.6E\n",
      nfvals, f, resmax);
    i__1 = iptem;
    fprintf(stderr, "cobyla: X =");
    for (i__ = 1; i__ <= i__1; ++i__) {
      if (i__>1) fprintf(stderr, "  ");
      fprintf(stderr, "%13.6E", x[i__]);
    }
    if (iptem < *n) {
      i__1 = *n;
      for (i__ = iptemp; i__ <= i__1; ++i__) {
        if (!((i__-1) % 4)) fprintf(stderr, "\ncobyla:  ");
        fprintf(stderr, "%15.6E", x[i__]);
      }
    }
    fprintf(stderr, "\n");
  }
  *maxfun = nfvals;
  return rc;
} /* cobylb */

/* ------------------------------------------------------------------------- */
int trstlp(int *n, int *m, double *a, 
    double *b, double *rho, double *dx, int *ifull, 
    int *iact, double *z__, double *zdota, double *vmultc,
     double *sdirn, double *dxnew, double *vmultd)
{
  /* System generated locals */
  int a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
  double d__1, d__2;

  /* Local variables */
  double alpha, tempa;
  double beta;
  double optnew, stpful, sum, tot, acca, accb;
  double ratio, vsave, zdotv, zdotw, dd;
  double sd;
  double sp, ss, resold = 0.0, zdvabs, zdwabs, sumabs, resmax, optold;
  double spabs;
  double temp, step;
  int icount;
  int iout, i__, j, k;
  int isave;
  int kk;
  int kl, kp, kw;
  int nact, icon = 0, mcon;
  int nactx = 0;


/* This subroutine calculates an N-component vector DX by applying the */
/* following two stages. In the first stage, DX is set to the shortest */
/* vector that minimizes the greatest violation of the constraints */
/*   A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K=2,3,...,M, */
/* subject to the Euclidean length of DX being at most RHO. If its length is */
/* strictly less than RHO, then we use the resultant freedom in DX to */
/* minimize the objective function */
/*      -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N) */
/* subject to no increase in any greatest constraint violation. This */
/* notation allows the gradient of the objective function to be regarded as */
/* the gradient of a constraint. Therefore the two stages are distinguished */
/* by MCON .EQ. M and MCON .GT. M respectively. It is possible that a */
/* degeneracy may prevent DX from attaining the target length RHO. Then the */
/* value IFULL=0 would be set, but usually IFULL=1 on return. */

/* In general NACT is the number of constraints in the active set and */
/* IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT */
/* contains a permutation of the remaining constraint indices. Further, Z is */
/* an orthogonal matrix whose first NACT columns can be regarded as the */
/* result of Gram-Schmidt applied to the active constraint gradients. For */
/* J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th */
/* column of Z with the gradient of the J-th active constraint. DX is the */
/* current vector of variables and here the residuals of the active */
/* constraints should be zero. Further, the active constraints have */
/* nonnegative Lagrange multipliers that are held at the beginning of */
/* VMULTC. The remainder of this vector holds the residuals of the inactive */
/* constraints at DX, the ordering of the components of VMULTC being in */
/* agreement with the permutation of the indices of the constraints that is */
/* in IACT. All these residuals are nonnegative, which is achieved by the */
/* shift RESMAX that makes the least residual zero. */

/* Initialize Z and some other variables. The value of RESMAX will be */
/* appropriate to DX=0, while ICON will be the index of a most violated */
/* constraint if RESMAX is positive. Usually during the first stage the */
/* vector SDIRN gives a search direction that reduces all the active */
/* constraint violations by one simultaneously. */

  /* Parameter adjustments */
  z_dim1 = *n;
  z_offset = 1 + z_dim1 * 1;
  z__ -= z_offset;
  a_dim1 = *n;
  a_offset = 1 + a_dim1 * 1;
  a -= a_offset;
  --b;
  --dx;
  --iact;
  --zdota;
  --vmultc;
  --sdirn;
  --dxnew;
  --vmultd;

  /* Function Body */
  *ifull = 1;
  mcon = *m;
  nact = 0;
  resmax = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
      z__[i__ + j * z_dim1] = 0.;
    }
    z__[i__ + i__ * z_dim1] = 1.;
    dx[i__] = 0.;
  }
  if (*m >= 1) {
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
      if (b[k] > resmax) {
        resmax = b[k];
        icon = k;
      }
    }
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
      iact[k] = k;
      vmultc[k] = resmax - b[k];
    }
  }
  if (resmax == 0.) {
    goto L480;
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sdirn[i__] = 0.;
  }

/* End the current stage of the calculation if 3 consecutive iterations */
/* have either failed to reduce the best calculated value of the objective */
/* function or to increase the number of active constraints since the best */
/* value was calculated. This strategy prevents cycling, but there is a */
/* remote possibility that it will cause premature termination. */

L60:
  optold = 0.;
  icount = 0;
L70:
  if (mcon == *m) {
    optnew = resmax;
  } else {
    optnew = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      optnew -= dx[i__] * a[i__ + mcon * a_dim1];
    }
  }
  if (icount == 0 || optnew < optold) {
    optold = optnew;
    nactx = nact;
    icount = 3;
  } else if (nact > nactx) {
    nactx = nact;
    icount = 3;
  } else {
    --icount;
    if (icount == 0) {
      goto L490;
    }
  }

/* If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to */
/* the active set. Apply Givens rotations so that the last N-NACT-1 columns */
/* of Z are orthogonal to the gradient of the new constraint, a scalar */
/* product being set to zero if its nonzero value could be due to computer */
/* rounding errors. The array DXNEW is used for working space. */

  if (icon <= nact) {
    goto L260;
  }
  kk = iact[icon];
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dxnew[i__] = a[i__ + kk * a_dim1];
  }
  tot = 0.;
  k = *n;
L100:
  if (k > nact) {
    sp = 0.;
    spabs = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = z__[i__ + k * z_dim1] * dxnew[i__];
      sp += temp;
      spabs += abs(temp);
    }
    acca = spabs + abs(sp) * .1;
    accb = spabs + abs(sp) * .2;
    if (spabs >= acca || acca >= accb) {
      sp = 0.;
    }
    if (tot == 0.) {
      tot = sp;
    } else {
      kp = k + 1;
      temp = sqrt(sp * sp + tot * tot);
      alpha = sp / temp;
      beta = tot / temp;
      tot = temp;
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        temp = alpha * z__[i__ + k * z_dim1] + beta * z__[i__ + kp * 
            z_dim1];
        z__[i__ + kp * z_dim1] = alpha * z__[i__ + kp * z_dim1] - 
            beta * z__[i__ + k * z_dim1];
        z__[i__ + k * z_dim1] = temp;
      }
    }
    --k;
    goto L100;
  }

/* Add the new constraint if this can be done without a deletion from the */
/* active set. */

  if (tot != 0.) {
    ++nact;
    zdota[nact] = tot;
    vmultc[icon] = vmultc[nact];
    vmultc[nact] = 0.;
    goto L210;
  }

/* The next instruction is reached if a deletion has to be made from the */
/* active set in order to make room for the new active constraint, because */
/* the new constraint gradient is a linear combination of the gradients of */
/* the old active constraints. Set the elements of VMULTD to the multipliers */
/* of the linear combination. Further, set IOUT to the index of the */
/* constraint to be deleted, but branch if no suitable index can be found. */

  ratio = -1.;
  k = nact;
L130:
  zdotv = 0.;
  zdvabs = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    temp = z__[i__ + k * z_dim1] * dxnew[i__];
    zdotv += temp;
    zdvabs += abs(temp);
  }
  acca = zdvabs + abs(zdotv) * .1;
  accb = zdvabs + abs(zdotv) * .2;
  if (zdvabs < acca && acca < accb) {
    temp = zdotv / zdota[k];
    if (temp > 0. && iact[k] <= *m) {
      tempa = vmultc[k] / temp;
      if (ratio < 0. || tempa < ratio) {
        ratio = tempa;
        iout = k;
      }
    }
    if (k >= 2) {
      kw = iact[k];
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        dxnew[i__] -= temp * a[i__ + kw * a_dim1];
      }
    }
    vmultd[k] = temp;
  } else {
    vmultd[k] = 0.;
  }
  --k;
  if (k > 0) {
    goto L130;
  }
  if (ratio < 0.) {
    goto L490;
  }

/* Revise the Lagrange multipliers and reorder the active constraints so */
/* that the one to be replaced is at the end of the list. Also calculate the */
/* new value of ZDOTA(NACT) and branch if it is not acceptable. */

  i__1 = nact;
  for (k = 1; k <= i__1; ++k) {
    d__1 = 0., d__2 = vmultc[k] - ratio * vmultd[k];
    vmultc[k] = max(d__1,d__2);
  }
  if (icon < nact) {
    isave = iact[icon];
    vsave = vmultc[icon];
    k = icon;
L170:
    kp = k + 1;
    kw = iact[kp];
    sp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sp += z__[i__ + k * z_dim1] * a[i__ + kw * a_dim1];
    }
    d__1 = zdota[kp];
    temp = sqrt(sp * sp + d__1 * d__1);
    alpha = zdota[kp] / temp;
    beta = sp / temp;
    zdota[kp] = alpha * zdota[k];
    zdota[k] = temp;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = alpha * z__[i__ + kp * z_dim1] + beta * z__[i__ + k * 
          z_dim1];
      z__[i__ + kp * z_dim1] = alpha * z__[i__ + k * z_dim1] - beta * 
          z__[i__ + kp * z_dim1];
      z__[i__ + k * z_dim1] = temp;
    }
    iact[k] = kw;
    vmultc[k] = vmultc[kp];
    k = kp;
    if (k < nact) {
      goto L170;
    }
    iact[k] = isave;
    vmultc[k] = vsave;
  }
  temp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    temp += z__[i__ + nact * z_dim1] * a[i__ + kk * a_dim1];
  }
  if (temp == 0.) {
    goto L490;
  }
  zdota[nact] = temp;
  vmultc[icon] = 0.;
  vmultc[nact] = ratio;

/* Update IACT and ensure that the objective function continues to be */
/* treated as the last active constraint when MCON>M. */

L210:
  iact[icon] = iact[nact];
  iact[nact] = kk;
  if (mcon > *m && kk != mcon) {
    k = nact - 1;
    sp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sp += z__[i__ + k * z_dim1] * a[i__ + kk * a_dim1];
    }
    d__1 = zdota[nact];
    temp = sqrt(sp * sp + d__1 * d__1);
    alpha = zdota[nact] / temp;
    beta = sp / temp;
    zdota[nact] = alpha * zdota[k];
    zdota[k] = temp;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = alpha * z__[i__ + nact * z_dim1] + beta * z__[i__ + k * 
          z_dim1];
      z__[i__ + nact * z_dim1] = alpha * z__[i__ + k * z_dim1] - beta * 
          z__[i__ + nact * z_dim1];
      z__[i__ + k * z_dim1] = temp;
    }
    iact[nact] = iact[k];
    iact[k] = kk;
    temp = vmultc[k];
    vmultc[k] = vmultc[nact];
    vmultc[nact] = temp;
  }

/* If stage one is in progress, then set SDIRN to the direction of the next */
/* change to the current vector of variables. */

  if (mcon > *m) {
    goto L320;
  }
  kk = iact[nact];
  temp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    temp += sdirn[i__] * a[i__ + kk * a_dim1];
  }
  temp += -1.;
  temp /= zdota[nact];
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sdirn[i__] -= temp * z__[i__ + nact * z_dim1];
  }
  goto L340;

/* Delete the constraint that has the index IACT(ICON) from the active set. */

L260:
  if (icon < nact) {
    isave = iact[icon];
    vsave = vmultc[icon];
    k = icon;
L270:
    kp = k + 1;
    kk = iact[kp];
    sp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sp += z__[i__ + k * z_dim1] * a[i__ + kk * a_dim1];
    }
    d__1 = zdota[kp];
    temp = sqrt(sp * sp + d__1 * d__1);
    alpha = zdota[kp] / temp;
    beta = sp / temp;
    zdota[kp] = alpha * zdota[k];
    zdota[k] = temp;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = alpha * z__[i__ + kp * z_dim1] + beta * z__[i__ + k * 
          z_dim1];
      z__[i__ + kp * z_dim1] = alpha * z__[i__ + k * z_dim1] - beta * 
          z__[i__ + kp * z_dim1];
      z__[i__ + k * z_dim1] = temp;
    }
    iact[k] = kk;
    vmultc[k] = vmultc[kp];
    k = kp;
    if (k < nact) {
      goto L270;
    }
    iact[k] = isave;
    vmultc[k] = vsave;
  }
  --nact;

/* If stage one is in progress, then set SDIRN to the direction of the next */
/* change to the current vector of variables. */

  if (mcon > *m) {
    goto L320;
  }
  temp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    temp += sdirn[i__] * z__[i__ + (nact + 1) * z_dim1];
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sdirn[i__] -= temp * z__[i__ + (nact + 1) * z_dim1];
  }
  goto L340;

/* Pick the next search direction of stage two. */

L320:
  temp = 1. / zdota[nact];
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sdirn[i__] = temp * z__[i__ + nact * z_dim1];
  }

/* Calculate the step to the boundary of the trust region or take the step */
/* that reduces RESMAX to zero. The two statements below that include the */
/* factor 1.0E-6 prevent some harmless underflows that occurred in a test */
/* calculation. Further, we skip the step if it could be zero within a */
/* reasonable tolerance for computer rounding errors. */

L340:
  dd = *rho * *rho;
  sd = 0.;
  ss = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    if ((d__1 = dx[i__], abs(d__1)) >= *rho * 1e-6f) {
      d__2 = dx[i__];
      dd -= d__2 * d__2;
    }
    sd += dx[i__] * sdirn[i__];
    d__1 = sdirn[i__];
    ss += d__1 * d__1;
  }
  if (dd <= 0.) {
    goto L490;
  }
  temp = sqrt(ss * dd);
  if (abs(sd) >= temp * 1e-6f) {
    temp = sqrt(ss * dd + sd * sd);
  }
  stpful = dd / (temp + sd);
  step = stpful;
  if (mcon == *m) {
    acca = step + resmax * .1;
    accb = step + resmax * .2;
    if (step >= acca || acca >= accb) {
      goto L480;
    }
    step = min(step,resmax);
  }

/* Set DXNEW to the new variables if STEP is the steplength, and reduce */
/* RESMAX to the corresponding maximum residual if stage one is being done. */
/* Because DXNEW will be changed during the calculation of some Lagrange */
/* multipliers, it will be restored to the following value later. */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dxnew[i__] = dx[i__] + step * sdirn[i__];
  }
  if (mcon == *m) {
    resold = resmax;
    resmax = 0.;
    i__1 = nact;
    for (k = 1; k <= i__1; ++k) {
      kk = iact[k];
      temp = b[kk];
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        temp -= a[i__ + kk * a_dim1] * dxnew[i__];
      }
      resmax = max(resmax,temp);
    }
  }

/* Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A */
/* device is included to force VMULTD(K)=0.0 if deviations from this value */
/* can be attributed to computer rounding errors. First calculate the new */
/* Lagrange multipliers. */

  k = nact;
L390:
  zdotw = 0.;
  zdwabs = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    temp = z__[i__ + k * z_dim1] * dxnew[i__];
    zdotw += temp;
    zdwabs += abs(temp);
  }
  acca = zdwabs + abs(zdotw) * .1;
  accb = zdwabs + abs(zdotw) * .2;
  if (zdwabs >= acca || acca >= accb) {
    zdotw = 0.;
  }
  vmultd[k] = zdotw / zdota[k];
  if (k >= 2) {
    kk = iact[k];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      dxnew[i__] -= vmultd[k] * a[i__ + kk * a_dim1];
    }
    --k;
    goto L390;
  }
  if (mcon > *m) {
    d__1 = 0., d__2 = vmultd[nact];
    vmultd[nact] = max(d__1,d__2);
  }

/* Complete VMULTC by finding the new constraint residuals. */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dxnew[i__] = dx[i__] + step * sdirn[i__];
  }
  if (mcon > nact) {
    kl = nact + 1;
    i__1 = mcon;
    for (k = kl; k <= i__1; ++k) {
      kk = iact[k];
      sum = resmax - b[kk];
      sumabs = resmax + (d__1 = b[kk], abs(d__1));
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        temp = a[i__ + kk * a_dim1] * dxnew[i__];
        sum += temp;
        sumabs += abs(temp);
      }
      acca = sumabs + abs(sum) * .1f;
      accb = sumabs + abs(sum) * .2f;
      if (sumabs >= acca || acca >= accb) {
        sum = 0.f;
      }
      vmultd[k] = sum;
    }
  }

/* Calculate the fraction of the step from DX to DXNEW that will be taken. */

  ratio = 1.;
  icon = 0;
  i__1 = mcon;
  for (k = 1; k <= i__1; ++k) {
    if (vmultd[k] < 0.) {
      temp = vmultc[k] / (vmultc[k] - vmultd[k]);
      if (temp < ratio) {
        ratio = temp;
        icon = k;
      }
    }
  }

/* Update DX, VMULTC and RESMAX. */

  temp = 1. - ratio;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dx[i__] = temp * dx[i__] + ratio * dxnew[i__];
  }
  i__1 = mcon;
  for (k = 1; k <= i__1; ++k) {
    d__1 = 0., d__2 = temp * vmultc[k] + ratio * vmultd[k];
    vmultc[k] = max(d__1,d__2);
  }
  if (mcon == *m) {
    resmax = resold + ratio * (resmax - resold);
  }

/* If the full step is not acceptable then begin another iteration. */
/* Otherwise switch to stage two or end the calculation. */

  if (icon > 0) {
    goto L70;
  }
  if (step == stpful) {
    goto L500;
  }
L480:
  mcon = *m + 1;
  icon = mcon;
  iact[mcon] = mcon;
  vmultc[mcon] = 0.;
  goto L60;

/* We employ any freedom that may be available to reduce the objective */
/* function before returning a DX whose length is less than RHO. */

L490:
  if (mcon == *m) {
    goto L480;
  }
  *ifull = 0;
L500:
  return 0;
} /* trstlp */












using namespace std;


FILE *stream;  
FILE *file,*outfile;  	  
 
double meanlikehood[500];    
double Sel_Int[11] = {2.326,2.603,1.755,1.554,1.4,1.271,1.159,1.058,0.966,0.88,0.798};  
double SelInt; 

int cantexp;  
int now;  
int vars;  
int auxMax;  
double Max;  
double  Trunc;  
int psize;  
int  Tour;  
int func;  
int ExperimentMode;  
int Ntrees;  
int Elit;  
int succexp;  
double meangen;   
int Nsteps;  
int InitTreeStructure;  
int VisibleChoiceVar;  
int Maxgen;  
int printvals;   
unsigned int Card;  
int seed;  
int* params;
int fun;  
int *timevector; 
char filedetails[30]; 
char MatrixFileName[30]; 
int BestElitism; 
double MaxMixtProb; 
double S_alpha;  
int StopCrit; //Stop criteria to stop the MT learning alg. 
int Prior; 
double Complex; 
int Coeftype;  
unsigned *Cardinalities;  
int Mutation; 
int CliqMaxLength; 
int MaxNumCliq; 
int OldWaySel; 
int LearningType;
int TypeMixture; 
int Cycles; 
 
double meaneval;  
double BestEval,AbsBestEval,AuxBest; 
int TruncMax; 
int NPoints;  
double  *BestInd, *AbsBestInd;  
ContPopul *pop,*selpop,*elitpop,*compact_pop; 
double *fvect; 
int  nsucc;

int Clock;

div_t ImproveStop;
 double auxtime, alltime,bestalltime;
 time_t ltime_init,ltime_end;
 struct tm *gmt;
 struct tm *gmtnew;

int LEARNEBNA=1;  
int EBNASCORE=K2_SCORE;
double  EBNA_ALPHA =0.05;
int  EBNA_SIMUL = PLS;
 

int TotEvaluations;
int EvaluationMode;
int currentexp;
int length;
long int explength;
 int MaxMPC;
 int TypeMPC;


// Global variables trajectory problem
const double pi = acos(-1.0);
double minbound[22];
double maxbound[22];
double currentminbound[22];
double currentmaxbound[22];

long int locusvalues[22];
int dim;
vector<int> sequence(22);
vector<double> t(22);
double *Delta_V;
	dsm_customobject c_o; // empty in this case
	double rp;
	double e;	
	double Isp;
	double mass;	
	double AUdist;
	double DVtotal;
	double DVonboard;
 int tract_vars;




void init_time()
{
 time( &ltime_init );
 gmt = localtime( &ltime_init );
 auxtime = - ( gmt->tm_mday * 86400 + gmt->tm_hour*3600 + gmt->tm_min*60 + gmt->tm_sec);
}


void end_time()
{
  time( &ltime_end );
  gmtnew = localtime( &ltime_end );
  auxtime = auxtime + gmtnew->tm_mday * 86400 + gmtnew->tm_hour*3600+gmtnew->tm_min*60+gmtnew->tm_sec;
}


 
int Selection() 
{ 
   int NPoints=0; 
 
   if (Tour==0)  
         {  
           pop->TruncSel(selpop,TruncMax); 
           selpop->UniformProb(TruncMax,fvect);
           //selpop->BotzmannDist(1.0,fvect);          
           NPoints = selpop->CompactPopNew(compact_pop,fvect); 
	   //NPoints = TruncMax;
           //compact_pop->CopyPop(selpop);           
	 } 
     else if(Tour==1) //Tournament selection 
	 {  
	   pop->TournSel(selpop,TruncMax); 
           selpop->UniformProb(psize,fvect); 
           NPoints = selpop->CompactPopNew(compact_pop,fvect); 
	 }  
    else if(Tour==2) //Proportional selection 
	 {  
	   pop->ProporDist(fvect);   
	   if (OldWaySel) 
           { 
            selpop->SUSSel(psize,pop,fvect);  
            selpop->UniformProb(psize,fvect);    
            NPoints = selpop->CompactPopNew(compact_pop,fvect);         
           } 
           else NPoints = pop->CompactPopNew(compact_pop,fvect);                         
          }  
     else if(Tour==3) //Boltzman selection 
	 {  
	   pop->BotzmannDist(1.0,fvect); 
	   if (OldWaySel) 
           { 
            selpop->SUSSel(psize,pop,fvect);  
            selpop->UniformProb(psize,fvect);    
            NPoints = selpop->CompactPopNew(compact_pop,fvect); 
           } 
           else NPoints = pop->CompactPopNew(compact_pop,fvect); 
           
	 }  
  if (Tour>0 || (Tour==0 && Elit>TruncMax)) pop->TruncSel(elitpop,Elit);  
    
   return NPoints; 
} 
 
inline void FindBestVal() 
{     
      if(Elit && Tour != 0)  
          { 
            BestEval =elitpop->Evaluations[0]; 
            BestInd = elitpop->S[0]; 
	  } 
      else if(Tour==0) 
      {  
        BestEval = selpop->Evaluations[0]; 
        BestInd = selpop->S[0]; 
      } 
      else  
          { 
	   int auxind =  pop->FindBestIndPos();  
           BestInd =  pop->S[auxind]; 
           BestEval = pop->Evaluations[auxind]; 
          } 
} 
 


inline void InitContPopulations() 
{
  int i; 
 if (Tour==0) 
   { 
     TruncMax = int(psize*Trunc);  
   
     if (BestElitism)  Elit = TruncMax;   //Only for Trunc Selection  
     selpop = new ContPopul(TruncMax,vars,Elit,minbound,maxbound);  
   }  
  else selpop = new ContPopul(psize,vars,Elit,minbound,maxbound);  
 
  if (Tour>0 || (Tour==0 && Elit>TruncMax)) 
  elitpop = new ContPopul(Elit,vars,Elit,minbound,maxbound); 
  pop = new ContPopul(psize,vars,Elit,minbound,maxbound);
  pop->RandInit();  
  compact_pop = new ContPopul(psize,vars,Elit,minbound,maxbound);  
  fvect = new double[psize];
 
 } 



 
inline void DeleteContPopulations() 
{ 
  delete compact_pop; 
  delete pop;  
  delete selpop;  
  if (Tour>0 || (Tour==0 && Elit>TruncMax)) delete elitpop; 
 delete[] fvect; 
} 


void ConvertFromCont(double* sol)
{
  int j;
  //for (j=0;j<tract_vars;j++) cout<<sol[j]<<" ";
  //cout<<endl;
  for (j=0;j<tract_vars;j++)
     {
       if(sol[j]>=currentminbound[j] && sol[j]<=currentmaxbound[j])
	 t[j] = sol[j];
       else 
         t[j] =  currentminbound[j] + ((currentmaxbound[j]-currentminbound[j]) *myrand());
       sol[j] = t[j];
     }
  //for (j=0;j<tract_vars;j++) cout<<t[j]<<" ";
  //cout<<endl;

//cout<<endl;
}


void ConvertFromContA(double* sol)
{
  int j;
  for (j=0;j<tract_vars;j++) t[j]=sol[j];
  
  //for (j=0;j<tract_vars;j++) cout<<t[j]<<" ";
  //cout<<endl;

  //cout<<endl;
}



void InitBounds()
{
  int i;
  for(i=0;i<tract_vars;i++)
    {
      currentminbound[i] =  minbound[i]; 
      currentmaxbound[i] =  maxbound[i]; 

    }
}


void UpdateBoundsFromSol(double* sol)
{
  int i;
  double ntmin, ntmax;
 
  ConvertFromCont(sol);

/*
  for(i=0;i<tract_vars;i++)
    {
      ntmin = (t[i] - (t[i]-minbound[i]) /2);
      ntmax = (t[i] + (maxbound[i]-t[i]) /2);
      if ( (ntmin < minbound[i]) || (ntmax==t[i]) || (ntmin==t[i]) ) 
      currentminbound[i] = minbound[i];
      else  currentminbound[i] = ntmin;
      if ( (ntmax > maxbound[i]) || (ntmax==t[i]) || (ntmin==t[i]))   currentmaxbound[i] = maxbound[i];
      else  currentmaxbound[i] = ntmax;
    }
*/   

  for(i=0;i<tract_vars;i++)
    {
      ntmin = (t[i] - (t[i]-currentminbound[i]) /2);
      ntmax = (t[i] + (currentmaxbound[i]-t[i]) /2);
      if ( ! ((ntmin < minbound[i]) || (ntmax==t[i]) || (ntmin==t[i])) ) 
      currentminbound[i] = ntmin;
      if (!( (ntmax > maxbound[i]) || (ntmax==t[i]) || (ntmin==t[i])))
      currentmaxbound[i] = ntmax;
    }
 

   cout<<"CURRENTBOUNDS: ";
    for(i=0;i<tract_vars;i++) cout<<currentminbound[i]<<"-"<<currentmaxbound[i]<<", ";
    cout<<endl;
    

}

 
void InitTrajectory()
{

 
 
 #if MGADSM_PROBLEM_TYPE == time2AUs  // SAGAS
	dim = 3;
	sequence[0] = 3;
	sequence[1] = 3;
	sequence[2] = 5;
	
	Delta_V = new double[dim+1];
	
      t[0] = 7020.49;
	t[1] = 5.34817;
	t[2] = 1;
	t[3] = 0.498915;
	t[4] = 788.763;
	t[5] = 484.349;
	t[6] = 0.4873;
	t[7] = 0.01;
	t[8] = 1.05;
	t[9] = 10.8516;
	t[10] = -1.57191;
	t[11] = -0.685429;
	
	//dsm_customobject c_o; // empty in this case
	rp = 3950;
	e = 0.98;	
	Isp = 0.0;
	mass = 0.0;	
	AUdist = 50.0;
	DVtotal = 6.782;
	DVonboard = 1.782;

     
        tract_vars = 12;

        minbound[0] = 7000.0;   maxbound[0] = 9100.0;
        minbound[1] = 0.0;      maxbound[1] = 7.0;
        minbound[2] = 0.0;      maxbound[2] = 1.0;
        minbound[3] = 0.0;      maxbound[3] = 1.0;
        minbound[4] = 50.0;     maxbound[4] = 2000.0;
        minbound[5] = 300.0;    maxbound[5] = 2000.0;
        minbound[6] = 0.01;     maxbound[6] = 0.9;
        minbound[7] = 0.01;     maxbound[7] = 0.9; 
        minbound[8] = 1.05;     maxbound[8] = 7.0;
        minbound[9] = 8.0;      maxbound[9] = 500.0;
        minbound[10] = -1*pi;   maxbound[10] = pi;
        minbound[11] = -1*pi;   maxbound[11] = pi;
      

#elif MGADSM_PROBLEM_TYPE == total_DV_rndv  // Cassini with DSM

#if total_DV_rndv_problem == cassini
	dim = 6;
	sequence[0] = 3;
	sequence[1] = 2;
	sequence[2] = 2;
	sequence[3] = 3;
	sequence[4] = 5;
	sequence[5] = 6;
	
	Delta_V = new double[dim+1];
	       
        t[0] = -815.144;
	t[1] = 3;
	t[2] = 0.623166;
	t[3] = 0.444834;
	t[4] = 197.334;
	t[5] = 425.171;
	t[6] = 56.8856;
	t[7] = 578.523;
	t[8] = 2067.98;
	t[9] = 0.01;
	t[10] = 0.470415;
	t[11] = 0.01;
	t[12] = 0.0892135;
	t[13] = 0.9;
	t[14] = 1.05044;
	t[15] = 1.38089;
	t[16] = 1.18824;
	t[17] = 76.5066;
	t[18] = -1.57225;
	t[19] = -2.01799;
	t[20] = -1.52153;
	t[21] = -1.5169;

	// all the rest is empty in this case
	//dsm_customobject c_o; 
	rp = 0;
	e = 0;	
	Isp = 0.0;
	mass = 0.0;
	AUdist = 0;
	DVtotal = 0;
	DVonboard = 0;


  tract_vars = 22;

        minbound[0] = -1000.0;   maxbound[0] = 0.0;
        minbound[1] = 3.0;      maxbound[1] = 5.0;
        minbound[2] = 0.0;      maxbound[2] = 1.0;
        minbound[3] = 0.0;      maxbound[3] = 1.0;
        minbound[4] = 100.0;     maxbound[4] = 400.0;
        minbound[5] = 100.0;    maxbound[5] =  500.0;
        minbound[6] = 30.0;     maxbound[6] = 300.0;
        minbound[7] = 400.0;     maxbound[7] = 1600.0; 
        minbound[8] = 800.0;     maxbound[8] = 2200.0;
        minbound[9] = 0.01;     maxbound[9] = 0.9;
        minbound[10] = 0.01;     maxbound[10] = 0.9; 
        minbound[11] = 0.01;     maxbound[11] = 0.9;
        minbound[12] = 0.01;     maxbound[12] = 0.9; 
        minbound[13] = 0.01;     maxbound[13] = 0.9;   
        minbound[14] = 1.05;     maxbound[14] = 6.0; 
        minbound[15] = 1.05;     maxbound[15] = 6.0; 
        minbound[16] = 1.15;     maxbound[16] = 6.5; 
        minbound[17] = 1.7;      maxbound[17] = 291.0;
        minbound[18] = -1*pi;   maxbound[18] = pi;
        minbound[19] = -1*pi;   maxbound[19] = pi;
        minbound[20] = -1*pi;   maxbound[20] = pi;
        minbound[21] = -1*pi;   maxbound[21] = pi;
#endif
#if total_DV_rndv_problem == messenger

	dim = 5;
	sequence[0] = 3;
	sequence[1] = 3;
	sequence[2] = 2;
	sequence[3] = 2;
	sequence[4] = 1;

	Delta_V = new double[dim+1];
	
        t[0] = 2363.36;
	t[1] = 1.68003;
	t[2] = 0.381885;
	t[3] = 0.512516;
	t[4] = 400;
	t[5] = 173.848;
	t[6] = 224.702;
	t[7] = 211.803;
	t[8] = 0.238464;
	t[9] = 0.265663;
	t[10] = 0.149817;
	t[11] = 0.485908;
	t[12] = 1.34411;
	t[13] = 3.49751;
	t[14] = 1.1;
	t[15] = 1.29892;
	t[16] = 2.49324;
	t[17] = 1.81426;

	// all the rest is empty in this case
	//dsm_customobject c_o; 
	rp = 0;
	e = 0;	
	Isp = 0.0;
	mass = 0.0;
	AUdist = 0;
	DVtotal = 0;
	DVonboard = 0;


        tract_vars = 18;

        minbound[0] = 1000.0;   maxbound[0] = 4100.0;
        minbound[1] = 1.0;      maxbound[1] = 5.0;
        minbound[2] = 0.0;      maxbound[2] = 1.0;
        minbound[3] = 0.0;      maxbound[3] = 1.0;
        minbound[4] = 200.0;    maxbound[4] = 400.0;
        minbound[5] = 30.0;    maxbound[5] =  400.0;
        minbound[6] = 30.0;    maxbound[6] =  400.0;
        minbound[7] = 30.0;    maxbound[7] =  400.0;
        minbound[8] = 0.01;     maxbound[8] = 0.99;
        minbound[9] = 0.01;     maxbound[9] = 0.99;
        minbound[10] = 0.01;     maxbound[10] = 0.99;
        minbound[11] = 0.01;     maxbound[11] = 0.99;
        minbound[12] = 1.1;      maxbound[12] = 6.0;
        minbound[13] = 1.1;      maxbound[13] = 6.0;
        minbound[14] = 1.1;      maxbound[14] = 6.0;
        minbound[15] = -1*pi;   maxbound[15] = pi;
        minbound[16] = -1*pi;   maxbound[16] = pi;
        minbound[17] = -1*pi;   maxbound[17] = pi;
    

#endif

#elif MGADSM_PROBLEM_TYPE == rndv

	dim = 6;
	
	sequence[0] = 3;
	sequence[1] = 3;
	sequence[2] = 4;
	sequence[3] = 3;
	sequence[4] = 3;
	sequence[5] = 10;

	Delta_V = new double[dim+1];
	
        t[0] = 1524.25;
	t[1] = 3.95107;
	t[2] = 0.738307;
	t[3] = 0.298318;
	t[4] = 365.123;
	t[5] = 728.902;
	t[6] = 256.049;
	t[7] = 730.485;
	t[8] = 1850;
	t[9] = 0.199885;
	t[10] = 0.883382;
	t[11] = 0.194587;
	t[12] = 0.0645205;
	t[13] = 0.493077;
	t[14] = 1.05;
	t[15] = 1.05;
	t[16] = 1.05;
	t[17] = 1.36925;
	t[18] = -1.74441;
	t[19] = 1.85201;
	t[20] = -2.61644;
	t[21] = -1.53468;
	
	//dsm_customobject c_o; 
	c_o.keplerian[0] = 3.50294972836275;
	c_o.keplerian[1] = 0.6319356;
	c_o.keplerian[2] =  7.12723;
	c_o.keplerian[3] = 	50.92302;
	c_o.keplerian[4] =  11.36788;
	c_o.keplerian[5] = 0.0;
	c_o.epoch = 52504.23754000012;
	c_o.mu = 0.0;

	// all the rest is empty in this case
	rp = 0;
	e = 0;	
	Isp = 0.0;
	mass = 0.0;
	AUdist = 0;
	DVtotal = 0;
	DVonboard = 0;

        tract_vars = 22;

        minbound[0] = 1460.0;   maxbound[0] = 1825.0;
        minbound[1] = 3.0;      maxbound[1] = 5.0;
        minbound[2] = 0.0;      maxbound[2] = 1.0;
        minbound[3] = 0.0;      maxbound[3] = 1.0;
        minbound[4] = 300.0;    maxbound[4] = 500.0;
        minbound[5] = 150.0;    maxbound[5] = 800.0;
        minbound[6] = 150.0;    maxbound[6] =  800.0;
        minbound[7] = 300.0;    maxbound[7] =  800.0;
        minbound[8] = 700.0;     maxbound[8] = 1850.0;
        minbound[9] = 0.01;     maxbound[9] = 0.9;
        minbound[10] = 0.01;     maxbound[10] = 0.9;
        minbound[11] = 0.01;     maxbound[11] = 0.9;
        minbound[12] = 0.01;      maxbound[12] = 0.9;
        minbound[13] = 0.01;      maxbound[13] = 0.9;
        minbound[14] = 1.05;      maxbound[14] = 9.0;
        minbound[15] = 1.05;      maxbound[15] = 9.0;
        minbound[16] = 1.05;      maxbound[16] = 9.0;
        minbound[17] = 1.05;      maxbound[17] = 9.0;
        minbound[18] = -1*pi;   maxbound[18] = pi;
        minbound[19] = -1*pi;   maxbound[19] = pi;
        minbound[20] = -1*pi;   maxbound[20] = pi;
	minbound[21] = -1*pi;   maxbound[21] = pi;

    #endif

	InitBounds();	
}







int ContRandom(int nvars, double* newsol, int trials, double* TBest, double* TheBestsol)
{
  double auxsol[500];  //  500 is maximum number of variables
  int eval_local,i,k,l,a,b,c, door,improve,norigvars;
  double obj_f,stepk,stepl,J,BestJ,oldval_k,oldval_l;
  double TheBest;

   ConvertFromCont(newsol);
   MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);	
   if (isnan(obj_f)) obj_f = 1000000.0; 
   eval_local = 1;      

   for(i=0;i<nvars;i++) TheBestsol[i] = newsol[i];
   TheBest = obj_f;
   BestJ = obj_f;

 improve = 1;
 a = 0;
 norigvars = nvars;

// while(a < trials) {
   b = 1;
   c = 0;
 while(a<trials && b<1000 && c<2000) //b<1000
  {
   a = a+1;
   k = randomint(norigvars);
   l = randomint(norigvars);

   for(i=0;i<nvars;i++) auxsol[i] = newsol[i];
   
   if(k != l)
     {      
       stepk = 0.00005*myrand()*((maxbound[k] - minbound[k])/b); 
       stepl = 0.00005*myrand()*((maxbound[k] - minbound[k])/b);
         
        door = randomint(5)+1;
        improve = 1;
	while(improve==1 && a<trials)
	{
         oldval_k = newsol[k];
         oldval_l = newsol[l];
         a++;
                      
         if(door==1)
	   {
            auxsol[k] = auxsol[k] +  stepk;
            auxsol[l] = auxsol[l] -  stepl;
           }
         else if(door == 2)
           {
            auxsol[k] = auxsol[k] +  stepk;
            auxsol[l] = auxsol[l] +  stepl;
           }
         else if(door == 3)
           {
            auxsol[k] = auxsol[k] -  stepk;
            auxsol[l] = auxsol[l] -  stepl;
           } 	  
         else if(door == 4)
	   {
            auxsol[l] = auxsol[l] -  stepl;         
           }
	 else
           {
            auxsol[k] = auxsol[k] +  stepk;
           }    
      
	 // Convert guarantees this check 
	 /*
         if (auxsol[k] < minbound[l] || auxsol[l] > maxbound[l])
	   {
             auxsol[k] =  minbound[l] + (maxbound[l] - minbound[l])*myrand();
           }
     
        if (auxsol[l] < minbound[l] || auxsol[l] > maxbound[l])
	   {
             auxsol[l] =  minbound[l] + (maxbound[l] - minbound[l])*myrand();
           }
	 */
        
       ConvertFromCont(auxsol);
       MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);	
       if (isnan(obj_f)) obj_f = 1000000.0; 
       eval_local++;      
       
      J = obj_f;

       if(a >= trials)
	 {  
          b = b*10;
          //a = 0;
          c = c+1;
         }
      
      
      if(J<TheBest)
       {
        for(i=0;i<nvars;i++) TheBestsol[i] = auxsol[i];
        TheBest = J;
        //cout<<"a "<<a<<" TheBest "<<TheBest<<endl;
        //for(i=0;i<vars;i++) cout<<TheBestsol[i]<<" ";
        //cout<<endl;
       }
      if (J<BestJ || ((J-BestJ)<0.0000001)) //b==1000 && 
       {  
        BestJ = J;
        for(i=0;i<nvars;i++) newsol[i] = auxsol[i];
        //a = 0;      
        //b = 1;
        
       }       
       else
       { 
        newsol[k] = oldval_k;
        newsol[l] = oldval_l;
        //auxsol[k] = oldval_k;
        //auxsol[l] = oldval_l;
        improve = 0;
       }
     }
    }  
  }

   *TBest = TheBest; 
   //cout<<setprecision(10)<<endl<<"MGA_DSM objective value = " << TheBest <<endl;
      
   return eval_local;
}




int calcfc(int n, int m, double *x, double *f, double *con, void *state_)
{
  double obj_f ;
  int i;

  /* Parameter adjustments */
  --con;
  --x;
  for(i=0;i<n;i++)
    {
      con[i+1] = x[i+1]- minbound[i];
      con[n+i+1] =  maxbound[i]- x[i+1];
    }
  ConvertFromContA(x+1);
  MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);
  if (isnan(obj_f)) obj_f = 1000000.0;
  *f=obj_f ;
  // cout<<setprecision(10)<<endl<<" eval  = " << obj_f<<endl;
  return 0;
} /* calcfc */


void LocalOpt(ContPopul* epop,int nelit, int epsize, int atgen)
{
 double obj_f;
 double* TheBestsol = new double[vars];
  int i,k,j,rc,eval_local,start;
 start = 0; 
 double rhobeg, rhoend;
 int maxfun = 50;
 int iprint = 0;
 rhobeg = 0.5;
 rhoend = 0.01;
 example_state state;
 state.nprob = 0;
 double f;

 eval_local = 0;

for(k=start; k < epsize;  k++)  
 {

   //   if(myrand()>0.5)
    
   ConvertFromCont(epop->S[k]);

   MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);	
   if (isnan(obj_f)) obj_f = 1000000.0; 

  //cout<<setprecision(10)<<endl<<"MGA_DSM objective value = " << obj_f<<endl;
       
   if(k<nelit) maxfun = 5 + fix(atgen/20);
   else maxfun = 10 + fix(atgen/10);

   //eval_local = eval_local + ContRandom(vars, epop->S[k], maxfun, &obj_f, TheBestsol);

   for(i=0; i < vars;  i++)  TheBestsol[i] = epop->S[k][i];

   rc = cobyla(vars, 2*vars, TheBestsol, rhobeg, rhoend, iprint, &maxfun, calcfc, &state);

   //for(i=0;i<vars;i++) cout<<epop->S[k][i]<<" ";
   //cout<<endl<<" After "<<endl;
   //for(i=0;i<vars;i++) cout<<TheBestsol[i]<<" ";
   //cout<<endl;

   eval_local = eval_local + maxfun;

   ConvertFromCont(TheBestsol);

   MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, f, Delta_V);
   if (isnan(f)) f = 1000000.0;
   
   if(f<obj_f)
     for(i=0; i < vars;  i++)  epop->S[k][i] = TheBestsol[i];
   else 
    {
     eval_local = eval_local + ContRandom(vars, epop->S[k], maxfun, &obj_f, TheBestsol);
    } 
          cout<<setprecision(10)<<k<<"  MGA_DSM objective value = " << obj_f<<"  "<<f<<endl;


   //if(k<nelit) eval_local = ContRandom(vars, epop->S[k], 10, &obj_f, TheBestsol);
   //else eval_local = ContRandom(vars, epop->S[k], 50, &obj_f, TheBestsol);

   /*
   cout<<" Before "<<endl;
   for(i=0;i<vars;i++) cout<<epop->S[k][i]<<" ";
   cout<<endl<<" After "<<endl;
   for(i=0;i<vars;i++) cout<<TheBestsol[i]<<" ";
   cout<<endl;
   for(i=0;i<vars;i++) epop->S[k][i]=TheBestsol[i];
   */   
   epop->SetVal(k,1000000.0 - obj_f); //Minimization is transformed in maximization 
   TotEvaluations += eval_local; 
 }
 delete[] TheBestsol;
}







void ImproveSols(ContPopul* epop,int nelit, int epsize, int atgen)
{
  double obj_f,newobj;
 double* TheBestsol = new double[vars+1];
 
 int i,k,j,eval_local,start,rc;
 //if (atgen==0) start=0;
 //else start=nelit;

 int nbett = 0;

 start = 0; 

for(k=start; k < epsize;  k++)  
 {
   ConvertFromCont(epop->S[k]);

   MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);	
   if (isnan(obj_f)) obj_f = 1000000.0; 
 
   //for(i=0;i<vars;i++) cout<<epop->S[k][i]<<" ";
   //  cout<<endl;
   cout<<setprecision(10)<<endl<<"MGA_DSM objective value = " << obj_f<<endl;   

  
   double rhobeg, rhoend;
   int maxfun = 50;
   int iprint = 0;
   rhobeg = 0.5;
   rhoend = 0.01;
   example_state state;
   state.nprob = 0;
   double f;

   //for(i=0; i < vars;  i++)  TheBestsol[i+1] = epop->S[k][i];

   rc = cobyla(vars, 2*vars, epop->S[k], rhobeg, rhoend, iprint, &maxfun, calcfc, &state);
 
   //for(i=0; i < vars;  i++)  epop->S[k][i] = TheBestsol[i+1];
  
   ConvertFromCont(epop->S[k]);
   MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, f, Delta_V);	  
   if (isnan(obj_f)) obj_f = 1000000.0; 
 
   if(f<=obj_f) nbett++;
   obj_f = f;

   //for(i=0;i<vars;i++) cout<<epop->S[k][i]<<" ";
   // cout<<endl;

    cout<<k+1<<" "<<nbett<<setprecision(10)<<" value = " << f<<endl;   


   /*
   if(obj_f<2.2) 
    {
     for(i=0;i<vars;i++) cout<<epop->S[k][i]<<" ";
     cout<<endl;
 
     cout<<k<<" "<<setprecision(10)<<obj_f;
     int j=0;
      while(j<10)
       {
	 eval_local = ContRandom(vars, epop->S[k], 100000, &newobj, TheBestsol);
        cout<<setprecision(10)<<"--->"<< newobj<<endl;
        j++;
       }
    }
   else
      eval_local = ContRandom(vars, epop->S[k], 50, &newobj, TheBestsol);
   */
   
   
  
   epop->SetVal(k,1000000.0 - obj_f); //Minimization is transformed in maximization 
   TotEvaluations += eval_local; 
 }
 delete[] TheBestsol;
}



double thefunct(double* x)
  
{
  double f,obj_f ;
   ConvertFromCont(x);
   MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);	
  f=obj_f ;


 return (f) ;
}

 

void ReadSols(ContPopul* cpop)  
 
{ 
  int i,j; 
  double auxfc;
  ifstream inFile;
  inFile.open("AllBests.txt");
  //stream = fopen( "AllBests.txt", "r" );  
  if (!inFile) {
    cerr << "Unable to open file datafile.txt";
    exit(1);   // call system to stop
 }
  
  for(i=0;i<6550;i++)
   {
     for(j=0;j<22;j++)
      {
	inFile >> auxfc;
      //cout<<auxfc<<" ";
      cpop->S[i][j] = auxfc;               
      } 
    // cout<<endl;    
   }


  cpop->RandInit();   
 ImproveSols(cpop,0, 6650, 0);
 inFile.close();
  
}






void evalcontfunction(ContPopul* epop,int nelit, int epsize, int atgen)
{
 double obj_f;
 
 int k,j,eval_local,start;
 if (atgen==0) start=0;
 else start=nelit;

 //cout<<"tract_vars is "<<tract_vars<<" length is "<<length<<endl;
 
for(k=start; k < epsize;  k++)  
 {
  
   ConvertFromCont(epop->S[k]);

   MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);	
   if (isnan(obj_f)) obj_f = 1000000.0; 
  //cout<<setprecision(10)<<endl<<"MGA_DSM objective value = " << obj_f<<endl;
   eval_local = 1;      
     
   epop->SetVal(k,1000000.0 - obj_f); //Minimization is transformed in maximization 
   TotEvaluations += eval_local; 
 }
 
}


int GaussianEDA()  //In this case, complexity is the threshold for chi-square 
{
  int i,l,ll,fgen;  
  double auxprob,sumprob;     
  IntTreeModel *IntTree;  
 
  init_time(); 
  InitContPopulations(); 
  GaussianModel* GaussEDA  = new GaussianModel(vars,psize);  
  
  i=0; auxprob =0; BestEval  = Max -1; fgen = -1;  

  TotEvaluations = 0;

  InitBounds();
  pop->RandInit();
  i = 0;
  NPoints = 10000;
   
  while (i<Maxgen && BestEval<Max && NPoints>10) 
  {   
    //evalcontfunction(pop,Elit,psize,i);
    LocalOpt(pop,Elit,psize,i);
 
     NPoints = Selection(); 
     GaussEDA->LearnModel(NPoints,selpop);       
     FindBestVal(); 
        

      if(printvals>1)
       {
	 for(ll=0;ll<printvals-1;ll++)// NPoints
	   {
	     cout<<"BestD"<<ll<<": ";
	     for(l=0;l<vars;l++) cout<<selpop->S[ll][l]<<" ";
	     cout<<endl;
	     cout<<"BestV"<<ll<<": "<<setprecision(10)<<" "<<1000000.0-selpop->Evaluations[ll]<<endl;
	     cout<<endl;
	   }
	 if(printvals)  cout<< " Gen : "<<i<<" Best: "<<1000000.0-BestEval<<" ProbBest: "<<auxprob<<" DifPoints: \
"<<NPoints<<" TotEvaluations  "<<TotEvaluations <<endl;
 
  
       }
    

      if (BestEval>=Max)   fgen  = i;	 
      else 
          { 
           if (Tour>0 || (Tour==0 && Elit>TruncMax)) elitpop->SetElit(Elit,pop);             
           else 
              {
		selpop->SetElit(Elit,pop);
	        for(ll=0;ll<Elit;ll++)   pop->Evaluations[ll]=selpop->Evaluations[ll];               
               }
          }

      if(NPoints>10) GaussEDA->GenPop(Elit,psize,pop);   
    
     i++;
   }  
  end_time(); 
  if(printvals>0)  cout<<"LastGen : "<<i<<" Best: "<<1000000.0-BestEval<<" ProbBest: "<<auxprob<<" DifPoints: "<<NPoints<<" TotEval: "<<TotEvaluations<<" time "<<auxtime<<endl;  //cout<<BestEval<<endl; 
  if(NPoints>10) NPoints = 10;
   for(ll=0;ll<printvals-1;ll++)// NPoints
	   {
	     cout<<"LastBestD"<<ll<<": ";
	     for(l=0;l<vars;l++) cout<<selpop->S[ll][l]<<" ";
	     cout<<endl;
	     cout<<"LastBestV"<<ll<<": "<<setprecision(10)<<" "<<1000000.0-selpop->Evaluations[ll]<<endl;
	     cout<<endl;
	   }
 
  DeleteContPopulations(); 
  delete GaussEDA;
  return fgen;  

}  



int MixtureGaussianEDA(int nclusters)  //In this case, complexity is the threshold for chi-square 
{
  int i,l,ll,fgen;  
  double auxprob,sumprob;     
  IntTreeModel *IntTree;  
 
  init_time(); 
  InitContPopulations();
 
  MixtureOfGaussians* MixtGaussEDA  = new MixtureOfGaussians(vars,psize,nclusters);  
 


  i=0; auxprob =0; BestEval  = Max -1; fgen = -1;  

  TotEvaluations = 0;

  InitBounds();
  pop->RandInit();
  i = 0;
  NPoints = 10000;
   
  //MixtGaussEDA->ExamplesClustering(pop->S);
  while (i<Maxgen && BestEval<Max && NPoints>=10) 
  {   
    //evalcontfunction(pop,Elit,psize,i);
    LocalOpt(pop,Elit,psize,i);
 
     NPoints = Selection(); 
     MixtGaussEDA->LearnModel(selpop,NPoints,10);       
     FindBestVal(); 
    
      if(printvals>1)
       {

     for(ll=0;ll<printvals-1;ll++)// NPoints
	   {
	     cout<<"BestD"<<ll<<": ";
	     for(l=0;l<vars;l++) cout<<selpop->S[ll][l]<<" ";
	     cout<<endl;
	     cout<<"BestV"<<ll<<": "<<setprecision(10)<<" "<<1000000.0-selpop->Evaluations[ll]<<endl;
	     cout<<endl;
	   }

	 if(printvals)  cout<< " Gen : "<<i<<" Best: "<<1000000.0-BestEval<<" ProbBest: "<<auxprob<<" DifPoints: \
"<<NPoints<<" TotEvaluations  "<<TotEvaluations <<endl;
 
  
       }
    

      if (BestEval>=Max)   fgen  = i;	 
      else 
          { 
           if (Tour>0 || (Tour==0 && Elit>TruncMax)) elitpop->SetElit(Elit,pop);             
           else 
              {
		selpop->SetElit(Elit,pop);
	        for(ll=0;ll<Elit;ll++)   pop->Evaluations[ll]=selpop->Evaluations[ll];               
               }
          }

     
      if(NPoints>10) MixtGaussEDA->GenPop(Elit,psize,pop);   

     
     i++;

   }  
  end_time(); 
  if(printvals>0)  cout<<"LastGen : "<<i<<" Best: "<<1000000.0-BestEval<<" ProbBest: "<<auxprob<<" DifPoints: "<<NPoints<<" TotEval: "<<TotEvaluations<<" time "<<auxtime<<endl;  //cout<<BestEval<<endl; 
 for(ll=0;ll<printvals-1;ll++)// NPoints
	   {
	     cout<<"LastBestD"<<ll<<": ";
	     for(l=0;l<vars;l++) cout<<selpop->S[ll][l]<<" ";
	     cout<<endl;
	     cout<<"LastBestV"<<ll<<": "<<setprecision(10)<<" "<<1000000.0-selpop->Evaluations[ll]<<endl;
	     cout<<endl;
	   }

 
if(NPoints>10) NPoints = 10;
  
 cout<<" Pass 01 "<<endl; 
  DeleteContPopulations(); 
 cout<<" Pass 02 "<<endl; 
  delete MixtGaussEDA;
 cout<<" Pass 03 "<<endl; 
  return fgen;  

}  





void PrintStatistics() 
{  
  int i;
  double auxmeangen,meanfit,sigma; 
 
  sigma = 0;
                   meaneval /=  cantexp; 
                   alltime  =  alltime/(1.0*cantexp); 
		   for (i=0;i<cantexp;i++) 
                   {
                    sigma += (meanlikehood[i] - meaneval)*(meanlikehood[i] - meaneval);
                    //cout<<sigma<<endl;
                   } 
                   sigma = sigma/(cantexp-1);
                   
                  if (succexp>0)  
                   {  
                    auxmeangen = meangen/succexp;
                    bestalltime = bestalltime/(1.0*succexp); 
                    if (BestElitism)  
                         meanfit = (auxmeangen+1)*(1-Trunc)*psize + psize*Trunc;     
                    else meanfit = (auxmeangen+1)*(psize-1) + 1; 
                    cout<<"TypeExp="<<ExperimentMode<<"  n="<<vars<<" T="<<Trunc<<" N="<<psize<<" Sel="<<Tour<<"  k="<<Cycles<<"  MaxGen="<<Maxgen<<"  Elit="<<Elit<<" Suc.="<<succexp<<"  g="<<(auxmeangen+1)<<"  ave="<<meanfit<<" meaneval "<<meaneval<<" sigma "<<sigma<<" timebest "<<bestalltime<<" fulltime "<<alltime<<endl;                   
                   } 
                  else  
                   {  
		     cout<<"TypeExp="<<ExperimentMode<<"  n="<<vars<<" T="<<Trunc<<" N="<<psize<<" Sel="<<Tour<<"  k="<<Cycles<<"  MaxGen="<<Maxgen<<"  Elit="<<Elit<<" Suc.="<<0<<"  g="<<0<<"  ave="<<0<<" meaneval "<<meaneval<<" sigma "<<sigma<<" fulltime "<<alltime<<" Eval "<<(TotEvaluations/(1.0*cantexp))<<endl; 
                   } 

		  //for(int ll=0;ll<Maxgen;ll++)  cout<<AllGen[ll]/(-1.0*cantexp)<<" ";
                  //cout<<endl;
} 



void runOptimizer(int algtype,int nrun)  
{  
    int succ=-1; 
        
  switch(algtype)  
                     {                     
                       case 5: succ =  GaussianEDA();break;  // Gaussian model
                       case 6: succ =  MixtureGaussianEDA(10);break;  // Mixture of Gaussian models 
                     }                      

  

   if (succ>-1)  
   { 
       succexp++; 
       meangen += succ;    
       bestalltime +=auxtime;      
        
   } 
   else nsucc++;
   alltime += auxtime;  
   meaneval += BestEval; 
   meanlikehood[nrun] = BestEval;  
} 



int main( int argc, char *argv[] )
{
   
 
  // ./trajectory 1 5 3 12 1000 50 250 1 1000000 3 0 1 1 2 2   
 
 int i,a;
  int protein_inst,modeprotein;
  int T,MaxMixtP,S_alph,Compl; 

  
 if( argc != 16 ) {
    std::cout << "Usage: " <<"cantexp  EDA{0:Markov, 1:Tree  2:Mixture, 4:AffEDA} modeprotein{2,3} prot_inst n psize Trunc max-gen" << std::endl;
    std::cout << "       Please read the README file." << std::endl;
    exit(1);
}

 params = new int[3];    
 cantexp = atoi(argv[1]);         // Number of experiments
 ExperimentMode = atoi(argv[2]); // Type of EDA
 
 length = atoi(argv[3]);             // Number of bits for  each variable   
 vars =  atoi(argv[4]);            //Number of variables (redundant because depends on instance)
 psize = atoi(argv[5]);          // Population size
 T = atoi(argv[6]);              // Percentage of truncation integer number (1:99)
 Maxgen =  atoi(argv[7]);        // Max number of generations 
 BestElitism = atoi(argv[8]);         // If there is or not BestElitism, if thereisnot BestElitism, Elitism = 1 by default;
 Max = atoi(argv[9]);
 CliqMaxLength = atoi(argv[10]);
 func  = atoi(argv[11]);   // Function to be used  0) Simple evaluation 1) Local optimizer 
 params[0]  = atoi(argv[12]);
 params[1]  = atoi(argv[13]);
 printvals  = atoi(argv[14]);
 Card =  atoi(argv[15]);
 
 


 Tour = 0;                       // Truncation Selection is used
 Ntrees = 2;                     // Number of Trees  for MT-EDA
 Elit = 1;                       // Elitism
 Nsteps = 50;                    // Learning steps of the Mixture Algorithm  
 InitTreeStructure = 1;    // 0 for a random init tree structures, 1 for a Chu&Liu learned init Tree Structure  
 VisibleChoiceVar = 0;     // 0 The Mixture choice var is hidden, 1 & 2 depends on a variable and the  unitation respectively  
 //printvals = 2;            // The printvals-1 best values in each generation are printed 
 MaxMixtP = 500;           // Maximum learning parameter mixture 
 S_alph = 0;               // Value alpha for smoothing 
 StopCrit = 1;             // Stop Criteria for Learning of trees alg.  
 Prior = 1;                // Type of prior. 
 Compl=75;                 // Complexities of the trees. 
 Coeftype=2;               // Type of coefficient calculation for Exact Learning. 
 //params[0] = 3 ;           //  Params for function evaluation 
 //params[1] = 3;  
 params[2] = 10;  
 
 
 //seed =  1243343896; 
 seed = (unsigned) time(NULL);  
 srand(seed); 
 cout<<"seed"<<seed<<endl; 

TypeMixture = 1; // Class of MT-FDA (1-Meila, 2-MutInf)
Mutation = 0 ; // Population based mutation  
//CliqMaxLength = 2; // Maximum size of the cliques for Markov  or maximum number of neighbors for MOA
MaxNumCliq = 300; // Maximum number of cliques for Markov 
OldWaySel = 0; // Selection with sel pop (1) or straight on Sel prob (0) 
LearningType = 6; // Learning for MNFDA (0-Markov, 1-JuntionTree) 
Cycles = 0 ; // Number of cycles for GS in the MNEDA or size for the clique in Markov EDA. 
Trunc = T/double(100);  
Complex  = Compl/double(100);  
MaxMixtProb =MaxMixtP/double(100); 
S_alpha = S_alph/double(100); 


int k,j,u;
double eval;



 InitTrajectory();  

 
  cout<<"Alg : "<<ExperimentMode<<", number codifying bits : "<<length<<", n : "<<vars<<", psize : "<<psize<<", Trunc : "<<T<<", max-gen : "<<Maxgen<<", BestElit. : "<<BestElitism<<", NNeighbors  : "<<CliqMaxLength<<", Card  : "<<Card<<", func : "<<func<<", params[0] : "<<params[0]<<", params[1] : "<<params[1]<<endl; 

        AbsBestInd = new double [vars];
        AbsBestEval = -1;
        TotEvaluations = 0;       
       	succexp = 0;  meangen = 0; meaneval = 0;  i =0;  nsucc =0; alltime = 0; bestalltime = 0;  
	while (i<cantexp) //&& nsucc<1
        { 
          currentexp = i;	  
	  runOptimizer(ExperimentMode,i);
          i++;
         //PrintStatistics();
        }  


	//ContPopul* cpop =  new ContPopul(7000,vars,Elit,minbound,maxbound); 
        //ReadSols(cpop);  
	//delete cpop;

         
         PrintStatistics();             
	delete[] AbsBestInd;   
     
      
 delete [] params; 
 delete[] Delta_V;
return 0;

}      




