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
#include "asa_cg.h"
//#include "asa_user.h"



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


/* prototypes for the function and gradient evaluation routines */










int asa_cg /*  return:
                      -2 (function value became nan in cg)
                      -1 (starting function value is nan in cg)
                       0 (convergence tolerance satisfied)
                       1 (change in func <= feps*|f| in cg)
                       2 (cg iterations in all passes or
                          in one pass exceeded their limit)
                       3 (slope always negative in line search in cg)
                       4 (number secant iterations exceed nsecant in cg)
                       5 (search direction not a descent direction in cg)
                       6 (line search fails in initial interval in cg)
                       7 (line search fails during bisection in cg)
                       8 (line search fails during interval update in cg)
                       9 (debugger is on and the function value increases in cg)
                      10 (out of memory)
                      11 (cbb iterations in all passes or
                          in one pass exceeded their limit)
                      12 (line search failed in cbb iteration)
                      13 (search direction in cbb is not descent direction)
                      14 (function value became nan in cbb) */
(
    double            *x, /* input: starting guess, output: the solution */
    double           *lo, /* lower bounds */
    double           *hi, /* upper bounds */
    INT                n, /* problem dimension */
    asa_stat       *Stat, /* structure with statistics (can be NULL) */
    asacg_parm    *CParm, /* user parameters, NULL = use default parameters */
    asa_parm      *AParm, /* user parameters, NULL = use default parameters */
    double      grad_tol, /* |Proj (x_k - g_k) - x_k|_inf <= grad_tol */
    double   (*value) (asa_objective *), /* evaluate objective function */
    void      (*grad) (asa_objective *), /* evaluate objective gradient */
    double (*valgrad) (asa_objective *), /* function and gradient
                                            NULL = use value & grad routines */
    double        *Work  /* either work array of size 7n + memory (m) or NULL */
)
{
    int gp, ident, j, nfree, status, *ifree ;
    INT cbb_totit, cg_totit ;
    double alpha, gj, pert_lo, pert_hi, t, tl, th, gnorm, ginorm, pgnorm, xnorm,
           xj, xg, xp, *work, *d, *g, *xtemp, *gtemp, *pg ;
    asacg_parm *cgParm, cgParmStruc ;
    asa_parm *asaParm, asaParmStruc ;
    asa_com Com ;
    asa_objective Objective ;

/* initialize the parameters */

    if ( CParm == NULL )
    {
        cgParm = &cgParmStruc ;
        asa_cg_default (cgParm) ;
    }
    else cgParm = CParm ;
    if ( cgParm->PrintParms ) asa_printcgParms (cgParm) ;

    if ( AParm == NULL )
    {
        asaParm = &asaParmStruc ;
        asa_default (asaParm) ;
    }
    else asaParm = AParm ;
    if ( asaParm->PrintParms ) asa_printParms (asaParm) ;

    /* abort after maxit iterations of cbb in one pass */
    if ( asaParm->maxit_fac == INF ) Com.pgmaxit = INT_INF ;
    else Com.pgmaxit = (INT) (((double) n)*asaParm->maxit_fac) ;

    /* abort after totit iterations of cbb in all passes */
    if ( asaParm->totit_fac == INF ) cbb_totit = INT_INF ;
    else cbb_totit = (INT) (((double) n)*asaParm->totit_fac) ;

    /* abort after maxfunc function evaluation in one pass of cbb */
    if ( asaParm->maxfunc_fac == INF ) Com.pgmaxfunc = INT_INF ;
    else Com.pgmaxfunc = (INT) (((double) n)*asaParm->maxfunc_fac) ;

    /* abort after totit iterations of cg in all passes */
    if ( cgParm->totit_fac == INF ) cg_totit = INT_INF ;
    else cg_totit = (INT) (((double) n)*cgParm->totit_fac) ;

    pert_lo = asaParm->pert_lo ;
    pert_hi = asaParm->pert_hi ;
    Com.user = &Objective ;
    Objective.n = n ;
    Com.tau1 = asaParm->tau1 ;
    Com.tau2 = asaParm->tau2 ;

    Com.cgParm = cgParm ;
    Com.asaParm = asaParm ;
    Com.x = x ;
    Com.n = n ;             /* problem dimension */
    Com.n5 = n % 5 ;
    Com.nf = (INT) 0 ;      /* number of function evaluations */
    Com.ng = (INT) 0 ;      /* number of gradient evaluations */
    Com.cbbiter = (INT) 0 ; /* number of cbb iterations evaluations */
    Com.cgiter = (INT) 0 ;  /* number of cg iterations */
    Com.AWolfe = cgParm->AWolfe ; /* do not touch user's AWolfe */
    Com.AArmijo = asaParm->AArmijo ; /* do not touch user's AArmijo */
    Com.value = value ;
    Com.grad = grad ;
    Com.valgrad = valgrad ;
    Com.DimReduce = FALSE ;
    ifree =  Com.ifree = (int*) malloc (n*sizeof (int)) ;

    if ( Work == NULL ) work = (double*) malloc ((5*n+asaParm->m)*sizeof (double)) ;
    else                work = Work ;
    if ( work == NULL )
    {
        printf ("Insufficient memory for specified problem dimension %e\n",
                 (double) n) ;
        status = 10 ;
        return (status) ;
    }
    d = Com.d = work ;
    g = Com.g = d+n ;
    xtemp = Com.xtemp = g+n ;
    gtemp = Com.gtemp = xtemp+n ;
    pg = Com.pg = gtemp+n ;
    Com.lastfvalues = pg+n ; /* size asaParm->m */
    Com.lo = lo ;
    Com.hi = hi ;
    Com.cbbiter = 0 ;
    Com.cbbfunc = 0 ;
    Com.cbbgrad = 0 ;
    Com.cgiter = 0 ;
    Com.cgfunc = 0 ;
    Com.cggrad = 0 ;

    ident = FALSE ;
    xnorm = ZERO ;
    for (j = 0; j < n; j++)
    {
        t = x [j] ;
        if      ( t > hi [j] ) t = hi [j] ;
        else if ( t < lo [j] ) t = lo [j] ;
        x [j] = t ;
        if ( xnorm < fabs (t) ) xnorm = fabs (t) ;
    }

    Com.f = asa_fg (g, x, &Com) ;
    pgnorm = ZERO ;
    gnorm = ZERO ;
    for (j = 0; j < n; j++)
    {
        xj = x [j] ;
        gj = g [j] ;
        xg = xj - gj ;
        if      ( xg > hi [j] ) xp = hi [j] - xj ;
        else if ( xg < lo [j] ) xp = lo [j] - xj ;
        else                    xp = -gj ;
        pg [j] = xp ;
        pgnorm = MAX (pgnorm, fabs (xp)) ;
        gnorm = MAX (gnorm, fabs (gj)) ;
    }
    if ( asaParm->StopRule ) Com.tol = MAX (pgnorm*asaParm->StopFac, grad_tol) ;
    else                     Com.tol = grad_tol ;

    Com.pgnorm = Com.pgnorm_start = pgnorm ;
    if ( asa_tol (pgnorm, &Com) )
    {
        status = 0 ;
        goto Exit ;
    }

    if ( xnorm != ZERO ) Com.alpha = alpha = xnorm/gnorm ;
    else                 Com.alpha = alpha = ONE/gnorm ;

    /* compute gradient norm for inactive variables */
    ginorm = ZERO ;
    nfree = 0 ;
    gp = FALSE ;
    for (j = 0; j < n; j++)
    {
        xj = x [j] ;
        tl = lo [j] ;
        th = hi [j] ;
        gj = g [j] ;
        xg = xj - alpha*gj ;
        if      ( (xg >= th) && (th-xj > pert_hi) ) gp = TRUE ;
        else if ( (xg <= tl) && (xj-tl > pert_lo) ) gp = TRUE ;
        if ( (xj-tl > pert_lo) && (th - xj > pert_hi) )
        {
            ginorm = MAX (ginorm, fabs (gj)) ;
            ifree [nfree] = j ;
            nfree++ ;
        }
    }
    Com.ginorm = ginorm ;
    Com.nfree = nfree ;

    if ( asaParm->PrintLevel >= 1 )
    {
        printf ("\ninitial f = %14.6e pgnorm = %14.6e ginorm = %14.6e\n",
                 Com.f, pgnorm, ginorm) ;
        printf ("            nfree = %i xnorm = %14.6e gp = %i\n",
                 nfree, xnorm, gp) ;
    }

    if ( (ginorm < Com.tau1*pgnorm) || gp || asaParm->GradProjOnly )
    {
        Com.cbbfunc = 1 ;
        Com.cbbgrad = 1 ;
        goto Grad_proj ;
    }
    else
    {
        Com.cgfunc = 1 ;
        Com.cggrad = 1 ;
        goto CG_descent ;
    }

    Grad_proj:
    if ( asaParm->PrintLevel >= 1 ) printf ("\nGradProj:\n") ;
    Com.DimReduce = FALSE ;
    status = asa_grad_proj(&Com) ;
    if ( asaParm->PrintLevel >= 1 )
    {
        printf ("exit Grad_proj\n") ;
    }
    if ( Com.cbbiter >= cbb_totit ) status = 11 ;
    if ( status >= 0 ) goto Exit ;

    /* extract free variable */
    nfree = 0 ;
    for (j = 0; j < n; j++)
    {
        xj = x [j] ;
        if ( (xj-lo [j] > pert_lo) && (hi [j] - xj > pert_hi) )
        {
            ifree [nfree] = j ;
            nfree++ ;
        }
    }
    Com.nfree = nfree ;

    CG_descent:
    if ( nfree != n )
    {
       asa_shrink_all (&Com) ;
       asa_copy (xtemp+nfree, x+nfree, n-nfree) ;
       Com.DimReduce = TRUE ;
    }
    else Com.DimReduce = FALSE ;

    if ( asaParm->PrintLevel >= 1 ) printf ("\nCG:\n") ;
    status = asa_descent (&Com) ;

    if ( asaParm->PrintLevel >= 1 )
    {
        printf ("exit the CG subroutine\n") ;
    }
    if ( Com.DimReduce ) asa_expand_all (&Com) ;
    if ( Com.cgiter >= cg_totit ) status = 2 ;

    if ( status >= -2 ) goto Exit ;

    /* ginorm < tau2* pgnorm without hitting boundary */
    if ( status == -5 )
    {
        Com.alpha = asa_init_bbstep (&Com) ;
        goto Grad_proj ;

    }
    /* ginorm >= tau2* pgnorm and many components of x hit boundary  */
    else if ( status == -4 )
    {
        ginorm = ZERO ;
        nfree = 0 ;
        for (j = 0 ; j < n; j++)
        {
            xj = x [j] ;
            if ( (xj-lo [j] > pert_lo) && (hi [j] - xj > pert_hi) )
            {
                t = fabs (g [j]) ;
                ginorm = MAX (ginorm, t) ;
                ifree [nfree] = j ;
                nfree++ ;
            }
        }
        Com.nfree = nfree ;
        Com.ginorm = ginorm ;

        if ( ginorm >= Com.tau1*Com.pgnorm ) goto CG_descent ;
        else
        {
           if ( asaParm->PrintLevel >= 1 ) printf ("ginorm < tau1* pgnorm\n") ;
           Com.alpha = asa_init_bbstep (&Com) ;
           goto Grad_proj ;
        }
    }
    /* ginorm >= tau2* pgnorm and only one component of x hits boundary */
    else if ( status == -3 )
    {
        if ( pgnorm < asaParm->pgdecay*MAX (Com.pgnorm_start, ONE) )
        {
            ident = asa_identify (x, g, Com.pgnorm, &Com) ;
        }
        if ( ident )
        {
            ident = FALSE ;
            ginorm = ZERO ;
            nfree = 0 ;
            for (j = 0 ; j < n; j++)
            {
                xj = x [j] ;
                if ( (xj-lo [j] > pert_lo) && (hi [j] - xj > pert_hi) )
                {
                    t = fabs (g [j]) ;
                    ginorm = MAX (ginorm, t) ;
                    ifree [nfree] = j ;
                    nfree++ ;
                }
            }
            Com.nfree = nfree ;
            Com.ginorm = ginorm ;
            if ( ginorm >= Com.tau1*Com.pgnorm ) goto CG_descent ;
            else
            {
               if ( asaParm->PrintLevel >= 1 )
                   printf ("ginorm < tau1* pgnorm\n" ) ;
               Com.alpha = asa_init_bbstep (&Com) ;
               goto Grad_proj ;
            }
        }
        else
        {
            Com.alpha = asa_init_bbstep (&Com) ;
            goto Grad_proj ;
        }
    }

    Exit:
    if ( (asaParm->PrintFinal) || (asaParm->PrintLevel >= 1) )
    {
        const char mess1 [] = "Possible causes of this error message:" ;
        const char mess2 [] = "   - your tolerance may be too strict: "
                              "grad_tol = " ;
        const char mess4 [] = "   - your gradient routine has an error" ;
        const char mess5 [] = "   - the parameter epsilon in "
                              "asa_descent_c.parm is too small" ;
        printf ("\nFinal convergence status = %d\n", status);
        if ( status == -2 )
        {
            printf ("Function value became nan at cg iteration %10.0e\n",
                     (double) Com.cgiter) ;
        }
        else if ( status == -1 )
        {
            printf ("Function value of starting point is nan at "
                     "cg iteration %10.0f\n", (double) Com.cgiter) ;
        }
        else if ( status == 0 )
        {
            printf ("Convergence tolerance for gradient satisfied\n") ;
        }
        else if ( status == 1 )
        {
            printf ("Terminating in cg since change in function value "
                    "<= feps*|f|\n") ;
        }
        else if ( status == 2 )
        {
            printf ("Number of iterations exceed specified limits "
                    "for cg routine\n") ;
            printf ("Iterations: %10.0f maxit: %10.0f totit: %10.0f\n",
                    (double) Com.cgiter, (double) Com.cgmaxit,
                    (double) cg_totit) ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
        }
        else if ( status == 3 )
        {
            printf ("Slope always negative in cg line search\n") ;
            printf ("%s\n", mess1) ;
            printf ("   - your cost function has an error\n") ;
            printf ("%s\n", mess4) ;
        }
        else if ( status == 4 )
        {
            printf ("Line search fails in cg, too many secant steps\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
        }
        else if ( status == 5 )
        {
            printf ("Search direction not a descent direction in cg\n") ;
        }
        else if ( status == 6 ) /* line search fails */
        {
            printf ("Line search fails in cg iteration\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
            printf ("%s\n", mess4) ;
            printf ("%s\n", mess5) ;
        }
        else if ( status == 7 ) /* line search fails */
        {
            printf ("Line search fails in cg iteration\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
        }
        else if ( status == 8 ) /* line search fails */
        {
            printf ("Line search fails in cg iteration\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
            printf ("%s\n", mess4) ;
            printf ("%s\n", mess5) ;
        }
        else if ( status == 9 )
        {
            printf ("Debugger is on, function value does not improve in cg\n") ;
            printf ("new value: %25.16e old value: %25.16e\n",
                Com.f_debug, Com.f0) ;
        }
        else if ( status == 10 )
        {
            printf ("Insufficient memory\n") ;
        }
        else if ( status == 11 )
        {
            printf ("Number of iterations or function evaluation exceed "
                          "specified limits for cbb routine\n") ;
            printf ("Iterations: %i maxit: %i totit: %i\n",
                     Com.cbbiter, Com.pgmaxit, cbb_totit) ;
            printf ("Total function evaluations: %i maxfunc: %i\n",
                     Com.nf, Com.pgmaxfunc);
        }
        if ( status == 12 ) /* line search fails in cbb iteration */
        {
            printf ("Line search fails in cbb iteration\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, Com.tol) ;
            printf ("%s\n", mess4) ;
        }

        if ( status == 13 )
        {
            printf ("Search direction not descent direction in "
                    "asa_grad_proj\n") ;
            printf ("directional derivative: %e\n", Com.gtd) ;
        }
        if ( status == 14 )
        {
             printf ("At cbb iteration %i function value became nan\n",
                      Com.cbbiter) ;
        }

        printf ("projected gradient max norm: %13.6e\n", Com.pgnorm) ;
        printf ("function value:              %13.6e\n", Com.f) ;
        printf ("\nTotal cg  iterations:           %10.0f\n",
                (double) Com.cgiter) ;
        printf ("Total cg  function evaluations: %10.0f\n",
                (double) Com.cgfunc) ;
        printf ("Total cg  gradient evaluations: %10.0f\n",
                (double) Com.cggrad) ;
        printf ("Total cbb iterations:           %10.0f\n",
                (double) Com.cbbiter) ;
        printf ("Total cbb function evaluations: %10.0f\n",
                (double) Com.cbbfunc) ;
        printf ("Total cbb gradient evaluations: %10.0f\n",
                    (double) Com.cbbgrad) ;
        printf ("------------------------------------------\n") ;
        printf ("Total function evaluations:     %10.0f\n",
                (double) Com.nf) ;
        printf ("Total gradient evaluations:     %10.0f\n",
                (double) Com.ng) ;
        printf ("==========================================\n\n") ;
    }
    free (ifree) ;
    if ( Work == NULL ) free (work) ;
    if ( Stat != NULL )
    {
        Stat->f = Com.f ;
        Stat->pgnorm = Com.pgnorm ;
        Stat->cgiter = Com.cgiter ;
        Stat->cgfunc = Com.cgfunc ;
        Stat->cggrad = Com.cggrad ;
        Stat->cbbiter = Com.cbbiter ;
        Stat->cbbfunc = Com.cbbfunc ;
        Stat->cbbgrad = Com.cbbgrad ;
    }
    return (status) ;
}

/* =========================================================================
   === asa_default ======================================================
   =========================================================================
   Set default parameter values for the ASA routine. The CG default
   parameter values are set by asa_cg_default.  If the parameter argument of
   asa_descent is NULL, this routine is called by asa_cg automatically.
   If the user wishes to set parameter values, then the asa_parameter structure
   should be allocated in the main program. The user could call asa_default
   to initialize the structure, and then individual elements in the structure
   could be changed, before passing the structure to asa_cg.
   =========================================================================*/
void asa_default
(
    asa_parm *Parm
)
{
    double eps, t ;

    /* T => print final statistics
       F => no printout of statistics */
    Parm->PrintFinal = TRUE ;

    /* Level 0  = no printing), ... , Level 4 = maximum printing */
    Parm->PrintLevel = 0 ;

    /* T => print parameters values
       F => do not display parameter values */
    Parm->PrintParms = FALSE ;

    /* T => use approximate nonmonotone Armijo line search
       F => use ordinary nonmonotone Armijo line search, switch to
            approximate Armijo when |f_r-f| < AArmijoFac*|min (f_r, f_{max})| */
    Parm->AArmijo = FALSE ;
    Parm->AArmijoFac = 1.e-8 ;

    /* Stop Rules (these override the corresponding cg parameters):
       T => ||proj_grad||_infty <= max(grad_tol,initial ||grad||_infty*StopFac)
       F => ||proj_grad||_infty <= grad_tol*(1 + |f_k|) */
    Parm->StopRule = TRUE ;
    Parm->StopFac = 0.e-12 ;

    /* T => estimated error in function value = eps*|min (f_r, f_{max}) |
       F => estimated error in function value = eps */
    Parm->PertRule = TRUE ;
    Parm->eps = 1.e-6 ;

    /* T => only use gradient projection algorithm
       F => let algorithm decide between grad_proj and cg_descent */
    Parm->GradProjOnly = FALSE ;

    /* maximum number of times the Armijo line search will perform
       backtracking steps */
    Parm->max_backsteps = (int) 50 ;

    /* abort cbb after maxit_fac*n iterations in one pass through cbb */
    Parm->maxit_fac = INF ;

    /* abort cbb after totit_fac*n iterations in all passes through cbb */
    Parm->totit_fac = INF ;

    /* abort cbb iteration after maxfunc_fac*n function evaluations */
    Parm->maxfunc_fac = INF ;

    /* perturbation in bounds based on machine epsilon, which we now compute */
    eps = ONE ;
    t = ONE ;
    while ( t > 0 )
    {
        eps /= TWO ;
        t = ONE + eps ;
        t -= ONE ;
    }
    eps *= 2 ;                   /* machine epsilon */
    Parm->pert_lo = 1.e3*eps ;   /* perturbation of lower bounds */
    Parm->pert_hi = 1.e3*eps ;   /* perturbation of upper bounds */

    /* search for non nan function value by shrinking search interval
       at most nshrink times */
    Parm->nshrink = (int) 50 ;

    /* factor by which interval shrinks when searching for non nan value */
    Parm->nan_fac = 2.e-1 ;

    /* update fr if fmin was not improved after L iterations */
    Parm->L = 3 ;

    /* fmax = max (f_{k-i}, i = 0, 1, ..., min (k, m-1) ) */
    Parm->m = 8 ;

    /* update fr if initial stepsize was accepted in previous P iterations */
    Parm->P = 40 ;

    /* CBB cycle length */
    Parm->nm = 4 ;

    /* Reinitialize BB stepsize, if (s^t y)/(||s|| ||y||) >= gamma
       and ||s|| <= min (parm3*|f_k+1|/||g_k+1||_infty, 1) */
    Parm->gamma = 0.975e0 ;

    /* update reference value fr if (fr-fmin)/(fc-fmin) > gamma1 */
    Parm->gamma1 = (double) Parm->m / (double) Parm->L ;

    /* update fr if (fr-f)/(fmax-f) > gamma2, np > P, and fmax > f */
    Parm->gamma2 = (double) Parm->P / (double) Parm->m ;

    /* terminate Armijo line search when
       phi(alpha) <= phi_r + alpha * delta * phi'(0) where phi_r = fr or fcomp*/
    Parm->delta = 1.0e-4 ;   /* Armijo line search parameter */

    /* stepsize s in the line search must satisfy lmin <= s <= lmax */
    Parm->lmin = 1.0e-20 ;
    Parm->lmax = 1.0e+20 ;

    /* attempt a quadratic interpolation step in cg_descent if the
       provisional stepsize times parm1 <= stepsize to boundary */
    Parm->parm1 = 1.e-1 ;

    /* if quadratic interpolation step is attempted, the provisional step
       is at most parm2*stepsize to boundary */
    Parm->parm2 = 9.e-1 ;

    /* used in the the criterion of reinitializing the BB stepsize */
    Parm->parm3 = 1.e-1 ;

    /* maximum number of previous BB steps used when s^t y <= ZERO */
    Parm->parm4 = 6 ;

    /* if ginorm < tau1*pgnorm, continue gradient projection steps  */
    Parm->tau1 = 1.e-1 ;

    /* decay factor for tau1 */
    Parm->tau1_decay = 5.e-1 ;

    /* ginorm < tau2*pgnorm implies subproblem solved in cgdescent */
    Parm->tau2 = 1.e-1 ;

    /* decay factor for tau2 */
    Parm->tau2_decay = 5.e-1 ;

    /* if pgnorm < pgdecay*MAX (pgnorm0, ONE), check the undecided index set
                                pgnorm0 = pgnorm at starting point */
    Parm->pgdecay = 1.e-4 ;

    /* backtracking decay factor in the Armijo line search */
    Parm->armijo_decay = 5.e-1 ;

    /* use quadratic interpolation to compute Armijo step if it
       lies in the interval [.1 alpha, .9 alpha] */
    Parm->armijo0 = 1.e-1 ;
    Parm->armijo1 = 9.e-1 ;
}

/* =========================================================================
   === asa_cg_default ======================================================
   =========================================================================
   Set default conjugate gradient parameter values. If the parameter argument
   of asa_cg is NULL, this routine is called by asa_cg automatically.
   If the user wishes to set parameter values, then the asa_parameter structure
   should be allocated in the main program. The user could call asa_cg_default
   to initialize the structure, and then individual elements in the structure
   could be changed, before passing the structure to asa_cg.
   =========================================================================*/
void asa_cg_default
(
    asacg_parm   *Parm
)
{
    /* Level 0 = no printing, ... , Level 4 = maximum printing */
    Parm->PrintLevel = 0 ;

    /* T => print parameters values
       F => do not display parameter values */
    Parm->PrintParms = FALSE ;

    /* T => use approximate Wolfe line search
       F => use ordinary Wolfe line search, switch to approximate Wolfe when
                |f_k+1-f_k| < AWolfeFac*C_k, C_k = average size of cost */
    Parm->AWolfe = FALSE ;
    Parm->AWolfeFac = 1.e-3 ;

    /* T => estimated error in function value is eps*Ck,
       F => estimated error in function value is eps */
    Parm->PertRule = TRUE ;
    Parm->eps = 1.e-6 ;

    /* T => attempt quadratic interpolation in line search when
                |f_k+1 - f_k|/f_k <= QuadCutOff
       F => no quadratic interpolation step */
    Parm->QuadStep = TRUE ;
    Parm->QuadCutOff = 1.e-12 ;

    /* T => check that f_k+1 - f_k <= debugtol*C_k
       F => no checking of function values */
    Parm->debug = FALSE ;
    Parm->debugtol = 1.e-10 ;

    /* factor in [0, 1] used to compute average cost magnitude C_k as follows:
       Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
    Parm->Qdecay = .7 ;

    /* if step is nonzero, it is the initial step of the initial line search */
    Parm->step = ZERO ;

    /* abort cg after maxit_fac*n iterations in one pass */
    Parm->maxit_fac = INF ;

    /* abort cg after totit_fac*n iterations in all passes */
    Parm->totit_fac = INF ;

    /* maximum number of times the bracketing interval grows or shrinks
       in the line search is nexpand */
    Parm->nexpand = (int) 50 ;

    /* maximum number of secant iterations in line search is nsecant */
    Parm->nsecant = (int) 50 ;

    /* conjugate gradient method restarts after (n*restart_fac) iterations */
    Parm->restart_fac = ONE ;

    /* stop when -alpha*dphi0 (estimated change in function value) <= feps*|f|*/
    Parm->feps = ZERO ;

    /* after encountering nan, growth factor when searching for
       a bracketing interval */
    Parm->nan_rho = 1.3 ;

    /* Wolfe line search parameter, range [0, .5]
       phi (a) - phi (0) <= delta phi'(0) */
    Parm->delta = .1 ;

    /* Wolfe line search parameter, range [delta, 1]
       phi' (a) >= sigma phi' (0) */
    Parm->sigma = .9 ;

    /* decay factor for bracket interval width in line search, range (0, 1) */
    Parm->gamma = .66 ;

    /* growth factor in search for initial bracket interval */
    Parm->rho = 5. ;

    /* conjugate gradient parameter beta_k must be >= eta*||d_k||_2 */
    Parm->eta = .01 ;

    /* starting guess for line search =
         psi0 ||x_0||_infty over ||g_0||_infty if x_0 != 0
         psi0 |f(x_0)|/||g_0||_2               otherwise */
    Parm->psi0 = .01 ;      /* factor used in starting guess for iteration 1 */

    /* for a QuadStep, function evaluated at psi1*previous step */
    Parm->psi1 = .1 ;

    /* when starting a new cg iteration, our initial guess for the line
       search stepsize is psi2*previous step */
    Parm->psi2 = 2. ;
}

/* =========================================================================
   === asa_descent =========================================================
   =========================================================================
   cg_descent conjugate gradient algorithm with modifications to handle the
   bound constraints.
   ========================================================================= */
int asa_descent /*  return:
                      -5 (ginorm < tau2*pgnorm without hitting boundary)
                      -4 (ginorm >=tau2*pgnorm, many x components hit boundary)
                      -3 (ginorm >=tau2*pgnorm, one x component hits boundary)
                      -2 (function value became nan)
                      -1 (starting function value is nan)
                       0 (convergence tolerance satisfied)
                       1 (change in func <= feps*|f|)
                       2 (total iterations exceeded maxit)
                       3 (slope always negative in line search)
                       4 (number secant iterations exceed nsecant)
                       5 (search direction not a descent direction)
                       6 (line search fails in initial interval)
                       7 (line search fails during bisection)
                       8 (line search fails during interval update)
                       9 (debugger is on and the function value increases)*/
(
    asa_com *Com
)
{
    int     i, iter, j, maxit, n, n5, nfree, nf, ng, nrestart,
            status ;
    double  delta2, eta_sq, Qk, Ck, pgnorm, ginorm,
            f, ftemp, gnorm, xnorm, gnorm2, dnorm2, denom,
            t, t1, t2, t3, t4, t5, dphi, dphi0, alpha, talpha,
            xj, gj, xg, xp, sts, sty, sk,
            yk, ykyk, ykgk, dkyk, yk1, yk2, yk3, yk4, yk5, beta,
            *x, *d, *g, *xtemp, *gtemp, *lo, *hi, *pg ;

    asacg_parm *Parm ;
    asa_parm *asaParm ;

/* initialization */

    x = Com->x ;
    lo = Com->lo ;
    hi = Com->hi ;
    n = Com->n ;
    d = Com->d ;
    g = Com->g ;
    xtemp = Com->xtemp ;
    gtemp = Com->gtemp ;
    pg = Com->pg ;
    nfree = Com->nfree ;
    nf = Com->nf ;
    ng = Com->ng ;
    pgnorm = Com->pgnorm ;
    ginorm = Com->ginorm ;
    Parm = Com->cgParm ;
    asaParm = Com->asaParm ;

    if ( Parm->PrintLevel >= 1 )
    {
        printf ("Dimension in CG, nfree = %i\n", nfree) ;
    }

    /* the conjugate gradient algorithm is restarted every nrestart iteration */
    nrestart = (INT) (((double) nfree)*Parm->restart_fac) ;

    /* abort when number of iterations reaches maxit in one pass through cg */
    if ( Parm->maxit_fac == INF ) Com->cgmaxit = maxit = INT_INF ;
    else Com->cgmaxit = maxit = (INT) (((double) n)*Parm->maxit_fac) ;

    n5 = nfree % 5 ;
    f = Com->f ;

    Ck = ZERO ;
    Qk = ZERO ;

/* initial function and gradient evaluations, initial direction */

    Com->f0 = f + f ;
    xnorm = asa_max (x, nfree) ;
    gnorm = ZERO ;
    gnorm2 = ZERO ;
    for (i = 0; i < n5; i++)
    {
        t = g [i] ;
        d [i] = -t ;
        gnorm2 += t*t ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
    }
    for (; i < nfree;)
    {
        t1 = g [i] ;
        d [i] = -t1 ;
        if ( gnorm < fabs (t1) ) gnorm = fabs (t1) ;
        i++ ;

        t2 = g [i] ;
        d [i] = -t2 ;
        if ( gnorm < fabs (t2) ) gnorm = fabs (t2) ;
        i++ ;

        t3 = g [i] ;
        d [i] = -t3 ;
        if ( gnorm < fabs (t3) ) gnorm = fabs (t3) ;
        i++ ;

        t4 = g [i] ;
        d [i] = -t4 ;
        if ( gnorm < fabs (t4) ) gnorm = fabs (t4) ;
        i++ ;

        t5 = g [i] ;
        d [i] = -t5 ;
        if ( gnorm < fabs (t5) ) gnorm = fabs (t5) ;
        i++ ;

        gnorm2 += t1*t1 + t2*t2 + t3*t3 + t4*t4 + t5*t5 ;
    }
    /* check that starting function value is nan */
    if ( f != f )
    {
        status = -1 ;
        goto Exit ;
    }

    if ( Parm->PrintLevel >= 2 )
    {
        printf ("iter: %5i f = %14.6e pgnorm = %14.6e ginorm = %14.6e\n\n",
          (int) 0, f, pgnorm, ginorm) ;
    }

    dphi0 = -gnorm2 ;
    delta2 = 2*Parm->delta - ONE ;
    eta_sq = Parm->eta*Parm->eta ;
    alpha = Parm->step ;
    if ( alpha == ZERO )
    {
        alpha = Parm->psi0*xnorm/gnorm ;
        if ( xnorm == ZERO )
        {
            if ( f != ZERO ) alpha = Parm->psi0*fabs (f)/gnorm2 ;
            else             alpha = ONE ;
        }
    }

/*  start the conjugate gradient iteration
    alpha starts as old step, ends as final step for current iteration
    f is function value for alpha = 0
    Com->QuadOK = TRUE means that a quadratic step was taken */

    for (iter = 1; iter <= maxit; iter++)
    {
        Com->QuadOK = FALSE ;
        alpha = Parm->psi2*alpha ;
        asa_maxstep (x, d, Com) ;
        if ( Parm->QuadStep )
        {
            if ( f != ZERO ) t = fabs ((f-Com->f0)/f) ;
            else             t = ONE ;
            if ( t > Parm->QuadCutOff )       /* take provisional step talpha */
            {
                talpha = Parm->psi1*alpha ;
                if ( Com->minstep >= asaParm->parm1*talpha )
                {
                    talpha = MIN (talpha, asaParm->parm2*Com->minstep) ;
                    asa_step (xtemp, x, d, talpha, nfree) ;
                    /*provisional function value*/
                    ftemp = asa_f (xtemp, Com) ;

                    /* check if function value is nan */
                    if ( ftemp != ftemp ) /* reduce stepsize */
                    {
                        for (i = 0; i < Parm->nexpand; i++)
                        {
                            talpha /= Parm->rho ;
                            asa_step (xtemp, x, d, talpha, nfree) ;
                            ftemp = asa_f (xtemp, Com) ;
                            if ( ftemp == ftemp ) break ;
                        }
                        if ( i == Parm->nexpand )
                        {
                            status = -2 ;
                            goto Exit ;
                        }
                    }

                    if ( ftemp < f )              /* check if QuadStep > 0 */
                    {
                       denom = TWO*(((ftemp-f)/talpha)-dphi0) ;
                       if ( denom > ZERO )    /* try a quadratic fit step */
                       {
                           Com->QuadOK = TRUE ;
                           alpha = -dphi0*talpha/denom ;
                       }
                    }
                }
            }
        }
        Com->f0 = f ;                          /* f0 saved as prior value */
        if ( Parm->PrintLevel >= 3 )
        {
            printf ("minstep =%14.6e, maxstep =%14.6e \n",
                     Com->minstep, Com->maxstep) ;
            printf ("QuadOK: %2i initial a: %14.6e f0: %14.6e dphi0: %14.6e\n",
                    Com->QuadOK, alpha, Com->f0, dphi0) ;
            if ( (alpha > Com->minstep) && Com->QuadOK )
            {
                printf("Quadratic step > minstep to boundary\n") ;
            }
        }

        /* parameters in Wolfe, approximate Wolfe conditions, and in update */
        Qk = Parm->Qdecay*Qk + ONE ;
        Ck = Ck + (fabs (f) - Ck)/Qk ;        /* average cost magnitude */

        if ( Parm->PertRule ) Com->fpert = f + Parm->eps*Ck ;
        else                  Com->fpert = f + Parm->eps ;

        Com->wolfe_hi = Parm->delta*dphi0 ;
        Com->wolfe_lo = Parm->sigma*dphi0 ;
        Com->awolfe_hi = delta2*dphi0 ;
        Com->alpha = alpha ;/* either double prior step or quadratic fit step */
        Com->f = f ;

        if ( Com->AWolfe )                  /* approximate Wolfe line search*/
        {
            if ( Parm->PrintLevel >= 3 )
            {
                printf ("Perform approximate Wolfe line search\n") ;
            }

            status = asa_line (dphi0, Com) ;
        }
        else                                  /* ordinary Wolfe line search */
        {
            if ( Parm->PrintLevel >= 3 )
            {
                 printf ("Perform ordinary Wolfe line search\n") ;
            }
            status = asa_lineW (dphi0, Com) ;
        }
        /* if ordinary Wolfe line search fails, possibly try approximate
           Wolfe line search*/
        if ( (status > 0) && !Com->AWolfe && (Parm->AWolfeFac > ZERO) )
        {
            Com->AWolfe = TRUE ;
            if ( Parm->PrintLevel >= 3 )
            {
                printf ("Ordinary Wolfe line search fails, "
                        "try approximate Wolfe line search\n") ;
            }

            status = asa_line (dphi0, Com) ;
        }

        alpha = Com->alpha ;
        f = Com->f ;
        dphi = Com->df ;

        if ( (status > 0) || (status == -1) || (status == -2) ) goto Exit ;

        /*Test for convergence to within machine epsilon
          [set feps to zero to remove this test] */

        if ( (-alpha*dphi0 <= Parm->feps*fabs (f)) && (status == 0) )
        {
            status = 1 ;
            goto Exit ;
        }

        /* compute beta, yk2, gnorm, gnorm2, dnorm2, update x and g */
        if ( iter % nrestart != 0 )
        {
            ginorm = ZERO ;
            pgnorm = ZERO ;
            for (j = 0; j < nfree; j++)
            {
                xj = xtemp [j] ;
                gj = gtemp [j] ;
                xg = xj - gj ;
                if      ( xg > hi [j] ) xp = hi [j] - xj ;
                else if ( xg < lo [j] ) xp = xj - lo [j] ;
                else                    xp = fabs (gj) ;
                pgnorm = MAX (pgnorm, xp) ;
                ginorm = MAX (ginorm, fabs (gj)) ;
            }
            for (; j < n; j++)
            {
                xj = x [j] ;
                gj = gtemp [j] ;
                xg = xj - gj ;
                if      ( xg > hi [j] ) xp = hi [j] - xj ;
                else if ( xg < lo [j] ) xp = xj - lo [j] ;
                else                    xp = fabs (gj) ;
                pgnorm = MAX (pgnorm, xp) ;
            }
            if ( asa_tol (pgnorm, Com) )
            {
                status = 0 ;
                for (j = 0; j < nfree; j++)
                {
                    xj = xtemp [j] ;
                    x [j] = xj ;
                    gj = gtemp [j] ;
                    xg = xj - gj ;
                    g [j] = gj ;
                    if      ( xg > hi [j] ) pg [j] = hi [j] - xj ;
                    else if ( xg < lo [j] ) pg [j] = lo [j] - xj ;
                    else                    pg [j] = -gj ;
                }
                for (; j < n; j++)
                {
                    xj = x [j] ;
                    gj = gtemp [j] ;
                    xg = xj - gj ;
                    g [j] = gj ;
                    if      ( xg > hi [j] ) pg [j] = hi [j] - xj ;
                    else if ( xg < lo [j] ) pg [j] = lo [j] - xj ;
                    else                    pg [j] = -gj ;
                }
                goto Exit1 ;
            }
            if ( ginorm < pgnorm*Com->tau2 ) status = -5 ;
            if ( status < -2 )
            {
                sts = ZERO ;
                sty = ZERO ;
                for (j = 0; j < nfree; j++)
                {
                    t = xtemp[j] ;
                    sk = t - x [j] ;
                    x [j] = t ;
                    sts += sk*sk ;

                    t = gtemp [j] ;
                    sty += sk*(t-g [j]) ;
                    g [j] = t ;
                }
                Com->sts = sts ;
                Com->sty = sty ;
                goto Exit ;
            }

            asa_copy (x, xtemp, nfree) ;
            dnorm2 = ZERO ;
            for (j = 0; j < n5; j++) dnorm2 = dnorm2 + d [j]*d [j] ;
            for (; j < nfree; j += 5)
            {
                dnorm2 = dnorm2 + d [j]*d [j] + d [j+1]*d [j+1]
                                              + d [j+2]*d [j+2]
                                              + d [j+3]*d [j+3]
                                              + d [j+4]*d [j+4] ;
            }
            ykyk = ZERO ;
            ykgk = ZERO ;
            for (j = 0; j < n5; j++)
            {
                t = gtemp [j] ;
                yk = t - g [j] ;
                ykyk += yk*yk ;
                ykgk += yk*t ;
                g [j] = t ;
            }
            for (j = n5; j < nfree; )
            {
                t1 = gtemp [j] ;
                yk1 = t1 - g [j] ;
                g [j] = t1 ;
                j++ ;

                t2 = gtemp [j] ;
                yk2 = t2 - g [j] ;
                g [j] = t2 ;
                j++ ;

                t3 = gtemp [j] ;
                yk3 = t3 - g [j] ;
                g [j] = t3 ;
                j++ ;

                t4 = gtemp [j] ;
                yk4 = t4 - g [j] ;
                g [j] = t4 ;
                j++ ;

                t5 = gtemp [j] ;
                yk5 = t5 - g [j] ;
                g [j] = t5 ;
                j++ ;

                ykyk += yk1*yk1 + yk2*yk2 + yk3*yk3 + yk4*yk4 + yk5*yk5 ;
                ykgk += yk1*t1 + yk2*t2 + yk3*t3 + yk4*t4 + yk5*t5 ;
            }

            dkyk = dphi - dphi0 ;
            beta = (ykgk - TWO*dphi*ykyk/dkyk)/dkyk ;
/*
    faster: initialize dnorm2 = gnorm2 at start, then
            dnorm2 = gnorm2 + beta**2*dnorm2 - 2.*beta*dphi
            gnorm2 = ||g_{k+1}||^2
            dnorm2 = ||d_{k+1}||^2
            dpi = g_{k+1}' d_k */

            t = -ONE/sqrt (dnorm2*MIN (eta_sq, gnorm2)) ;
            beta = MAX (beta, t) ;

/*    update search direction d = -g + beta*dold */

            gnorm2 = ZERO ;
            for (i = 0; i < n5; i++)
            {
                t = g [i] ;
                d [i] = -t + beta*d [i] ;
                gnorm2 += t*t ;
            }
            for (; i < nfree; )
            {
                t1 = g [i] ;
                d [i] = -t1 + beta*d [i] ;
                i++ ;

                t2 = g [i] ;
                d [i] = -t2 + beta*d [i] ;
                i++ ;

                t3 = g [i] ;
                d [i] = -t3 + beta*d [i] ;
                i++ ;

                t4 = g [i] ;
                d [i] = -t4 + beta*d [i] ;
                i++ ;

                t5 = g [i] ;
                d [i] = -t5 + beta*d [i] ;
                i++ ;

                gnorm2 += t1*t1 + t2*t2 + t3*t3 + t4*t4 + t5*t5 ;
            }
            dphi0 = -gnorm2 + beta*dphi ;
            if ( Parm->debug ) /* Check the dphi0 = d'g */
            {
                t = ZERO ;
                for (j=0; j<nfree; j++)  t = t + d[j]*g[j] ;
                if ( fabs(t-dphi0) > Parm->debugtol*fabs(dphi0) )
                {
                    printf("Warning, dphi0 != d'g!\n");
                    printf("dphi0:%14.6e, d'g:%14.6e\n",dphi0, t) ;
                }
            }
        }
        else
        {
            /* search direction d = -g */
            if ( Parm->PrintLevel >= 3 ) printf ("RESTART CG\n") ;
            ginorm = ZERO ;
            pgnorm = ZERO ;
            gnorm2 = ZERO ;
            for (j = 0; j < nfree; j++)
            {
                xj = xtemp [j] ;
                gj = gtemp [j] ;
                d [j] = -gj ;
                ginorm = MAX (ginorm, fabs (gj)) ;
                gnorm2 += gj*gj ;
                xg = xj - gj ;
                if      ( xg > hi [j] ) xp = hi [j] - xj ;
                else if ( xg < lo [j] ) xp = xj - lo [j] ;
                else                    xp = fabs (gj) ;
                pgnorm = MAX (pgnorm, xp) ;
            }
            for (; j < n; j++)
            {
                xj = x [j] ;
                gj = gtemp [j] ;
                xg = xj - gj ;
                if      ( xg > hi [j] ) xp = hi [j] - xj ;
                else if ( xg < lo [j] ) xp = xj - lo [j] ;
                else                    xp = fabs (gj) ;
                pgnorm = MAX (pgnorm, xp) ;
            }
            if ( asa_tol (pgnorm, Com) )
            {
                status = 0 ;
                for (j = 0; j < nfree; j++)
                {
                    xj = xtemp [j] ;
                    x [j] = xj ;
                    gj = gtemp [j] ;
                    xg = xj - gj ;
                    g [j] = gj ;
                    if      ( xg > hi [j] ) pg [j] = hi [j] - xj ;
                    else if ( xg < lo [j] ) pg [j] = lo [j] - xj ;
                    else                    pg [j] = -gj ;
                }
                for (; j < n; j++)
                {
                    xj = x [j] ;
                    gj = gtemp [j] ;
                    xg = xj - gj ;
                    g [j] = gj ;
                    if      ( xg > hi [j] ) pg [j] = hi [j] - xj ;
                    else if ( xg < lo [j] ) pg [j] = lo [j] - xj ;
                    else                    pg [j] = -gj ;
                }
                goto Exit1 ;
            }
            if ( ginorm < pgnorm*Com->tau2 ) status = -5 ;
            if ( status < -2 )
            {
                sts = ZERO ;
                sty = ZERO ;
                for (j = 0; j < nfree; j++)
                {
                    t = xtemp[j] ;
                    sk = t - x [j] ;
                    x [j] = t ;
                    sts += sk*sk ;

                    t = gtemp [j] ;
                    sty += sk*(t-g [j]) ;
                    g [j] = t ;
                }
                Com->sts = sts ;
                Com->sty = sty ;
                goto Exit ;
            }

            dphi0 = -gnorm2 ;
            asa_copy (x, xtemp, nfree) ;
            asa_copy (g, gtemp, nfree) ;
        }
        if ( !Com->AWolfe )
        {
            if ( fabs (f-Com->f0) <= Parm->AWolfeFac*Ck ) Com->AWolfe = TRUE ;
        }

        if ( Parm->PrintLevel >= 2 )
        {
            printf ("iter: %5i f = %14.6e pgnorm = %14.6e ginorm = %14.6e\n\n",
              (int) iter, f, pgnorm, ginorm) ;
        }

        if ( Parm->debug )
        {
            if ( f > Com->f0 + Ck*Parm->debugtol )
            {
                status = 9 ;
                goto Exit ;
            }
        }

        if ( dphi0 > ZERO )
        {
           status = 5 ;
           goto Exit ;
        }
    }
    status = 2 ;

Exit:
    if ( status < -2 )
    {
        for (j = nfree; j < n; j++) g [j] = gtemp [j] ;
    }
    else
    {
        pgnorm = ZERO ;
        for (j = 0; j < n; j++)
        {
            xj = xtemp [j] ;
            x [j] = xj ;
            gj = gtemp [j] ;
            g [j] = gj ;
            xg = xj - gj ;
            if      ( xg > hi [j] ) xp = hi [j] - xj ;
            else if ( xg < lo [j] ) xp = lo [j] - xj ;
            else                    xp = -gj ;
            pgnorm = MAX (pgnorm, fabs (xp)) ;
            pg [j] = xp ;
        }
    }

Exit1:
    Com->pgnorm = pgnorm ;
    Com->ginorm = ginorm ;
    Com->f = f ;
    Com->f_debug = f ;
    Com->cgfunc += Com->nf - nf ;
    Com->cggrad += Com->ng - ng ;
    Com->cgiter += iter ;
    if ( Parm->PrintLevel >= 2 )
    {
        printf ("iter: %5i f = %14.6e pgnorm = %14.6e ginorm = %14.6e\n\n",
                (int) iter, f, pgnorm, ginorm) ;
    }
    if ( Parm->PrintLevel >= 1 )
    {
        printf ("\nCG Termination status: %i\n", status) ;
        if ( status == -5 )
        {
            printf ("ginorm < tau2*pgnorm without hitting boundary\n") ;
        }
        if ( status == -4 )
        {
            printf ("ginorm >= tau2*pgnorm, many x components hit boundary\n") ;
        }
        else if ( status == -3 )
        {
            printf ("ginorm >= tau2*pgnorm, one x component hits boundary\n") ;
        }
        printf ("proj gradient max norm: %13.6e\n", pgnorm) ;
        printf ("function value:         %13.6e\n", f) ;
        printf ("cg iterations:          %13.6e\n", (double) iter) ;
        printf ("function evaluations:   %13.6e\n", (double) Com->nf - nf) ;
        printf ("gradient evaluations:   %13.6e\n", (double) Com->ng - ng) ;
    }
    return (status) ;
}

/* =========================================================================
   === asa_Wolfe ===========================================================
   =========================================================================
   Check whether the Wolfe or the approximate Wolfe conditions are satisfied
   ========================================================================= */
int asa_Wolfe
(
    double       alpha , /* stepsize */
    double           f , /* function value associated with stepsize alpha */
    double        dphi , /* derivative value associated with stepsize alpha */
    asa_com        *Com  /* cg com */
)
{
    if ( dphi >= Com->wolfe_lo )
    {

        /* test original Wolfe conditions */
        if ( f - Com->f0 <= alpha*Com->wolfe_hi )
        {
            if ( Com->cgParm->PrintLevel >= 4 )
            {
                printf ("wolfe f: %14.6e f0: %14.6e dphi: %14.6e\n",
                         f, Com->f0, dphi) ;
            }
            return (1) ;
        }
        /* test approximate Wolfe conditions */
        else if ( Com->AWolfe )
        {
            if ( (f <= Com->fpert) && (dphi <= Com->awolfe_hi) )
            {
                if ( Com->cgParm->PrintLevel >= 4 )
                {
                    printf ("f: %14.6e fpert: %14.6e dphi: %14.6e awolf_hi: "
                            "%14.6e\n", f, Com->fpert, dphi, Com->awolfe_hi) ;
                }
                return (1) ;
            }
        }
    }
    return (0) ;
}

/* =========================================================================
   === asa_tol =============================================================
   =========================================================================
   Check for convergence
   ========================================================================= */
int asa_tol
(
    double      pgnorm, /* projected gradient sup-norm */
    asa_com       *Com
)
{
    /*StopRule = T => |grad|_infty <=max (tol, |grad|_infty*StopFac)
                 F => |grad|_infty <= tol*(1+|f|)) */
    if ( Com->asaParm->StopRule )
    {
        if ( pgnorm <= Com->tol ) return (1) ;
    }
    else if ( pgnorm <= Com->tol*(ONE + fabs (Com->f)) ) return (1) ;
    return (0) ;
}

/* =========================================================================
   === asa_step ============================================================
   =========================================================================
   Compute xtemp = x + alpha d
   ========================================================================= */
void asa_step
(
    double *xtemp , /*output vector */
    double     *x , /* initial vector */
    double     *d , /* search direction */
    double  alpha , /* stepsize */
    INT         n   /* length of the vectors */
)
{
    INT n5, i ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) xtemp [i] = x[i] + alpha*d[i] ;
    for (; i < n; i += 5)
    {
        xtemp [i]   = x [i]   + alpha*d [i] ;
        xtemp [i+1] = x [i+1] + alpha*d [i+1] ;
        xtemp [i+2] = x [i+2] + alpha*d [i+2] ;
        xtemp [i+3] = x [i+3] + alpha*d [i+3] ;
        xtemp [i+4] = x [i+4] + alpha*d [i+4] ;
    }
}

/* =========================================================================
   === asa_line ============================================================
   =========================================================================
   Approximate Wolfe line search routine
   ========================================================================= */
int asa_line
(
    double       dphi0, /* function derivative at starting point (alpha = 0) */
    asa_com       *Com  /* cg com structure */
)
{
    int i, iter, nfree, nsecant, nshrink, ngrow, status ;
    double a, dphia, b, dphib, c, alpha, phi, dphi, alphaold, phiold,
           a0, da0, b0, db0, width, fquad, rho, minstep, maxstep,
           *x, *xtemp, *d, *gtemp ;
    asacg_parm *Parm ;

    nfree = Com->nfree ;
    x = Com->x ;         /* current iterate */
    d = Com->d ;         /* current search direction */
    xtemp = Com->xtemp ; /* x + alpha*d */
    gtemp = Com->gtemp ; /* gradient at x + alpha*d */
    minstep = Com->minstep ;
    maxstep = Com->maxstep ;
    alpha = Com->alpha ;
    if ( alpha > minstep )
    {
        alpha = minstep ;
        Com->QuadOK = FALSE ;
    }
    phi = Com->f ;
    Parm = Com->cgParm ;
    rho = Parm->rho ;
    asa_step (xtemp, x, d, alpha, nfree) ;
    asa_g (gtemp, xtemp, Com) ;
    dphi = asa_dot (gtemp, d, nfree) ;

    /* check if gradient is nan; if so, reduce stepsize */
    if ( dphi != dphi )
    {
        for (i = 0; i < Parm->nexpand; i++)
        {
            alpha /= rho ;
            asa_step (xtemp, x, d, alpha, nfree) ;
            asa_g (gtemp, xtemp, Com) ;
            dphi = asa_dot (gtemp, d, nfree) ;
            if ( dphi == dphi ) break ;
        }
        if ( i == Parm->nexpand )
        {
            status = -2 ;
            goto Exit ;
        }
        Com->QuadOK = FALSE ;
        rho = Parm->nan_rho ;
    }

/*Find initial interval [a,b] such that dphia < 0, dphib >= 0,
         and phia <= phi0 + feps*Ck */

    a = ZERO ;
    dphia = dphi0  ;
    ngrow = 0 ;
    nshrink = 0 ;
    while ( dphi < ZERO )
    {
        phi = asa_f (xtemp, Com) ;

/* if QuadStep in effect and quadratic conditions hold, check wolfe condition*/

        if ( Com->QuadOK )
        {
            if ( ngrow == 0 ) fquad = MIN (phi, Com->f0) ;
            if ( phi <= fquad )
            {
                if ( Parm->PrintLevel >= 4 )
                {
                    printf ("alpha: %14.6e phi: %14.6e fquad: %14.6e\n",
                            alpha, phi, fquad) ;
                }
                if ( asa_Wolfe (alpha, phi, dphi, Com) )
                {
                    status = 0 ;
                    goto Exit ;
                }
            }
        }
        if ( phi <= Com->fpert )
        {
            a = alpha ;
            dphia = dphi ;
        }
        else
        {
            /* contraction phase, only break at termination or Secant step */
            b = alpha ;
            while ( TRUE )
            {
                alpha = .5*(a+b) ;
                nshrink++ ;
                if ( nshrink > Parm->nexpand )
                {
                    status = 6 ;
                    goto Exit ;
                }
                asa_step (xtemp, x, d, alpha, nfree) ;
                asa_g (gtemp, xtemp, Com) ;
                dphi = asa_dot (gtemp, d, nfree) ;
                if ( dphi >= ZERO ) goto Secant ;
                phi = asa_f (xtemp, Com) ;
                if ( Parm->PrintLevel >= 4 )
                {
                    printf ("contract, a: %14.6e b: %14.6e alpha: %14.6e phi: "
                            "%14.6e dphi: %14.6e\n", a, b, alpha, phi, dphi) ;
                }
                if ( Com->QuadOK && (phi <= fquad) )
                {
                    if ( asa_Wolfe (alpha, phi, dphi, Com) )
                    {
                        status = 0 ;
                        goto Exit ;
                    }
                }
                if ( phi <= Com->fpert )
                {
                    a = alpha ;
                    dphia = dphi ;
                }
                else
                {
                    b = alpha ;
                }
            }
        }

/* expansion phase */

        ngrow++ ;
        if ( ngrow > Parm->nexpand )
        {
            status = 3 ;
            goto Exit ;
        }
        alphaold = alpha ;
        alpha = MIN (rho*alpha, minstep) ;
        if ( alpha != alphaold )
        {
            asa_step (xtemp, x, d, alpha, nfree) ;
            asa_g (gtemp, xtemp, Com) ;
            dphi = asa_dot (gtemp, d, nfree) ;
            if ( Parm->PrintLevel >= 4 )
            {
                printf ("expand,   a: %14.6e alpha: %14.6e phi: "
                         "%14.6e dphi: %14.6e\n", a, alpha, phi, dphi) ;
            }
        }
        else /* a new constraint is active */
        {
            do /* while statement */
            {
                alphaold = alpha ;
                phiold = phi ;
                if ( alpha < maxstep )
                {
                    alpha = rho*alphaold ;
                    asa_project (xtemp, x, d, alpha, Com) ;
                    phi = asa_f (xtemp, Com) ;
                }
            } while ( phi < phiold ) ;
            if ( alphaold == minstep )
            {
                asa_step (xtemp, x, d, minstep, nfree) ;
                status = -3 ;
            }
            else
            {
                asa_project (xtemp, x, d, alphaold, Com) ;
                asa_g (gtemp, xtemp, Com) ;
                status = -4 ;
            }
            phi = phiold ;
            goto Exit ;
        }
    }

Secant:
    b = alpha ;
    dphib = dphi ;
    if ( Com->QuadOK )
    {
        phi = asa_f (xtemp, Com) ;
        if ( ngrow + nshrink == 0 ) fquad = MIN (phi, Com->f0) ;
        if ( phi <= fquad )
        {
            if ( asa_Wolfe (alpha, phi, dphi, Com) )
            {
                status = 0 ;
                goto Exit ;
            }
        }
    }
    nsecant = Parm->nsecant ;
    for (iter = 1; iter <= nsecant; iter++)
    {
        if ( Parm->PrintLevel >= 4 )
        {
            printf ("secant, a: %14.6e b: %14.6e da: %14.6e db: %14.6e\n",
                     a, b, dphia, dphib) ;
        }
        width = Parm->gamma*(b - a) ;
        if ( -dphia <= dphib ) alpha = a - (a-b)*(dphia/(dphia-dphib)) ;
        else                   alpha = b - (a-b)*(dphib/(dphia-dphib)) ;
        c = alpha ;
        a0 = a ;
        b0 = b ;
        da0 = dphia ;
        db0 = dphib ;
        status = asa_update (&a, &dphia, &b, &dphib, &alpha, &phi,
                    &dphi, Com) ;
        if ( status >= 0 ) goto Exit ;
        else if ( status == -2 )
        {
            if ( c == a )
            {
                if ( dphi > da0 ) alpha = c - (c-a0)*(dphi/(dphi-da0)) ;
                else              alpha = a ;
            }
            else
            {
                if ( dphi < db0 ) alpha = c - (c-b0)*(dphi/(dphi-db0)) ;
                else              alpha = b ;
            }
            if ( (alpha > a) && (alpha < b) )
            {
                if ( Parm->PrintLevel >= 4 ) printf ("2nd secant\n") ;
                status = asa_update (&a, &dphia, &b, &dphib, &alpha, &phi,
                          &dphi, Com) ;
                if ( status >= 0 ) goto Exit ;
            }
        }

/* bisection iteration */

        if ( b-a >= width )
        {
            alpha = .5*(b+a) ;
            if ( Parm->PrintLevel >= 4 ) printf ("bisection\n") ;
            status = asa_update (&a, &dphia, &b, &dphib, &alpha, &phi,
                        &dphi, Com) ;
            if ( status >= 0 ) goto Exit ;
        }
        else if ( b <= a )
        {
            status = 7 ;
            goto Exit ;
        }
    }
    status = 4 ;

Exit:
    Com->alpha = alpha ;
    Com->f = phi ;
    Com->df = dphi ;
    return (status) ;
}

/* =========================================================================
   === asa_lineW ===========================================================
   =========================================================================
   Ordinary Wolfe line search routine.
   This routine is identical to asa_line except that the function
   psi [a] = phi [a] - phi [0] - a*delta*dphi [0] is minimized instead of
   the function phi
   ========================================================================= */
int asa_lineW
(
    double       dphi0 , /* function derivative at starting point (alpha = 0) */
    asa_com       *Com   /* cg com structure */
)
{
    int i, iter, nfree, nsecant, nshrink, ngrow, status ;
    double a, dpsia, b, dpsib, c, alpha, phi, dphi, alphaold, phiold,
           a0, da0, b0, db0, width, fquad, rho, psi, dpsi, minstep, maxstep,
           *x, *d, *xtemp, *gtemp ;
    asacg_parm *Parm ;

    nfree = Com->nfree ;
    x = Com->x ;         /* current iterate */
    d = Com->d ;         /* current search direction */
    xtemp = Com->xtemp ; /* x + alpha*d */
    gtemp = Com->gtemp ; /* gradient at x + alpha*d */
    minstep = Com->minstep ;
    maxstep = Com->maxstep ;
    alpha = Com->alpha ;
    if ( alpha > minstep )
    {
        alpha = minstep ;
        Com->QuadOK = FALSE ;
    }
    phi = Com->f ;
    Parm = Com->cgParm ;
    rho = Parm->rho ;
    asa_step (xtemp, x, d, alpha, nfree) ;
    asa_g (gtemp, xtemp, Com) ;
    dphi = asa_dot (gtemp, d, nfree) ;

    /* check if gradient is nan; if so, reduce stepsize */
    if ( dphi != dphi )
    {
        for (i = 0; i < Parm->nexpand; i++)
        {
            alpha /= rho ;
            asa_step (xtemp, x, d, alpha, nfree) ;
            asa_g (gtemp, xtemp, Com) ;
            dphi = asa_dot (gtemp, d, nfree) ;
            if ( dphi == dphi ) break ;
        }
        if ( i == Parm->nexpand )
        {
            status = -2 ;
            goto Exit ;
        }
        Com->QuadOK = FALSE ;
        rho = Parm->nan_rho ;
    }
    dpsi = dphi - Com->wolfe_hi ;

    /*Find initial interval [a,b] such that dphia < 0, dphib >= 0,
         and phia <= phi0 + feps*Ck */

    a = ZERO ;
    dpsia = dphi0 - Com->wolfe_hi ;
    ngrow = 0 ;
    nshrink = 0 ;
    while ( dpsi < ZERO )
    {
        phi = asa_f (xtemp, Com) ;
        psi = phi - alpha*Com->wolfe_hi ;

        /* if QuadStep in effect and quadratic conditions hold,
           check Wolfe condition*/

        if ( Com->QuadOK )
        {
            if ( ngrow == 0 ) fquad = MIN (phi, Com->f0) ;
            if ( phi <= fquad )
            {
                if ( Parm->PrintLevel >= 4 )
                {
                    printf ("alpha: %14.6e phi: %14.6e fquad: %14.6e\n",
                            alpha, phi, fquad) ;
                }
                if ( asa_Wolfe (alpha, phi, dphi, Com) )
                {
                    status = 0 ;
                    goto Exit ;
                }
            }
        }
        if ( psi <= Com->fpert )
        {
            a = alpha ;
            dpsia = dphi ;
        }
        else
        {
            /* contraction phase, only break at termination or Secant step */
            b = alpha ;
            while ( TRUE )
            {
                alpha = .5*(a+b) ;
                nshrink++ ;
                if ( nshrink > Parm->nexpand )
                {
                    status = 6 ;
                    goto Exit ;
                }
                asa_step (xtemp, x, d, alpha, nfree) ;
                asa_g (gtemp, xtemp, Com) ;
                dphi = asa_dot (gtemp, d, nfree) ;
                dpsi = dphi - Com->wolfe_hi ;
                if ( dpsi >= ZERO ) goto Secant ;
                phi = asa_f (xtemp, Com) ;
                psi = phi - alpha*Com->wolfe_hi ;
                if ( Parm->PrintLevel >= 4 )
                {
                    printf ("contract, a: %14.6e b: %14.6e alpha: %14.6e phi: "
                            "%14.6e dphi: %14.6e\n", a, b, alpha, phi, dphi) ;
                }
                if ( Com->QuadOK && (phi <= fquad) )
                {
                    if ( asa_Wolfe (alpha, phi, dphi, Com) )
                    {
                        status = 0 ;
                        goto Exit ;
                    }
                }
                if ( psi <= Com->fpert )
                {
                    a = alpha ;
                    dpsia = dpsi ;
                }
                else
                {
                    b = alpha ;
                }
            }
        }

/* expansion phase */

        ngrow++ ;
        if ( ngrow > Parm->nexpand )
        {
            status = 3 ;
            goto Exit ;
        }
        alphaold = alpha ;
        alpha = MIN (rho*alpha, minstep) ;
        if ( alpha != alphaold )
        {
            asa_step (xtemp, x, d, alpha, nfree) ;
            asa_g (gtemp, xtemp, Com) ;
            dphi = asa_dot (gtemp, d, nfree) ;
            dpsi = dphi - Com->wolfe_hi ;
            if ( Parm->PrintLevel >= 4 )
            {
                printf ("expand,   a: %14.6e alpha: %14.6e phi: "
                         "%14.6e dphi: %14.6e\n", a, alpha, phi, dphi) ;
            }
        }
        else /* a new constraint is active */
        {
            do /* while statement */
            {
                alphaold = alpha ;
                phiold = phi ;
                if ( alpha < maxstep )
                {
                    alpha = rho*alphaold ;
                    asa_project (xtemp, x, d, alpha, Com) ;
                    phi = asa_f (xtemp, Com) ;
                }
            } while ( phi < phiold ) ;
            if ( alphaold == minstep )
            {
                asa_step (xtemp, x, d, minstep, nfree) ;
                status = -3 ;
            }
            else
            {
                asa_project (xtemp, x, d, alphaold, Com) ;
                asa_g (gtemp, xtemp, Com) ;
                status = -4 ;
            }
            phi = phiold ;
            goto Exit ;
        }
    }

Secant:
    b = alpha ;
    dpsib = dpsi ;
    if ( Com->QuadOK )
    {
        phi = asa_f (xtemp, Com) ;
        if ( ngrow + nshrink == 0 ) fquad = MIN (phi, Com->f0) ;
        if ( phi <= fquad )
        {
            if ( asa_Wolfe (alpha, phi, dphi, Com) )
            {
                status = 0 ;
                goto Exit ;
            }
        }
    }
    nsecant = Parm->nsecant ;
    for (iter = 1; iter <= nsecant; iter++)
    {
        if ( Parm->PrintLevel >= 4 )
        {
            printf ("secant, a: %14.6e b: %14.6e da: %14.6e db: %14.6e\n",
                     a, b, dpsia, dpsib) ;
        }
        width = Parm->gamma*(b - a) ;
        if ( -dpsia <= dpsib ) alpha = a - (a-b)*(dpsia/(dpsia-dpsib)) ;
        else                   alpha = b - (a-b)*(dpsib/(dpsia-dpsib)) ;
        c = alpha ;
        a0 = a ;
        b0 = b ;
        da0 = dpsia ;
        db0 = dpsib ;
        status = asa_updateW (&a, &dpsia, &b, &dpsib, &alpha, &phi, &dphi,
                   &dpsi, Com) ;
        if ( status >= 0 ) goto Exit ;
        else if ( status == -2 )
        {
            if ( c == a )
            {
                if ( dpsi > da0 ) alpha = c - (c-a0)*(dpsi/(dpsi-da0)) ;
                else              alpha = a ;
            }
            else
            {
                if ( dpsi < db0 ) alpha = c - (c-b0)*(dpsi/(dpsi-db0)) ;
                else              alpha = b ;
            }
            if ( (alpha > a) && (alpha < b) )
            {
                if ( Parm->PrintLevel >= 4 ) printf ("2nd secant\n") ;
                status = asa_updateW (&a, &dpsia, &b, &dpsib, &alpha, &phi,
                   &dphi, &dpsi, Com) ;
                if ( status >= 0 ) goto Exit ;
            }
        }

/* bisection iteration */

        if ( b-a >= width )
        {
            alpha = .5*(b+a) ;
            if ( Parm->PrintLevel >= 4 ) printf ("bisection\n") ;
            status = asa_updateW (&a, &dpsia, &b, &dpsib, &alpha, &phi, &dphi,
                       &dpsi, Com) ;
            if ( status >= 0 ) goto Exit ;
        }
        else if ( b <= a )
        {
            status = 7 ;
            goto Exit ;
        }
    }
    status = 4 ;

Exit:
    Com->alpha = alpha ;
    Com->f = phi ;
    Com->df = dphi ;
    return (status) ;
}

/* =========================================================================
   === asa_update ==========================================================
   =========================================================================
   update returns: 8 if too many iterations
                   0 if Wolfe condition is satisfied
                  -1 if interval is updated and a search is done
                  -2 if the interval updated successfully
   ========================================================================= */
int asa_update
(
    double          *a , /* left side of bracketing interval */
    double      *dphia , /* derivative at a */
    double          *b , /* right side of bracketing interval */
    double      *dphib , /* derivative at b */
    double      *alpha , /* trial step (between a and b) */
    double        *phi , /* function value at alpha (returned) */
    double       *dphi , /* function derivative at alpha (returned) */
    asa_com       *Com   /* cg com structure */
)
{
    int nfree, nshrink, status ;
    double *x, *d, *xtemp, *gtemp ;
    asacg_parm *Parm ;

    nfree = Com->nfree ;
    x = Com->x ;         /* current iterate */
    d = Com->d ;         /* current search direction */
    xtemp = Com->xtemp ; /* x + alpha*d */
    gtemp = Com->gtemp ; /* gradient at x + alpha*d */
    Parm = Com->cgParm ;
    asa_step (xtemp, x, d, *alpha, nfree) ;
    *phi = asa_fg (gtemp, xtemp, Com) ;
    *dphi = asa_dot (gtemp, d, nfree) ;
    if ( Parm->PrintLevel >= 4 )
    {
        printf ("update alpha: %14.6e phi: %14.6e dphi: %14.6e\n",
                 *alpha, *phi, *dphi) ;
    }
    if ( asa_Wolfe (*alpha, *phi, *dphi, Com) )
    {
        status = 0 ;
        goto Exit2 ;
    }
    status = -2 ;
    if ( *dphi >= ZERO )
    {
        *b = *alpha ;
        *dphib = *dphi ;
        goto Exit2 ;
    }
    else
    {
        if ( *phi <= Com->fpert )
        {
            *a = *alpha ;
            *dphia = *dphi ;
            goto Exit2 ;
        }
    }
    nshrink = 0 ;
    *b = *alpha ;
    while ( TRUE )
    {
        *alpha = .5*(*a + *b) ;
        nshrink++ ;
        if ( nshrink > Parm->nexpand )
        {
            status = 8 ;
            goto Exit2 ;
        }
        asa_step (xtemp, x, d, *alpha, nfree) ;
        *phi = asa_fg (gtemp, xtemp, Com) ;
        *dphi = asa_dot (gtemp, d, nfree) ;
        if ( Parm->PrintLevel >= 4 )
        {
            printf ("contract, a: %14.6e alpha: %14.6e "
                    "phi: %14.6e dphi: %14.6e\n", *a, *alpha, *phi, *dphi) ;
        }
        if ( asa_Wolfe (*alpha, *phi, *dphi, Com) )
        {
            status = 0 ;
            goto Exit2 ;
        }
        if ( *dphi >= ZERO )
        {
            *b = *alpha ;
            *dphib = *dphi ;
            goto Exit1 ;
        }
        if ( *phi <= Com->fpert )
        {
            if ( Parm->PrintLevel >= 4 )
            {
                printf ("update a: %14.6e dphia: %14.6e\n", *alpha, *dphi) ;
            }
            *a = *alpha ;
            *dphia = *dphi ;
        }
        else *b = *alpha ;
    }
Exit1:
    status = -1 ;
Exit2:
    if ( Parm->PrintLevel >= 3 )
    {
        printf ("UP a: %14.6e b: %14.6e da: %14.6e db: %14.6e status: %i\n",
                 *a, *b, *dphia, *dphib, status) ;
    }
    return (status) ;
}

/* =========================================================================
   === asa_updateW =========================================================
   =========================================================================
   This routine is identical to asa_update except that the function
   psi [a] = phi [a] - phi [0] - a*delta*dphi [0] is minimized instead of
   the function phi. The return int has the following meaning:
                   8 if too many iterations
                   0 if Wolfe condition is satisfied
                  -1 if interval is updated and a search is done
                  -2 if the interval updated successfully
   ========================================================================= */
int asa_updateW
(
    double          *a , /* left side of bracketing interval */
    double      *dpsia , /* derivative at a */
    double          *b , /* right side of bracketing interval */
    double      *dpsib , /* derivative at b */
    double      *alpha , /* trial step (between a and b) */
    double        *phi , /* function value at alpha (returned) */
    double       *dphi , /* derivative of phi at alpha (returned) */
    double       *dpsi , /* derivative of psi at alpha (returned) */
    asa_com       *Com   /* cg com structure */
)
{
    double psi ;
    int nfree, nshrink, status ;
    double *x, *d, *xtemp, *gtemp ;
    asacg_parm *Parm ;

    nfree = Com->nfree ;
    x = Com->x ;         /* current iterate */
    d = Com->d ;         /* current search direction */
    xtemp = Com->xtemp ; /* x + alpha*d */
    gtemp = Com->gtemp ; /* gradient at x + alpha*d */
    Parm = Com->cgParm ;
    asa_step (xtemp, x, d, *alpha, nfree) ;
    *phi = asa_fg (gtemp, xtemp, Com) ;
    psi = *phi - *alpha*Com->wolfe_hi ;
    *dphi = asa_dot (gtemp, d, nfree) ;
    *dpsi = *dphi - Com->wolfe_hi ;
    if ( Parm->PrintLevel >= 4 )
    {
        printf ("update alpha: %14.6e psi: %14.6e dpsi: %14.6e\n",
                 *alpha, psi, *dpsi) ;
    }
    if ( asa_Wolfe (*alpha, *phi, *dphi, Com) )
    {
        status = 0 ;
        goto Exit2 ;
    }
    status = -2 ;
    if ( *dpsi >= ZERO )
    {
        *b = *alpha ;
        *dpsib = *dpsi ;
        goto Exit2 ;
    }
    else
    {
        if ( psi <= Com->fpert )
        {
            *a = *alpha ;
            *dpsia = *dpsi ;
            goto Exit2 ;
        }
    }
    nshrink = 0 ;
    *b = *alpha ;
    while ( TRUE )
    {
        *alpha = .5*(*a + *b) ;
        nshrink++ ;
        if ( nshrink > Parm->nexpand )
        {
            status = 8 ;
            goto Exit2 ;
        }
        asa_step (xtemp, x, d, *alpha, nfree) ;
        *phi = asa_fg (gtemp, xtemp, Com) ;
        *dphi = asa_dot (gtemp, d, nfree) ;
        *dpsi = *dphi - Com->wolfe_hi ;
        psi = *phi - *alpha*Com->wolfe_hi ;
        if ( Parm->PrintLevel >= 4 )
        {
            printf ("contract, a: %14.6e alpha: %14.6e "
                    "phi: %14.6e dphi: %14.6e\n", *a, *alpha, *phi, *dphi) ;
        }
        if ( asa_Wolfe (*alpha, *phi, *dphi, Com) )
        {
            status = 0 ;
            goto Exit2 ;
        }
        if ( *dpsi >= ZERO )
        {
            *b = *alpha ;
            *dpsib = *dpsi ;
            goto Exit1 ;
        }
        if ( psi <= Com->fpert )
        {
            if ( Parm->PrintLevel >= 4 )
            {
                printf ("update a: %14.6e dpsia: %14.6e\n", *alpha, *dpsi) ;
            }
            *a = *alpha ;
            *dpsia = *dpsi ;
        }
        else *b = *alpha ;
    }
Exit1:
    status = -1 ;
Exit2:
    if ( Parm->PrintLevel >= 3 )
    {
        printf ("UP a: %14.6e b: %14.6e da: %14.6e db: %14.6e status: %i\n",
                 *a, *b, *dpsia, *dpsib, status) ;
    }
    return (status) ;
}

/* =========================================================================
   === asa_project =========================================================
   =========================================================================
   Project a vector into the feasible set
   ========================================================================= */
void asa_project
(
    double  *xnew,
    double     *x,
    double     *d,
    double  alpha,
    asa_com  *Com   /* cg com structure */
)
{
    int j, n ;
    double t, *lo, *hi ;
    lo = Com->lo ;
    hi = Com->hi ;
    n = Com->nfree ;
    for (j = 0; j < n; j++)
    {
        t = x [j] + alpha*d [j] ;
        if      ( t > hi [j] ) t = hi [j] ;
        else if ( t < lo [j] ) t = lo [j] ;
        xnew [j] = t ;
    }
}

/* =========================================================================
   === asa_maxstep =========================================================
   =========================================================================
   Compute maximum step in the search direction until hitting boundary
   ========================================================================= */
void asa_maxstep
(
    double       *x, /* current iterate */
    double       *d, /* direction */
    asa_com    *Com
)
{
    double *lo, *hi, step, minstep, maxstep ;
    int j, n ;

    n = Com->nfree ;
    minstep = INF ;
    maxstep = ZERO ;
    lo = Com->lo ;
    hi = Com->hi ;

    for (j = 0;  j < n; j++)
    {
        if ( d [j] > ZERO )
        {
            if ( hi [j] < INF )
            {
                step = (hi [j] - x [j])/d [j] ;
                minstep = MIN (minstep, step) ;
                maxstep = MAX (maxstep, step) ;
            }
        }
        else if ( d [j] < ZERO )
        {
            if ( lo [j] >-INF )
            {
                step = (lo [j] - x [j])/ d [j] ;
                minstep = MIN (minstep, step) ;
                maxstep = MAX (maxstep, step) ;
            }
        }
    }
    Com->minstep = minstep ;
    Com->maxstep = maxstep ;
    return ;
}

/* =========================================================================
   === asa_grad_proj =======================================================
   =========================================================================
   Nonmonotone gradient projection algorithm modified to combine with cg_descent
   for bound constraints problem.

   Notation:
   fmin  = min {f(x_i),  0 <= i <= k}, k is current iteration (best value)
   fmax  = max {f_{k-i}, i = 0, 1, ... , min(k,m-1)}, m = memory
   fr    = reference function value
   f     = f(x_k), current function value
   fc    = maximum objective function value since the last minimum
          (best) function value was found. In other words, if
          k is the current iteration f(k_1) = fmin, then
          fc = max {f(x_i), k_1 <= i <= k}
   fcomp = min {fr, fmax}
   fr_pert = fr + pert, pert = Parm->eps*|fcomp| or Parm->eps
                        depending on PertRule
   fcomp_pert = fcomp + pert
   ftemp = f(xtemp), temporary (or a trial) function value at xtemp
   ========================================================================= */
int asa_grad_proj /*return:
                      -1 (give active constraints to cg routine)
                       0 (convergence tolerance satisfied)
                      11 (number of iterations or function evaluations
                          exceed limit)
                      12 (line search fails)
                      13 (search direction in linesearch is not descent)
                      14 (function value became nan) */
(
    asa_com *Com
)
{
    int count, i, ident, index, iter, j, ll, ll0, mcount, mm, n, nf, nf_line,
        ng, nl, np, status, hitbound, getbound, freebound ;
    double alpha, armijo_decay, armijo0, armijo1, f, fmin, fc, fr, sts, gtd,
           fmax, lambda, pgnorm, ginorm, gnorm, xnorm,
           xj, gj, xp, xg, t, th, tl,
           pert_lo, pert_hi, atemp, ftemp, fcomp, sty, yty, s, y, cosine,
           *lo, *hi, *x, *d, *g, *xtemp, *gtemp, *pg, *lastfvalues ;
    double fr_pert, fcomp_pert, dphia, Armijo_hi, AArmijo_hi ;
    asa_parm *Parm ;

    n = Com->n ;
    x = Com->x ;
    lo = Com->lo ;
    hi = Com->hi ;
    d = Com->d ;
    g = Com->g ;
    xtemp = Com->xtemp ;
    gtemp = Com->gtemp ;
    pg = Com->pg ;
    lastfvalues = Com->lastfvalues ;
    pgnorm = Com->pgnorm ;
    ginorm = Com->ginorm ;
    nf = Com->nf ;
    ng = Com->ng ;
    Parm = Com->asaParm ;
    pert_lo = Parm->pert_lo ;
    pert_hi = Parm->pert_hi ;
    armijo_decay = Parm->armijo_decay ;
    armijo0 = Parm->armijo0 ;
    armijo1 = Parm->armijo1 ;
    f = Com->f ;

    iter = 0 ;
    status = 0 ;
    count = 0 ;
    mcount = 2 ;
    ident = FALSE ;
    lambda = Com->alpha ;
    lastfvalues [0] = f ;
    for (i = 1; i < Parm->m; i++) lastfvalues [i] = -INF ;
    mm = 0 ;
    ll = 0 ;
    nl = 0 ;
    np = 0 ;
    fmin = f ;
    fr = f ;
    fc = f ;

    if ( Parm->PrintLevel >= 3 )
    {
        printf ("Initial stepsize in cbb: %14.6e\n", lambda) ;
    }

    while ( TRUE )
    {
        if ( Parm->PrintLevel >= 2 )
        {
            printf ("cbb iter: %5i f: %14.6e pgnorm: %14.6e\n\n",
                     iter, f, pgnorm) ;
        }

        if ( !Parm->GradProjOnly )
        {
            if ( ginorm >= Com->tau1*pgnorm )
            {
               if ( ident || (count >= mcount) )
               {
                   status = -1 ;
                   goto Exit ;
               }
            }
            else
            {
                if ( ident )
                {
                    ident = FALSE ;
                    Com->tau1 *= Parm->tau1_decay ;
                    Com->tau2 *= Parm->tau2_decay ;
                }
            }
        }
        iter++ ;
        hitbound = FALSE ;
        getbound = FALSE ;
        freebound = FALSE ;
        sts = ZERO ;
        gtd = ZERO ;
        for (j = 0; j < n; j++)
        {
            xj = x [j] ;
            gj = g [j] ;
            xp = -lambda*gj ;
            xg = xj + xp ;
            th = hi [j] ;
            tl = lo [j] ;
            if ( xg >= th )
            {
                xp = th - xj ;
                xtemp [j] = th ;
                if ( xp > pert_hi ) getbound = TRUE ;
            }
            else if ( xg <= tl )
            {
                xp = tl - xj ;
                xtemp [j] = tl ;
                if ( -xp > pert_lo ) getbound = TRUE ;
            }
            else
            {
                xtemp [j] = xg ;
                if ( (xj == th) || (xj == tl) ) freebound = TRUE ;
            }
            d [j] = xp ;
            gtd += gj*xp ; /* g'd (derivative in search direction) */
            sts += xp*xp ;
        }
        if ( getbound ) ll++ ;
        nf_line = Com->nf ;

        if (gtd >= ZERO)
        {
            status = 13 ;
            Com->gtd = gtd ;
            goto Exit ;
        }

       /* start of cbb line search */
        if ( Parm->PrintLevel >= 4 )
        {
            printf ("Linesearch in cbb, f: %14.6e gtd: %14.6e\n", f, gtd) ;
        }
        fmax = lastfvalues [0] ;
        for (i = 1; i < Parm->m; i++) fmax = MAX (fmax, lastfvalues [i]) ;
        alpha = ONE ;
        ftemp = asa_f (xtemp, Com) ;

        if ( nl == Parm->L )
        {
            fr = fmax ;
            t = (fr-fmin)/(fc-fmin) ;
            if ( t > Parm->gamma1 ) fr = fc ;
            nl = 0 ;
        }

        if ( (np > Parm->P) && (fmax > f) )
        {
           t = (fr-f)/(fmax-f) ;
           if ( t > Parm->gamma2 ) fr = fmax ;
        }

        fcomp = MIN (fmax, fr) ;

        if ( Parm->PrintLevel >= 4 )
        {
            printf ("fr: %14.6e fcomp: %14.6e\n", fr, fcomp) ;
        }

        /* Approximate nonmonotone Armijo line search, decrease alpha until:
           phi'(alpha) <= [2(phi_r - phi(0))/alpha] + (2 delta - 1) phi'(0) and
           phi(alpha) <= phi_r, where phi_r = fr_pert or fcomp_pert. */
        if ( Com->AArmijo)
        {
            if ( Parm->PertRule ) t = Parm->eps*fabs(fcomp) ;
            else                  t = Parm->eps ;
            fr_pert = fr + t ;
            fr = fr_pert ;
            fcomp_pert = fcomp + t ;
            if ( Parm->PrintLevel >= 3 )
            {
                printf ("Perform approximate Armijo line search\n") ;
                if ( Parm->PrintLevel >= 4 )
                {
                    printf ("fr_pert: %14.6e fcomp_pert: %14.6e\n",
                             fr_pert, fcomp_pert) ;
                }
            }

            AArmijo_hi = (TWO*Parm->delta - ONE)*gtd ;
            if ( ftemp != ftemp ) /* function value is nan, reduce stepsize */
            {
                for (i = 0; i < Parm->nshrink; i++)
                {
                    ll++ ;
                    alpha *= Parm->nan_fac ;
                    asa_step (xtemp, x, d, alpha, n) ;
                    ftemp = asa_f (xtemp, Com) ;
                    if ( ftemp == ftemp ) break ;
                }
                if ( (i == Parm->nshrink) || (alpha == ZERO) )
                {
                    status = 14 ;
                    goto exit_with_error ;
                }

                if ( ftemp <= fcomp_pert)
                {
                    asa_g (gtemp, xtemp, Com) ;
                    dphia = asa_dot (gtemp, d, n) ;
                    if (dphia <= TWO*(fcomp_pert - f)/alpha + AArmijo_hi )
                        goto exit_cbbls ; /* unit step is valid */
                }
            }
            else
            {
                if ( mm == 0 )
                {
                    if ( ftemp <= fr_pert )
                    {
                        asa_g (gtemp, xtemp, Com) ;
                        dphia = asa_dot (gtemp, d, n) ;
                        if (dphia <= TWO*(fr_pert - f) + AArmijo_hi )
                        {
                            mm++ ;
                            goto exit_cbbls ;
                        }
                    }
                }
                else
                {
                    if ( ftemp <= fcomp_pert )
                    {
                        asa_g (gtemp, xtemp, Com) ;
                        dphia = asa_dot (gtemp, d, n) ;
                        if (dphia <= TWO*(fcomp_pert - f) + AArmijo_hi )
                        {
                            mm++ ;
                            goto exit_cbbls ;
                        }
                    }
                }
            }

            /* backtracking approximate nonmonotone line search */
            ll0 = ll ;
            while ( TRUE )
            {
                /* Modified Raydan's quadratic interpolation line search */
                t = TWO*(ftemp-f-alpha*gtd) ;
                if ( t != ZERO )
                {
                    atemp = (-gtd*alpha*alpha)/t ;
                    if ( (atemp < armijo0*alpha) || (atemp > armijo1*alpha ) )
                    {
                        atemp = armijo_decay*alpha ;
                    }
                    alpha = atemp ;
                }
                else alpha *= armijo_decay ;

                asa_saxpy (xtemp, x, d, alpha, n) ; /* xtemp = x + alpha*d */
                ftemp = asa_f (xtemp, Com) ;
                ll++ ;

                if ( Parm->PrintLevel >= 4 )
                {
                    printf ("alpha: %14.6e ftemp: %14.6e\n", alpha, ftemp) ;
                }

                if ( ftemp <= fcomp_pert )
                {
                    asa_g (gtemp, xtemp, Com) ;
                    dphia = asa_dot (gtemp, d, n) ;
                    if (dphia <= TWO*(fcomp_pert - f)/alpha + AArmijo_hi )
                        goto exit_cbbls ;
                }

                if ( (alpha <= ZERO) || (ll-ll0 >= Parm->max_backsteps) )
                {
                   status = 12 ;
                   goto exit_with_error ;
                }
            }
            /* End of approximate Armijo line search */
        }

        /* Ordinary nonmonotone Armijo line search, decrease alpha until
           phi(alpha) <= phi_r + alpha * delta * phi'(0)
           where phi_r = fr or fcomp. */
        else
        {
            if ( Parm->PrintLevel >= 3 )
            {
                printf ("Perform ordinary Armijo line search\n") ;
            }

            Armijo_hi = Parm->delta*gtd ;
            if ( ftemp != ftemp ) /* function value is nan, reduce stepsize */
            {
                for (i = 0; i < Parm->nshrink; i++)
                {
                    ll++ ;
                    alpha *= Parm->nan_fac ;
                    asa_step (xtemp, x, d, alpha, n) ;
                    ftemp = asa_f (xtemp, Com) ;
                    if ( ftemp == ftemp ) break ;
                }
                if ( (i == Parm->nshrink) || (alpha == ZERO) )
                {
                    status = 14 ;
                    goto exit_with_error ;
                }
                if ( ftemp <= fcomp+alpha*Armijo_hi ) goto exit_cbbls ;
            }
            else
            {
                if ( mm == 0 ) t = fr ;
                else           t = fcomp ;
                if ( ftemp <= t+Armijo_hi )
                {
                    mm++ ;
                    goto exit_cbbls ;
                }
            }

            ll0 = ll ;
            while ( TRUE )
            {
                /* Modified Raydan's quadratic interpolation line search */
                t = TWO*(ftemp-f-alpha*gtd) ;
                if ( t != ZERO )
                {
                    atemp = (-gtd*alpha*alpha)/t ;
                    if ( (atemp < armijo0*alpha) || (atemp > armijo1*alpha ) )
                    {
                        atemp = armijo_decay*alpha ;
                    }
                    alpha = atemp ;
                }
                else alpha *= armijo_decay ;

                asa_saxpy (xtemp, x, d, alpha, n) ; /* xtemp = x + alpha*d */
                ftemp = asa_f (xtemp, Com) ;
                ll++ ;

                if ( Parm->PrintLevel >= 4 )
                {
                    printf ("alpha: %14.6e ftemp: %14.6e\n", alpha, ftemp) ;
                }

                if ( ftemp <= fcomp+alpha*Armijo_hi ) break ;

                if ( (alpha <= ZERO) || (ll-ll0 >= Parm->max_backsteps) )
                {
                    /* try approximate Armijo line search  */
                    if ( Parm->AArmijoFac > ZERO ) fr = fcomp ;
                    else                                 /* line search fails */
                    {
                        status = 12 ;
                        goto exit_with_error ;
                    }
                }
            }
            /* End of ordinary Armijo line search */
        }

        exit_cbbls:

        if ( ftemp <= fmin )
        {
             fmin = ftemp ;
             fc = ftemp ;
             nl = 0 ;
        }
        else nl++ ;
        if ( ftemp > fc ) fc = ftemp ;

        exit_with_error:
        /* end of cbbls */

        if ( getbound && (alpha == ONE) ) hitbound = TRUE ;
        if ( Parm->PrintLevel >= 3 )
        {
            printf ("hitbound = %i freebound = %i alpha = %14.6e\n",
                     hitbound, freebound, alpha) ;
        }

        if ( hitbound || freebound ) count = 0 ;
        else                         count++ ;

        sts *= alpha*alpha ;
        if ( Com->nf == nf_line + 1 ) np++ ;
        else                          np = 0 ;

        if ( !Com->AArmijo)  asa_g (gtemp, xtemp, Com) ;

        /* linesearch fails */
        if ( status > 0 )
        {
            if ( ftemp < f )
            {
                f = ftemp ;
                asa_copy(x, xtemp, n) ;
                asa_copy(g, gtemp, n) ;
            }
            pgnorm = ZERO ;
            for (j = 0; j < n; j++)
            {
                xj = x [j] ;
                gj = g [j] ;
                xg = xj - gj ;
                if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                else if ( xg <= lo [j] ) xp = lo [j] - xj ;
                else                     xp = -gj ;
                pgnorm = MAX (pgnorm, fabs (xp)) ;
                pg [j] = xp ;
            }
            goto Exit ;
        }

        index = 0 ;

        if ( (ll >= 1) || (mm >= Parm->nm) || (iter <= 1) )
        {
            index = 1 ;
            sty = ZERO ;
            pgnorm = ZERO ;
            ginorm = ZERO ;
            for (j = 0; j < n; j++)
            {
                xj = xtemp [j] ;
                gj = gtemp [j] ;
                xg = xj - gj ;
                if ( (xj - lo [j] > pert_lo) && (hi [j] - xj > pert_hi) )
                    ginorm = MAX (ginorm, fabs (gj));
                if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                else if ( xg <= lo [j] ) xp = xj - lo [j] ;
                else                     xp = fabs (gj) ;
                pgnorm = MAX (pgnorm, xp) ;
                sty += (xj - x [j])*(gj - g [j]) ;
                x [j] = xj ;
                g [j] = gj ;
            }

            if ( asa_tol (pgnorm, Com) )
            {
                f = ftemp ;
                for (j  = 0; j < n; j++)
                {
                    xj = xtemp [j] ;
                    gj = gtemp [j] ;
                    xg = xj - gj ;
                    if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                    else if ( xg <= lo [j] ) xp = lo [j] - xj ;
                    else                     xp = -gj ;
                    pg [j] = xp ;
                }
                status = 0 ;

                goto Exit ;
            }
        }

        else
        {
            pgnorm = ZERO ;
            ginorm = ZERO ;
            gnorm = ZERO ;
            sty = ZERO ;
            yty = ZERO ;
            for (j = 0; j < n; j++)
            {
                xj = xtemp [j] ;
                gj = gtemp [j] ;
                xg = xj - gj ;
                t = fabs (gj) ;
                gnorm = MAX (gnorm, t) ;
                if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                else if ( xg <= lo [j] ) xp = xj - lo [j] ;
                else                     xp = t ;
                pgnorm = MAX (pgnorm, xp) ;
                if ( (xj - lo [j] > pert_lo) && (hi [j] - xj > pert_hi) )
                {
                    ginorm = MAX (ginorm, t) ;
                }
                s = xj - x [j] ;
                y = gj - g [j] ;
                sty += s*y ;
                yty += y*y ;
                x [j] = xj ;
                g [j] = gj ;
            }
            if ( asa_tol (pgnorm, Com) )
            {
                f = ftemp ;
                for (j  = 0; j < n; j++)
                {
                    xj = xtemp [j] ;
                    gj = gtemp [j] ;
                    xg = xj - gj ;
                    if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                    else if ( xg <= lo [j] ) xp = lo [j] - xj ;
                    else                     xp = -gj ;
                    pg [j] = xp ;
                }
                status = 0 ;
                goto Exit ;
            }
            s = Parm->parm3*fabs (ftemp)/gnorm ;
            t = MAX (s, ONE) ;
            if ( sts > t*t ) index = 1 ;
            else
            {
                t = MIN (s, ONE) ;
                if ( sts <= t*t )
                {
                    cosine = fabs (sty)/sqrt (sts*yty) ;
                    if ( cosine >= Parm->gamma ) index = 1 ;
                }
            }
        }

        if ( index == 1 )
        {
            ll = 0 ;
            if ( sty <= ZERO)
            {
                if ( mm >= Parm->parm4 )
                {
                    xnorm = asa_max (x, n) ;
                    t = MIN (ONE/pgnorm, xnorm/pgnorm) ;
                    lambda = MAX (t, lambda) ;
                    mm = 0 ;
                }
            }
            else
            {
                t = MAX (Parm->lmin, sts/sty) ;
                lambda = MIN (Parm->lmax, t) ;
                mm = 0 ;
            }
        }

        /* If not GradProjOnly, check if the active constraints are identified*/
        if ( !Parm->GradProjOnly &&
              pgnorm < Parm->pgdecay*MAX (ONE, Com->pgnorm_start) )
        {
            ident = asa_identify(x, g, pgnorm, Com) ;
        }

        f = ftemp ;
        lastfvalues [iter % Parm->m] = f ;

        /* check for excessive iterations/function evaluations */
        if ( (iter >= Com->pgmaxit) || (Com->nf - nf  >= Com->pgmaxfunc) )
        {
            for (j = 0; j < n; j++)
            {
                xj = x [j] ;
                gj = g [j] ;
                xg = xj - gj ;
                if      ( xg >= hi [j] ) xp = hi [j] - xj ;
                else if ( xg <= lo [j] ) xp = lo [j] - xj ;
                else                     xp = -gj ;
                pg [j] = xp ;
            }
            status = 11 ;
            goto Exit ;
        }

        if ( !Com->AArmijo )
        {
            if ( fabs(fr - f) <= Parm->AArmijoFac*fabs(fcomp) )
                Com->AArmijo = TRUE ;
        }
    }
    Exit:
    Com->f = f ;
    Com->ginorm = ginorm ;
    Com->pgnorm = pgnorm ;
    Com->cbbiter += iter ;
    Com->cbbfunc += Com->nf - nf ;
    Com->cbbgrad += Com->ng - ng ;
    if ( Parm->PrintLevel >= 2 )
    {
        if(status != -1) printf ("cbb iter: %5i f: %14.6e pgnorm: %14.6e\n\n",
                                  iter, f, pgnorm) ;
    }
    if ( Parm->PrintLevel >= 1 )
    {
        printf ("\nCBB Termination status: %i\n", status) ;
        if ( status == -1 )
            printf ("terminate cbb iteration, branch to cg iteration\n") ;

        printf ("proj gradient max norm: %13.6e\n", pgnorm) ;
        printf ("function value:         %13.6e\n", f) ;
        printf ("cbb iterations:         %13.6e\n", (double) iter) ;
        printf ("function evaluations:   %13.6e\n", (double) Com->nf - nf) ;
        printf ("gradient evaluations:   %13.6e\n", (double) Com->ng - ng) ;
    }
    return (status) ;
}

/* =========================================================================
   === asa_init_bbstep =====================================================
   =========================================================================
   Calculate initial BB stepsize
   ========================================================================= */
double asa_init_bbstep
(
    asa_com *Com
)
{
    int n ;
    double alpha, lmax, lmin, pgnorm, xnorm, sts, sty, t, *x ;
    x = Com->x ;
    sts = Com->sts ;
    sty = Com->sty ;
    pgnorm = Com->pgnorm ;
    n = Com->n ;
    lmin = Com->asaParm->lmin ;
    lmax = Com->asaParm->lmax ;

    if ( sty > ZERO )
    {
        t = MIN (sts/sty, lmax) ;
        alpha = MAX (lmin, t) ;
    }
    else
    {
        xnorm = asa_max (x, n) ;
        if ( xnorm > ZERO ) alpha = MIN (ONE, xnorm)/pgnorm ;
        else                alpha = ONE/pgnorm ;
    }
    return (alpha) ;
}

/* =========================================================================
   ==== asa_f ==============================================================
   Evaluate the function
   =========================================================================*/
double asa_f
(
    double    *x,
    asa_com *Com
)
{
    double f ;
    asa_objective *user ;
    user = Com->user ;
    user->x = x ;
    Com->nf++ ;
    if ( Com->DimReduce )
    {
        /* Expand x to the full space*/
        asa_expandx (x, Com) ;

        /* Evaluate function */
        user->ifree = Com->ifree ;
        user->nfree = Com->nfree ;
        f = Com->value (user) ;

        /* Shrink x to the reduced space */
        asa_shrinkx (x, Com) ;
    }
    else
    {
        /* Evaluate function */
        user->ifree = NULL ;
        user->nfree = Com->n ;
        f = Com->value (user) ;
    }
    return (f) ;

}

/* =========================================================================
   ==== asa_g ==============================================================
   Evaluate the gradient
   =========================================================================*/
void asa_g
(
    double    *g,
    double    *x,
    asa_com *Com
)
{
    asa_objective *user ;
    user = Com->user ;
    user->x = x ;
    user->g = g ;
    Com->ng++ ;
    if ( Com->DimReduce )
    {
        /* Expand x to the full space*/
        asa_expandx (x, Com) ;

        /* Evaluate gradient */
        user->ifree = Com->ifree ;
        user->nfree = Com->nfree ;
        Com->grad (user) ;

        /* Shrink x and g to the reduced space */
        asa_shrinkxg (x, g, Com) ;
    }
    else
    {
        /* Evaluate gradient */
        user->ifree = NULL ;
        user->nfree = Com->n ;
        Com->grad (user) ;
    }
}


/* =========================================================================
   ==== asa_fg =============================================================
   Evaluate the function and gradient
   =========================================================================*/
double asa_fg
(
    double    *g,
    double    *x,
    asa_com *Com
)
{
    asa_objective *user ;
    double f ;
    Com->nf++ ;
    Com->ng++ ;
    user = Com->user ;
    user->x = x ;
    user->g = g ;
    if ( Com->DimReduce )
    {
        /* Expand x to the full space*/
        asa_expandx (x, Com) ;

        /* Evaluate function and gradient */
        user->ifree = Com->ifree ;
        user->nfree = Com->nfree ;
        if ( Com->valgrad != NULL )
        {
            f = Com->valgrad (user) ;
        }
        else
        {
            Com->grad (user) ;
            f = Com->value (user) ;
        }

        /* Shrink x and g to the reduced space */
        asa_shrinkxg (x, g, Com) ;
    }
    else
    {
        /* Evaluate function and gradient */
        user->ifree = NULL ;
        user->nfree = Com->n ;
        if ( Com->valgrad != NULL )
        {
            f = Com->valgrad (user) ;
        }
        else
        {
            Com->grad (user) ;
            f = Com->value (user) ;
        }
    }
    return (f) ;
}

/* =========================================================================
   ==== asa_identify =======================================================
   Check whether the bounds with strict complementarity
   are approximately identified
   =========================================================================*/
int asa_identify
(
   double     *x,
   double     *g,
   double pgnorm,
   asa_com  *Com
)
{
    int ident, j, n ;
    double t, t1, xj, *lo, *hi ;
    n = Com->n ;
    lo = Com->lo ;
    hi = Com->hi ;
    ident = TRUE ;
    t = sqrt (pgnorm) ;
    t1 = t*t*t ;
    for (j = 0; j < n; j++)
    {
        xj = x [j] ;
        if ( ((xj -lo [j] >= t1) || (hi [j] - xj >= t1))
             && (fabs (g [j]) > t) ) ident = FALSE ;
    }
    return (ident) ;
}

/* =========================================================================
   === asa_expandx =========================================================
   =========================================================================
   Expand x array from size nfree to full size of dimension n based on
   indices of free variables
   ========================================================================= */
void asa_expandx
(
    double    *x,
    asa_com *Com
)
{
    int i, j, nfree, *ifree ;
    double t ;
    ifree = Com->ifree ;
    nfree = Com->nfree ;
    for (j = nfree-1; j >= 0; j--)
    {
        i = ifree [j] ;
        if ( j != i )
        {
            t = x [i] ;
            x [i] = x [j] ;
            x [j] = t ;
        }
    }
}

/* =========================================================================
   === asa_shrinkx =========================================================
   =========================================================================
   Compress x array to dimension nfree based on indices of free variables
   ========================================================================= */
void asa_shrinkx
(
    double    *x,
    asa_com *Com
)
{
    int i, j, nfree, *ifree ;
    double t ;
    ifree = Com->ifree ;
    nfree = Com->nfree ;
    for (j = 0; j < nfree; j++)
    {
        i = ifree [j] ;
        if ( j != i )
        {
            t = x [i] ;
            x [i] = x [j] ;
            x [j] = t ;
        }
    }
}

/* =========================================================================
   === asa_shrinkxg ========================================================
   =========================================================================
   Compress x and g arrays based on indices of free variables
   ========================================================================= */
void asa_shrinkxg
(
    double    *x,
    double    *g,
    asa_com *Com
)
{
    int i, j, nfree, *ifree ;
    double t ;
    ifree = Com->ifree ;
    nfree = Com->nfree ;
    for (j = 0; j < nfree; j++)
    {
        i = ifree [j] ;
        if ( j != i )
        {
            t = x [i] ;
            x [i] = x [j] ;
            x [j] = t ;

            t = g [i] ;
            g [i] = g [j] ;
            g [j] = t ;
        }
    }
}

/* =========================================================================
   === asa_expand_all ======================================================
   =========================================================================
   Expand vectors x, g, pg, lo and hi from the reduced space (dimension nfree)
   to the full space (dimension n).
   ========================================================================= */
void asa_expand_all
(
    asa_com *Com
)
{
    int i, j, nfree, *ifree ;
    double t, *x, *g, *pg, *lo, *hi ;
    x = Com->x ;
    g = Com->g ;
    pg = Com->pg ;
    lo = Com->lo ;
    hi = Com->hi ;
    ifree = Com->ifree ;
    nfree = Com->nfree ;
    for (j = nfree-1; j >= 0; j--)
    {
        i = ifree [j] ;
        if ( j != i )
        {
            t = x [i] ;
            x [i] = x [j] ;
            x [j] = t ;

            t = g [i] ;
            g [i] = g [j] ;
            g [j] = t ;

            t = pg [i] ;
            pg [i] = pg [j] ;
            pg [j] = t ;

            t = lo [i] ;
            lo [i] = lo [j] ;
            lo [j] = t ;

            t = hi [i] ;
            hi [i] = hi [j] ;
            hi [j] = t ;
        }
    }
}

/* =========================================================================
   === asa_shrink_all ======================================================
   =========================================================================
   Shrink vectors x, g, lo and hi from the full space (dimension n)
   to the reduced space (dimension nfree).
   ========================================================================= */

void asa_shrink_all
(
    asa_com *Com
)
{
    int i, j, nfree, *ifree ;
    double t, *lo, *hi, *g, *x ;
    x = Com->x ;
    g = Com->g ;
    lo = Com->lo ;
    hi = Com->hi ;
    ifree = Com->ifree ;
    nfree = Com->nfree ;
    for (j = 0; j < nfree; j++)
    {
        i = ifree [j] ;
        if ( i != j )
        {
            t = x [i] ;
            x [i] = x [j] ;
            x [j] = t ;

            t = g [i] ;
            g [i] = g [j] ;
            g [j] = t ;

            t = lo [i] ;
            lo [i] = lo [j] ;
            lo [j] = t ;

            t = hi [i] ;
            hi [i] = hi [j] ;
            hi [j] = t ;
        }
    }
}

/* =========================================================================
   === asa_dot =============================================================
   =========================================================================
   Compute dot product of x and y, vectors of length n
   ========================================================================= */
double asa_dot
(
    double *x , /* first vector */
    double *y , /* second vector */
    int     n   /* length of vectors */
)
{
    int i, n5 ;
    double t ;
    t = ZERO ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) t += x [i]*y [i] ;
    for (; i < n; i += 5)
    {
        t += x [i]*y[i] + x [i+1]*y [i+1] + x [i+2]*y [i+2]
                        + x [i+3]*y [i+3] + x [i+4]*y [i+4] ;
    }
    return (t) ;
}

/* =========================================================================
   === asa_copy ============================================================
   =========================================================================
   Copy vector x into vector y
   ========================================================================= */
void asa_copy
(
    double *y, /* target vector */
    double *x, /* given vector */
    int     n  /* dimension */
)
{
    int j, n10 ;
    n10 = n % 10 ;
    for (j = 0; j < n10; j++) y [j] = x [j] ;
    for (; j < n; j += 10)
    {
        y [j] = x [j] ;
        y [j+1] = x [j+1] ;
        y [j+2] = x [j+2] ;
        y [j+3] = x [j+3] ;
        y [j+4] = x [j+4] ;
        y [j+5] = x [j+5] ;
        y [j+6] = x [j+6] ;
        y [j+7] = x [j+7] ;
        y [j+8] = x [j+8] ;
        y [j+9] = x [j+9] ;
    }
}


/* =========================================================================
   === asa_saxpy ===========================================================
   =========================================================================
   Compute z = y + ax
   ========================================================================= */
void asa_saxpy
(
    double *z,
    double *y,
    double *x,
    double  a,
    int     n
)
{
    int j, n5 ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) z [j] = y [j] + a*x [j] ;
    for (; j < n; j += 5)
    {
        z [j]   = y [j]   + a*x [j] ;
        z [j+1] = y [j+1] + a*x [j+1] ;
        z [j+2] = y [j+2] + a*x [j+2] ;
        z [j+3] = y [j+3] + a*x [j+3] ;
        z [j+4] = y [j+4] + a*x [j+4] ;
    }
}
/* =========================================================================
   === asa_max =============================================================
   =========================================================================
   Return max {fabs (x [j]) : 1 <= j < n}
   ========================================================================= */
double asa_max
(
    double *x,
    int     n
)
{
    double xnorm ;
    int j, n5 ;
    n5 = n % 5 ;
    xnorm = ZERO ;
    for (j = 0; j < n5; j++) if ( xnorm < fabs (x [j]) ) xnorm = fabs (x [j]) ;
    for (; j < n; j += 5)
    {
        if ( xnorm < fabs (x [j]  ) ) xnorm = fabs (x [j]) ;
        if ( xnorm < fabs (x [j+1]) ) xnorm = fabs (x [j+1]) ;
        if ( xnorm < fabs (x [j+2]) ) xnorm = fabs (x [j+2]) ;
        if ( xnorm < fabs (x [j+3]) ) xnorm = fabs (x [j+3]) ;
        if ( xnorm < fabs (x [j+4]) ) xnorm = fabs (x [j+4]) ;
    }
    return (xnorm) ;
}

/* =========================================================================
   === asa_printcgParms ====================================================
   =========================================================================
   Print the contents of the asacg_parm structure
   ========================================================================= */
void asa_printcgParms
(
    asacg_parm  *Parm
)
{
    printf ("\nCG PARAMETERS:\n") ;
    printf ("\n") ;
    printf ("Wolfe line search parameter ..................... delta: %e\n",
             Parm->delta) ;
    printf ("Wolfe line search parameter ..................... sigma: %e\n",
             Parm->sigma) ;
    printf ("decay factor for bracketing interval ............ gamma: %e\n",
             Parm->gamma) ;
    printf ("growth factor for bracket interval ................ rho: %e\n",
             Parm->rho) ;
    printf ("growth factor for bracket interval after nan .. nan_rho: %e\n",
             Parm->nan_rho) ;
    printf ("truncation factor for cg beta ..................... eta: %e\n",
             Parm->eta) ;
    printf ("perturbation parameter for function value ......... eps: %e\n",
             Parm->eps) ;
    printf ("factor for computing average cost .............. Qdecay: %e\n",
             Parm->Qdecay) ;
    printf ("relative change in cost to stop QuadStep ... QuadCutOff: %e\n",
             Parm->QuadCutOff) ;
    printf ("stop when cost change <= feps*|f| ................. eps: %e\n",
             Parm->feps) ;
    printf ("cost change factor, approx Wolfe transition . AWolfeFac: %e\n",
             Parm->AWolfeFac) ;
    printf ("restart cg every restart_fac*n iterations . restart_fac: %e\n",
             Parm->restart_fac) ;
    printf ("starting guess parameter in first iteration ...... psi0: %e\n",
             Parm->psi0) ;
    printf ("factor multiply starting guess in quad step ...... psi1: %e\n",
             Parm->psi1) ;
    printf ("initial guess factor for general iteration ....... psi2: %e\n",
             Parm->psi2) ;
    printf ("starting step in first iteration if nonzero ...... step: %e\n",
             Parm->step) ;
    printf ("max expansions in line search ................. nexpand: %i\n",
             Parm->nexpand) ;
    printf ("max secant iterations in line search .......... nsecant: %i\n",
             Parm->nsecant) ;
    printf ("max cg iterations is n*maxit_fac ............ maxit_fac: %e\n",
             Parm->maxit_fac) ;
    printf ("total max cg iterations is n*totit_fac ...... totit_fac: %e\n",
             Parm->totit_fac) ;
    printf ("error tolerance when debugger turned on ..... .debugtol: %e\n",
             Parm->debugtol) ;
    printf ("print level (0 = none, 4 = maximum) ........ PrintLevel: %i\n",
             Parm->PrintLevel) ;
    printf ("\nLogical parameters:\n") ;
    if ( Parm->PertRule )
        printf ("    Error estimate for function value is eps*Ck\n") ;
    else
        printf ("    Error estimate for function value is eps\n") ;
    if ( Parm->QuadStep )
        printf ("    Use quadratic interpolation step\n") ;
    else
        printf ("    No quadratic interpolation step\n") ;
    if ( Parm->PrintParms )
        printf ("    Print the parameter structure\n") ;
    else
        printf ("    Do not print parameter structure\n") ;
    if ( Parm->AWolfe)
        printf ("    Approximate Wolfe line search\n") ;
    else
        printf ("    Wolfe line search") ;
        if ( Parm->AWolfeFac > ZERO )
            printf (" ... switching to approximate Wolfe\n") ;
        else
            printf ("\n") ;
    if ( Parm->debug)
        printf ("    Check for decay of cost, debugger is on\n") ;
    else
        printf ("    Do not check for decay of cost, debugger is off\n") ;
}

/* =========================================================================
   === asa_printParms ====================================================
   =========================================================================
   Print the contents of the asa_parm structure
   ========================================================================= */
void asa_printParms
(
    asa_parm  *Parm
)
{
    printf ("\nASA PARAMETERS:\n") ;
    printf ("\n") ;
    printf ("update fr if fmin not improved after L iterations.... L: %i\n",
             Parm->L) ;
    printf ("fmax = max (f_{k-i}, i = 0, 1, ..., min (k, m-1) )... m: %i\n",
             Parm->m) ;
    printf ("update fr if P previous initial stepsizes accepted... P: %i\n",
             Parm->P) ;
    printf ("CBB cycle length.................................... nm: %i\n",
             Parm->nm) ;
    printf ("criterion for updating reference value fr....... gamma1: %e\n",
             Parm->gamma1) ;
    printf ("criterion for updating reference value fr....... gamma2: %e\n",
             Parm->gamma2) ;
    printf ("max tries to find non NAN function value ...... nshrink: %i\n",
             Parm->nshrink) ;
    printf ("interval decay factor in NAN search ............nan_fac: %e\n",
             Parm->nan_fac) ;
    printf ("perturbation parameter for function value.......... eps: %e\n",
             Parm->eps) ;
    printf ("cost change factor, approx Armijo transition,AArmijoFac: %e\n",
             Parm->AArmijoFac) ;
    printf ("Armijo line search parameter .................... delta: %e\n",
             Parm->delta) ;
    printf ("Armijo decay factor .......................armijo_decay: %e\n",
             Parm->armijo_decay) ;
    printf ("criterion for Q interpolation, cbb line search,.armijo0: %e\n",
             Parm->armijo0) ;
    printf ("criterion for Q interpolation, cbb line search,.armijo1: %e\n",
             Parm->armijo1) ;
    printf ("criterion for reinitializing BB stepsize ........ gamma: %e\n",
             Parm->gamma) ;
    printf ("Lower bound for initial stepsize ................. lmin: %e\n",
             Parm->lmin) ;
    printf ("Upper bound for initial stepsize ................. lmax: %e\n",
             Parm->lmax) ;
    printf ("used when trying a quadratic interpolation step.. parm1: %e\n",
             Parm->parm1) ;
    printf ("used when trying a quadratic interpolation step.. parm2: %e\n",
             Parm->parm2) ;
    printf ("criterion for reinitializing the BB stepsize..... parm3: %e\n",
             Parm->parm3) ;
    printf ("maximum previous BB steps used when s^t y <= 0... parm4: %i\n",
             Parm->parm4) ;
    printf ("if ginorm < tau1*pgnorm, continue grad_proj ...... tau1: %e\n",
             Parm->tau1) ;
    printf ("decay factor for tau1 ...................... tau1_decay: %e\n",
             Parm->tau1_decay) ;
    printf ("ginorm < tau2*pgnorm => subproblem solved in cg... tau2: %e\n",
             Parm->tau2) ;
    printf ("decay factor for tau2 ...................... tau2_decay: %e\n",
             Parm->tau2_decay) ;
    printf ("max number of Armijo backtracking steps . max_backsteps: %i\n",
             Parm->max_backsteps) ;
    printf ("max cbb iterations in 1 pass is n*maxit_fac . maxit_fac: %e\n",
             Parm->maxit_fac) ;
    printf ("total number cbb iterations is n*totit_fac .. totit_fac: %e\n",
             Parm->totit_fac) ;
    printf ("max func evals in cbb is n*maxfunc_fac .... maxfunc_fac: %e\n",
             Parm->maxfunc_fac) ;
    printf ("criterion for checking undecided index set..... pgdecay: %e\n",
             Parm->pgdecay) ;
    printf ("perturbation of lower bounds .................. pert_lo: %e\n",
             Parm->pert_lo) ;
    printf ("perturbation of upper bounds .................. pert_hi: %e\n",
             Parm->pert_hi) ;
    printf ("factor multiplying gradient in stop condition . StopFac: %e\n",
             Parm->StopFac) ;
    printf ("print level (0 = none, 4 = maximum) ........ PrintLevel: %i\n",
             Parm->PrintLevel) ;
    printf ("\nLogical parameters:\n") ;
    if ( Parm->PertRule )
        printf ("    Error estimate for function value is eps*|fcomp|\n") ;
    else
        printf ("    Error estimate for function value is eps\n") ;
    if ( Parm->PrintFinal )
        printf ("    Print final cost and statistics\n") ;
    else
        printf ("    Do not print final cost and statistics\n") ;
    if ( Parm->PrintParms )
        printf ("    Print the parameter structure\n") ;
    else
        printf ("    Do not print parameter structure\n") ;
    if ( Parm->AArmijo)
        printf ("    Approximate nonmonotone Armijo line search\n") ;
    else
        printf ("    Nonmonotone Armijo line search") ;
        if ( Parm->AArmijoFac > ZERO )
            printf (" ... switching to approx nonmonotone Armijo\n") ;
        else
            printf ("\n") ;
    if ( Parm->StopRule )
    {
        if ( Parm->StopFac == ZERO )
        {
            printf ("    Stopping condition based on gradient tolerance\n") ;
        }
        else
        {
            printf ("    Stopping condition uses initial grad tolerance\n") ;
        }
    }
    else
        printf ("    Stopping condition weighted by absolute cost\n") ;
    if ( Parm->GradProjOnly )
        printf ("    Only use the gradient projection algorithm\n") ;
    else
        printf ("    Apply gradient projection algorithm and cg_descent\n") ;
}
/*
Version 1.1 Change:

    1. Pass a structure asa_objective to the user evaluation routines.
       This allows asa_cg to pass more information to the user which
       might be used to speedup his routines to evaluate the objective
       function and its gradient.  Two elements of the structure are
       ifree and nfree.  If ifree is not NULL, then ifree is a pointer
       to an integer array containing the indices of the free variables
       while nfree is the number of free variables.

    2. Halt the Armijo backtracking line search in cbb when the number
       of backtracking steps reaching Parm->max_backsteps

Version 1.2 Change:

    1. Correct the Armijo line search in cbb by dividing by including the
       factor "/alpha" in the termination condition
*/




















double myvalue
(
    asa_objective *asa
) ;

void mygrad
(
    asa_objective *asa
) ;

double myvalgrad
(
    asa_objective *asa
) ;



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

void LocalOpt(ContPopul* epop,int nelit, int epsize, int atgen)
{
 double obj_f;
 double* TheBestsol = new double[vars];
 
 int i,k,j,eval_local,start;


   
 start = 0; 
 //cout<<"tract_vars is "<<tract_vars<<" length is "<<length<<endl;
 
for(k=start; k < epsize;  k++)  
 {
   ConvertFromCont(epop->S[k]);

   MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);	
   if (isnan(obj_f)) obj_f = 1000000.0; 
  //cout<<setprecision(10)<<endl<<"MGA_DSM objective value = " << obj_f<<endl;
      
  
   if(k<nelit) eval_local = ContRandom(vars, epop->S[k], 10
, &obj_f, TheBestsol);
   else eval_local = ContRandom(vars, epop->S[k], 50, &obj_f, TheBestsol);

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
 double* TheBestsol = new double[vars];
 
 int k,j,eval_local,start;
 //if (atgen==0) start=0;
 //else start=nelit;



  asacg_parm cgParm ;
  asa_parm asaParm ;
  asa_cg_default (&cgParm) ;
  asa_default (&asaParm) ;


     cgParm.PrintParms = TRUE ;
    //cgParm.PrintLevel = 0 ;
    //asaParm.PrintParms = TRUE ;
    //asaParm.PrintLevel = 0 ;

 
 start = 0; 

for(k=start; k < epsize;  k++)  
 {
   ConvertFromCont(epop->S[k]);

   MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);	
   if (isnan(obj_f)) obj_f = 1000000.0; 
 
   cout<<setprecision(10)<<endl<<"MGA_DSM objective value = " << obj_f<<endl;   

   asa_cg(epop->S[k], minbound, maxbound, vars, NULL, NULL, NULL, 1.e-8, myvalue, mygrad, NULL, NULL) ;
   
   ConvertFromCont(epop->S[k]);
   MGA_DSM(t, dim, sequence, c_o, rp, e, Isp, mass, AUdist, DVtotal, DVonboard, obj_f, Delta_V);	
   cout<<setprecision(10)<<endl<<" value = " << obj_f<<endl;   


   if (isnan(obj_f)) obj_f = 1000000.0; 
  
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

 
double myvalue /* evaluate the objective function */
(
    asa_objective *asa
)
{
  double f, t,objf ;
  double* x;
 
    INT i, n ;
    x = asa->x ;
    n = asa->n ;
     
     
    f = thefunct(x);

   return (f) ;
}

void mygrad /* evaluate the gradient of the objective function */
(
    asa_objective *asa
)
{
    double tt, *g, *x ;
    INT i, n ;
    x = asa->x ;
    g = asa->g ;
    n = asa->n ;
    for (i = 0; i < n; i++)
    {
        tt = i + 1 ;
        tt = sqrt (tt) ;
        //g[i] = exp(x[i]) -  tt ;
        g[i] = 0; //myrand();
    }
    g[randomint(n)]=1;
    return ;
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
   
  while (i<Maxgen && BestEval<Max && NPoints>Elit/5) 
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
  while (i<Maxgen && BestEval<Max && NPoints>=Elit) 
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


 /*
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
 */

 ContPopul* cpop =  new ContPopul(7000,vars,Elit,minbound,maxbound); 

 ReadSols(cpop);  

delete cpop;

         
//PrintStatistics();             
	delete[] AbsBestInd;   
     
      
 delete [] params; 
 delete[] Delta_V;
return 0;

}      




