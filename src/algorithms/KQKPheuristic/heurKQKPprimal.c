#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <memory.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <sys/times.h>
#include <unistd.h>
#include <float.h>

/* Bring in the declarations for the string functions */
#define TESTS       10   /* number of test to run in each series */
#define MSIZE      645   /* maximum number of 0-1 variables */
#define INFTY    FLT_MAX  /* very large float used in subgradient optimization */
#define NBSOL  200    /* Number of heuristic solutions */
#define PRIMAL_HEUR  1    /* Code for Primal Heuristic */
#define ROUNDING_HEUR  2    /* Code for Rounding Heuristic */

/* ======================================================================
				   macros
   ====================================================================== */

#define srand(x)     srand48(x)
#define randm(x)    (lrand48() % (long) (x))
#define MIN(a, b)                  ((a) < (b) ? (a) : (b))
#define MAX(a, b)                  ((a) > (b) ? (a) : (b))
#define ABS(a)                  ((a) < 0 ? -(a) : (a))
#define DET(a1, a2, b1, b2)     ((a1) * (ptype) (b2) - (a2) * (ptype) (b1))
#define SWAP(a, b)  { register lpitem t; t = *(a); *(a) = *(b); *(b) = t; }


/* ======================================================================
         type declarations
   ====================================================================== */

typedef int     boolean;
typedef int     ntype;   /* number of items          */
typedef int     itype;   /* item profits and weights */
typedef int     stype;   /* sum of profit or weight  */
typedef double  ptype;   /* product type             */
typedef float   etype;   /* efficiency type          */


/* ======================================================================
 			    global variables
   ====================================================================== */

/**
 * Name of the heuristic printed in the output file
 */
char heur_name[] = "Generic";

int n, c, nbobj;
float z, BORNESDP;
float p[MSIZE][MSIZE];
int w[MSIZE];
int x[MSIZE];
double xsdp[MSIZE];
boolean xstar[MSIZE];
boolean xsol[NBSOL][MSIZE]; /* Solutions heuristiques */
int nbsolheur; /* nb of heuristic solutions */
ntype nbobjets;
itype ERREUR;
FILE *trace;
FILE *matq;
FILE *matc;
FILE *matl; // todo: add l.txt
FILE *matA;
FILE *matb;
FILE *matAbis;
FILE *matbbis;
double heur_val = 0; // initialization value

void terminate_qkp(char *str);
void maketest(int n1, int r, int pct);
void checksolution(int c, float z);
int checksolutionunit(int c, float z);
void primalHeuristic(int *x, int N, double alpha);
void roundingHeuristic(int * x, int N, double alpha);


void load_data(int *n_) ; 


/* ======================================================================
				main
   ====================================================================== */

int main(int argc, char *argv[])
{
  int r, pct, v, i, j, fix, k, l;//, z_H_Pri;
  float z, z_H_Sur;
  double time;
  FILE *sol;
  double a=0.0;
  int no_nul;
  int fin=0;
  float zbest=0.0;
  int max_problem = 1; //1 if maximization, 0 if minimzation
  int heuristic_code = PRIMAL_HEUR; //PRIMAL_HEUR or ROUNDING_HEUR
//  int heuristic_code = ROUNDING_HEUR; //PRIMAL_HEUR or ROUNDING_HEUR

  load_data(&n) ;


  if (max_problem) zbest = -1e9; else zbest = 1e9;


    nbsolheur=0;
    for(l=0; l<NBSOL ; l++)
      for (i=0; i< n; i++)
        xsol[l][i]=0;

	switch (heuristic_code) {

	// called at the beginning
	case PRIMAL_HEUR:
    // printf("case 1 \n ") ; 

		for (i = 0; i < 100; i++) {

			primalHeuristic(xstar, n, ((double) i) / (double) 100);

			// copy the result only if it's better than the current and it's a feasible solution
    		memcpy(xsol[nbsolheur], xstar, MSIZE*sizeof(boolean));
//			if (checksolutionunit(c,zbest) && 
//				((max_problem && zbest < heur_val)
//					|| (!max_problem && zbest > heur_val))) {
			if (checksolutionunit(c,zbest) && zbest!= heur_val) {
		    nbsolheur++;
				zbest = heur_val;
			}

		}

		break;

	case ROUNDING_HEUR:

      // printf("case 2 \n "); 
      /*
    
      matq = fopen("q.txt", "w");
      matc = fopen("c.txt", "w");
      matA = fopen("A.txt", "w");
      matb = fopen("b.txt", "w");
      matAbis = fopen("Abis.txt", "w");
      matbbis = fopen("bbis.txt", "w");

      // Ecriture dans q.txt
      no_nul=(n*(n-1))/2;
      fprintf(matq,"%d %d\n",n,no_nul);
      for ( i=1;i<n;i++)
        for (j=i+1;j<n+1;j++){
          fprintf(matq,"%d %d ",i,j);
          if (i==n-1 && j==n)
            fprintf(matq,"%f ",-p[i-1][j-1]);
          else
            fprintf(matq,"%f\n",-p[i-1][j-1]);
        }

      // Ecriture dans c.txt 
      fprintf(matc,"%d\n",1);
      for (i=1;i<n;i++)
        fprintf(matc,"%f\n",-p[i-1][i-1]);
      fprintf(matc,"%f ",-p[n-1][n-1]);

      // Ecriture dans A.txt
      fprintf(matA,"%d\n",1);
      for (i=1;i<n;i++){
        fprintf(matA,"%d %d ",1,i);
        fprintf(matA,"%d\n",1);
      }
      fprintf(matA,"%d %d ",1,n);
      fprintf(matA,"%d ",1);

      // Ecriture dans b.txt
      fprintf(matb,"%d ",nbobj);

      // Ecriture dans Abis.txt
      fprintf(matAbis,"%d\n",1);
      for (i=1;i<n;i++){
        fprintf(matAbis,"%d %d ",1,i);
        fprintf(matAbis,"%d\n",w[i-1]);
      }
      fprintf(matAbis,"%d %d ",1,n);
      fprintf(matAbis,"%d ",w[n-1]);

      // Ecriture dans bbis.txt
      fprintf(matbbis,"%d ",c);

      fclose(matq);
      fclose(matA);
      fclose(matc);
      fclose(matb);
      fclose(matAbis);
      fclose(matbbis);

//       fflush(NULL);
//       printf("\n\n%d: c %d\n", v, c);
//       fflush(NULL);

//       cpu_time(NULL);

      system("./QCR_E-kQKP");
      system("grep Primal sol_chr.txt > valeur.txt");
//       system("grep Primal sol.txt > valeur.txt");
      FILE *tmp_file = fopen("valeur.txt", "r");
      char tmp[100];
      char *aie = fgets(tmp, 99, tmp_file);
      fclose(tmp_file);
      if(aie == NULL)
      {
        fflush(NULL);
        printf("Infaisable ou diverge\n");
        fflush(NULL);
//        break;
      }
      else
      {
        float bound;
        sscanf(tmp, "Primal objective value: %f", &bound);
        BORNESDP = bound;
	
//	        fflush(NULL);
//        printf("valeur bound: %f\n", bound);
//        fflush(NULL);
        sol = fopen("prob.sol","r");
        fin=0;
        while(fin < 2)
        {
          fscanf(sol,"%lf",&a);
          if (a==2)
          {
            fscanf(sol,"%lf",&a);
            if (a==1)
              fin++;
          }
        }
        fscanf(sol,"%d",&i);
        fscanf(sol,"%d",&j);
        fscanf(sol,"%lf",&a);
        for (l=0;l<=n;l++)
        {
          fscanf(sol,"%d",&i);
          fscanf(sol,"%d",&j);
          fscanf(sol,"%d",&i);
          fscanf(sol,"%d",&j);
          fscanf(sol,"%lf",&a);
          xsdp[l]=a;
        }
        fclose(sol);
      }

		for (i = 0; i < 100; i++) {
			roundingHeuristic(xstar, n, ((double) i) / (double) 100); // change the sign because minimization problem

    		memcpy(xsol[nbsolheur], xstar, MSIZE*sizeof(boolean));
//			if (checksolutionunit(c,zbest) && 
//				((max_problem && zbest < heur_val)
//					|| (!max_problem && zbest > heur_val))
			if (checksolutionunit(c,zbest) && zbest!= heur_val) {
		        nbsolheur++;
				zbest = heur_val;
			}
		}

    */
		break;

	  default:
		printf("Choosen heuristic doesn't exist\n");
		exit(1);
	}

  checksolution(c,zbest);

    // printf(" nbsolheur = %d \n", nbsolheur);

    FILE* f = fopen("heur_sol.data", "w") ; 
    fprintf(f, "%d \n", nbsolheur) ;

    for (i = 0; i < nbsolheur; i++)
    {
      for (j = 0 ; j < n; j++ )
        fprintf(f, "%d ", xsol[i][j]) ; 
      fprintf(f, "\n") ;
    }
    
    fclose(f) ; 


    system("rm -f *.txt") ;
    system("rm -f *.sb") ; 
    system("rm -f *.sdpa") ; 
    system("rm -f *.sol") ; 

  return EXIT_SUCCESS;
}  /* END main */



void load_data(int *n_ ){
  FILE *f = fopen("inst.data", "r") ; 
  int l;
  int i, j;
  // read n
  fscanf(f, "%d\n", n_) ;

  // read matrix p
  for(i=0; i< *n_; i++)
    for(j=0; j< *n_; j++){
      p[i][j] = 0.0 ;
    }
  
  fscanf(f, "%d\n", &l) ;

  for(int iter =0; iter<l; iter++){
    float val;
    fscanf(f, "%d %d %f\n", &i, &j, &val) ;
    p[i-1][j-1] = val ;
    p[j-1][i-1] = val; 
  }


  for(i = 0; i< *n_; i++){
    int val ;
    fscanf(f, "%d %d\n", &j, &val) ; 
    w[j-1] = val;
  }

  // read capa
  fscanf(f, "%d\n", &c) ;

  // read k-item
  fscanf(f, "%d\n", &nbobj) ;
  fclose(f) ;
}




void primalHeuristic(int *x, int N, double alpha) {
	int i;
	double random[N];
//	srand(time(NULL));

	for (i = 0; i < N; i++) {
		random[i] = ((double) rand() / (double) RAND_MAX);
	}

	for (i = 0; i < N; i++) {
		if (random[i] > alpha) {
			x[i] = 1;
		} else {
			x[i] = 0;
		}
	}
}

void roundingHeuristic(int *x,int N, double alpha)
//	int *x;int N; // Number of variables
//	double alpha; //
{
	int i;

	for (i = 0; i < N; i++) {
		if (xsdp[i] > alpha)
				x[i] = 1;
			else
				x[i] = 0;
	}
}

/* =======================================================================
                                  terminate_qkp
   ======================================================================= */

void terminate_qkp(char *str)
{
  printf("%s\n", str);
  printf("PROGRAM IS TERMINATED !!!\n\n");
  exit(-1);
}


/* ======================================================================
        maketest
   ====================================================================== */

void maketest(int n1, int r, int pct)
{
  int i, j;
  float psum;
  int wsum;

  n = n1;
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      p[i][j] = p[j][i] = (float) (randm(100) >= pct ? 0 : randm(r)+1);
    }
    w[i] = randm(r/2)+1;
//    printf("w[%d]=%d\n", i, w[i]);
  }
  psum = 0; wsum = 0;
  for (i = 0; i < n; i++) {
    x[i] = 0; wsum += w[i];
    for (j = 0; j < n; j++) psum += p[i][j];
  }
  if (wsum - 50 <= 0) terminate_qkp("too small weight sum");
  nbobj = randm(n/4) + 1;
  if (nbobj==1)
    nbobj++;
  nbobjets = nbobj;
//   printf("\n\n nbobj: %d\n", nbobj);
//   c    = randm(wsum-50) + 50;
//   c    = randm((nbobj*100)-50) + 50;
  c    = randm((nbobj*30)-50) + 50;

//  printf("MAKETEST %d: psum %f wsum %d c %d nbobjets %d nbobj %d\n", n, psum, wsum, c, nbobjets, nbobj);
}

/* ======================================================================
        checksolution
   ====================================================================== */

void checksolution(int c, float z)
{
  int i, j, l;
  float psum;
  int wsum, ksum;

  for(l=0;l<nbsolheur;l++)
  {
    psum = 0.0;
    wsum = ksum = 0;
    for (i = 0; i < n; i++) {
      if (!xsol[l][i]) continue;
//    printf("w[%d]=%d \n", i, w[i]);
      wsum += w[i];
      ksum++;
//    printf("i=%d x=%d\n", i, x[i]);
      for (j = 0; j < n; j++) {
        if (xsol[l][j])
        {
          psum += p[i][j];
//         printf("psum=%f et p[%d][%d]=%f, \n", psum, i, j, p[i][j]);
        }
      }
    }
//    printf("CHECK %d: nbsolheur %d psum %f ksum %d \n", n, l, psum, ksum);
    // printf("%f \n", psum);
    //if (ksum != nbobj) terminate_qkp("bad nb item");
    //if (psum != z) terminate_qkp("bad solution");
  }
  // printf("\n");
// printf("CHECK %d: psum %f z %f wsum %d c %d ksum %d nbobj %d\n", n, psum,z,wsum,c,ksum,nbobj);
// if (wsum > c) terminate_qkp("excess weight");
}

/* ======================================================================
        checksolutionunit
   ====================================================================== */

int checksolutionunit(int c, float z)
{
  int i, j, l;
  float psum;
  int wsum, ksum;

  psum = 0.0;
  wsum = ksum = 0;
  for (i = 0; i < n; i++) {
    if (!xsol[nbsolheur][i]) continue;
//    printf("w[%d]=%d \n", i, w[i]);
    wsum += w[i];
    ksum++;
//    printf("i=%d x=%d\n", i, x[i]);
    for (j = 0; j < n; j++) {
      if (xsol[nbsolheur][j])
      {
        psum += p[i][j];
//         printf("psum=%f et p[%d][%d]=%f, \n", psum, i, j, p[i][j]);
      }
    }
  }
  heur_val = psum;
//  printf("CHECKUNIT %d: nbsolheur %d psum %f ksum %d \n", n, nbsolheur, psum, ksum);
    //if (ksum != nbobj) terminate_qkp("bad nb item");
    //if (psum != z) terminate_qkp("bad solution");
//  printf("\n");
// printf("CHECK %d: psum %f z %f wsum %d c %d ksum %d nbobj %d\n", n, psum,z,wsum,c,ksum,nbobj);
// if (wsum > c) terminate_qkp("excess weight");
  if (ksum==nbobj && wsum <= c)
  {
//    printf("CHECKUNIT %d: nbsolheur %d psum %f ksum %d \n", n, nbsolheur, psum, ksum);
    return 1;
  }
  else
    return 0;
}

