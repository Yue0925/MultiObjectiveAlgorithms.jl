/*------------------------------------------------------------------------*/
/*  File: BestHybrid.c                                            */
/*  Version 1.0                                                           */
/*------------------------------------------------------------------------*/

/* ----- insta.data format
* n
* lines for matrix p
* i j p[i][j]
* ...
* i w[i]
* ...
* c
* nbobj
*/

/* ----- heur_sol.data
* val heuristic_value
* i 1
* ...
*/

/* BestHybrid.c - Entering and optimizing an EkQKP problem with real cofficients in the objective */
/*
 * The present code is written in ANSI-C, and has been compiled with
 * the GNU-C compiler under MAC OS X 14.6.1. To obtain an exacutable code
 * compile with
 *
 *   gcc -O3 BestHybrid_objreel_yue.c -o BestHybrid_objreel_yue
 *
 * The program reads three arguments
 *
 *   n     number of items 
 *   r     range of coefficients i.e. profits weights are in [1,r]
 *   pct   density of instance, i.e. expected frequency of nonzero
 *         profits
 *
 * having read the parameters, 10 instances are constructed and 
 * solved using the sdp_heur algorithm.
 *
 * (c) Copyright, January 2024,
 *
 * This code can be used free of charge for research and academic purposes
 * only.*/

/* This algorithm gives a lower bound for the exact k-item quadratic knapsack problem */

/* Bring in the declarations for the string functions */
#define TESTS       10   /* number of test to run in each series */
#define MSIZE      645   /* maximum number of 0-1 variables */
#define INFTY    FLT_MAX  /* very large float used in subgradient optimization */
#define EPSILON  1E-6    /* small tolerance used in subgradient optimization */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
//#include <values.h>
#include <limits.h>
#include <string.h>
#include <math.h>
//#include <malloc.h>
#include <sys/times.h>
#include <unistd.h>
#include <float.h>


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

/* The problem we are optimizing will have ? rows, ? columns and 
   ? nonzeros, and ? nonzeros in the quadratic coefficient matrix.  */

// #define NUMROWS    ?
// #define NUMCOLS    ?
// #define NUMNZ      ?
// #define NUMQNZ     ?

/* ======================================================================
         type declarations
   ====================================================================== */

typedef int     boolean;
typedef int     ntype;   /* number of items          */
typedef int     itype;   /* item profits and weights */
typedef int     stype;   /* sum of profit or weight  */
typedef double  ptype;   /* product type             */
typedef float   etype;   /* efficiency type          */

typedef struct {               /* item for solving lp-relaxation in row */
  etype e; /* profit to weight ratio */
  etype p; /* profit */
  etype w; /* weight */
  etype *x; /* pointer to solution variable */
} lpitem;

typedef struct ilist {
  struct ilist *prev;
  struct ilist *next;
  struct ilist *b;/* break item, before this item was changed */
  etype  psumb;   /* psumb before this item was changed */
  etype  wsumb;   /* wsumb before this item was changed */
  etype  p;       /* the items profit */
  etype  w;       /* the items weight */
  ntype  no;      /* the items position in a p/w ordering */
} itemlist;

typedef struct rk {
  struct ilist *b;/* current break item */
  etype  psumb;   /* current profit sum up to break item */
  etype  wsumb;   /* current weight sum up to break item */
  itemlist it[MSIZE];
  itemlist head;
  itemlist tail;
} rowknap;
 
typedef int (*funcptr) (const void *, const void *);


/* ======================================================================
 			    global variables
   ====================================================================== */

int n, c, nbobj;
float z, BORNESDP;
float p[MSIZE][MSIZE];
int w[MSIZE];
int x[MSIZE];
etype pbis[MSIZE][MSIZE];
itype wbis[MSIZE];
itype xbis[MSIZE];
int xred[MSIZE];
etype pover[MSIZE];
etype minw[MSIZE];
boolean xstar[MSIZE];
boolean xsol[MSIZE]; /* Solution heuristique */
rowknap rowk[MSIZE];
etype z, z0, ptot, heur, impr, zbrute, lagbound, inibound;
etype wtot;
etype fixp;
etype fixw;
long maxlagr;
long iterations;
stype compt; /* modification LL: ajout de compt */
etype NKR; /* modification LL: ajout de NKR */
etype ustar;
itype indicemin;
etype epsilon=0.0001;
ntype nbobjets;
itype ERREUR;
int opt_H_Sur=0;
int opt_H_Sur_bis;
FILE *trace;
FILE *matq;
FILE *matl; // todo: add l.txt
FILE *matc;
FILE *matA;
FILE *matb;
FILE *matAbis;
FILE *matbbis;
double t1;

void error(char *str);
int povercomp(ntype *a, ntype *b);
int icomp(ntype *a, ntype *b);
int ecomp(lpitem *a, lpitem *b);
void pfree(void *p);
void *palloc(long size, long no);
void partsort(lpitem *f, lpitem *l, etype c);
etype lpsolve(lpitem *a, etype c, int n);
void findpover(void);
void findorder(ntype *t);
void findminw(void);
void improve(int *xprime);
void greedy(void);
void reorder(ntype *new, ntype *old);
static void invertorder(int *rev, int *org);
void iterate(int loopno);
etype findbound(int f, float c, int h, int val);
ntype reduce(ntype *ord);
void decrease(ntype nf, etype *ps, etype *ws);
void increase(int nf);
void initsets(void);
void removeitem(int t);
void insertitem(int t);
etype breakbound(etype ps, etype ws, int t);
void quadbranch(etype ps, etype ws, int t);
int quadknap_Sur(int no, float *cap, float *ptab, float *wtab, int *xtab, float *zqkp, int *stop, int nbobj, int capa, int *winit, float *zinit);
void improve_EkQKP(int *xprime);
void greedy_EkQKP(void);
void sous_gradient(int loopnb);
float quadknap_EkQKP_Sur(int no, int cap, float *ptab, int *wtab, int *xtab, int nbobj, int *H_Sur_opt);//, int z_glob)
float quadknap_EkQKP_Sur_bis(int no, int cap, float *ptab, int *wtab, int *xtab, int nbobj, int *H_Sur_opt, double Epsilon);//, int z_glob)
void improve_EkQKP_end(int *xprime);
float improve_EkQKP_hybrid(int no, int cap, float *ptab, int *wtab, int *xtab, int nbobj, float z_glob);
/*void cpu_time(double *t);*/
void terminate_qkp(char *str);
void printitems(void);
void maketest(int n1, int r, int pct);
void sumdata(int n1, int r1, int pct1, float z, long c, double tottime);
void checksolution(int c, float z);

void load_data(int *n_) ; 

/* ======================================================================
				main
   ====================================================================== */

int main(int argc, char *argv[])
{
  int n, r, pct, v, i, j, fix, k, l;//, z_H_Pri;
  float z, z_H_Sur;
  double time;
  FILE *sol;
  double a=0.0;
  int no_nul;
  static int diagx[MSIZE];
  int fin=0;
  int nred, zsur=0;
  static float pred[MSIZE][MSIZE];
  static int wred[MSIZE];
  static float psdp[MSIZE][MSIZE];
  static int wsdp[MSIZE];
  static double xsdp[MSIZE];
  static int tab_indice[MSIZE];
  int cred, wsum, kred;
  float psum, zred;
  double psumfrac;
  int tempo=0;
  int nbfix1=0;
  double Epsilon=0.05;
  double Epsilon_bis=0.0;
  int opt_H_Sur, opt_H_Sur_bis=0;
  int diff_fix, prev_fix;
  float zbest=0.0;

  load_data(&n) ;

  // printf("n = %d\n", n) ;
  // printf("matrix p : \n");
  // for(i = 0; i<n; i++)
  //   for(j = 0; j<n; j++)
  //       printf("p[%d][%d] = %lf\n", i, j, p[i][j]) ; 
  // printf("w = [");
  // for(i = 0; i<n; i++) printf("%d ", w[i]) ; 
  // printf("] \n") ; 

  // printf("c = %d\n", c) ; 
  // printf("nbobj = %d\n", nbobj) ; 


  tempo=0;
  Epsilon_bis=0.00;
  opt_H_Sur_bis=0;
  zbest=0.0;

  z_H_Sur = quadknap_EkQKP_Sur_bis(n, c, p, w, x, nbobj, &opt_H_Sur_bis, Epsilon_bis);
  zbest=z_H_Sur;
  memcpy(xsol, xstar, MSIZE*sizeof(boolean));

  zred=0.0;
  diff_fix=0; prev_fix=0;
  
    if (z_H_Sur != -1)
    {
      matq = fopen("q.txt", "w");
      matc = fopen("c.txt", "w");
      matA = fopen("A.txt", "w");
      matb = fopen("b.txt", "w");
      matAbis = fopen("Abis.txt", "w");
      matbbis = fopen("bbis.txt", "w");
      matl = fopen("l.txt", "w");
      fprintf(matl, "%d \n", 0) ; 
      fclose(matl);

      /* Ecriture dans q.txt*/
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

      /* Ecriture dans c.txt */
      fprintf(matc,"%d\n",1);
      for (i=1;i<n;i++)
        fprintf(matc,"%f\n",-p[i-1][i-1]);
      fprintf(matc,"%f ",-p[n-1][n-1]);

      /* Ecriture dans A.txt*/
      fprintf(matA,"%d\n",1);
      for (i=1;i<n;i++){
        fprintf(matA,"%d %d ",1,i);
        fprintf(matA,"%d\n",1);
      }
      fprintf(matA,"%d %d ",1,n);
      fprintf(matA,"%d ",1);

      /* Ecriture dans b.txt*/
      fprintf(matb,"%d ",nbobj);

      /* Ecriture dans Abis.txt*/
      fprintf(matAbis,"%d\n",1);
      for (i=1;i<n;i++){
        fprintf(matAbis,"%d %d ",1,i);
        fprintf(matAbis,"%d\n",w[i-1]);
      }
      fprintf(matAbis,"%d %d ",1,n);
      fprintf(matAbis,"%d ",w[n-1]);

      /* Ecriture dans bbis.txt*/
      fprintf(matbbis,"%d ",c);

      fclose(matq);
      fclose(matA);
      fclose(matc);
      fclose(matb);
      fclose(matAbis);
      fclose(matbbis);
      cred=c; nred=n; kred=nbobj;
      for(i=0;i<n;i++)
      {
        for(j=0;j<n;j++)
	      {
          psdp[i][j]=pred[i][j]=p[i][j];
	      }
        wsdp[i]=wred[i]=w[i];
        xsdp[i]=2.0;
        tab_indice[i]=i;
        diagx[i+1]=2;
      }

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

      while (nred)
//        if (!opt_H_Sur_bis)
      {
       cred=c; nred=n; kred=nbobj;
       for(i=0;i<n;i++)
       {
         for(j=0;j<n;j++)
	       {
           psdp[i][j]=pred[i][j]=p[i][j];
	       }
         wsdp[i]=wred[i]=w[i];
         xsdp[i]=2.0;
         tab_indice[i]=i;
         diagx[i+1]=2;
       }
       Epsilon_bis+=0.01;

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
       fix = k = 0;
       nbfix1=0;

       for (l=0;;)
       {
         fscanf(sol,"%d",&i);
         fscanf(sol,"%d",&j);
         fscanf(sol,"%d",&i);
         fscanf(sol,"%d",&j);
         fscanf(sol,"%lf",&a);
         while (xsdp[k]==1.0 || xsdp[k]==0.0)
           k++;
         if (i==j)
         {
           if (a<Epsilon_bis)
           {
             diagx[i-1]=0;
             xsdp[k]=0.0;
             tab_indice[i-1]=k;
	           fix++;
           }
           else
           {
             xsdp[k]=a;
             tab_indice[i-1]=k;
           }
           k++; l++;
         }
         if (l==nred) break;
       }

       fclose(sol);
       fix = k = 0;
       psum = 0.0;
       wsum = 0;
       for (i = 1; i < nred+1; i++)
       {
         if (diagx[i]==2)
         {
           wred[k]=wsdp[i-1];
           l=0;
           for (j = 1; j < nred+1; j++)
           {
             if (diagx[j]==2)
             {
               pred[k][l] = psdp[i-1][j-1];
               l++;
             }
           }
           k++;
           continue;
         }
         fix++;

         if (!diagx[i])
         {
           continue;
         }
         wsum += wsdp[i-1];
         for (j = 1; j < nred+1; j++)
         {
           if (diagx[j]==1)
           {
             psum += psdp[i-1][j-1];
           }
         }  
       }  
       cred-=wsum;
       zred+=psum;
       nred-=fix;
       diff_fix=fix-prev_fix;
       prev_fix=fix;
       if (nred==0) break;
      
      //Calculer la valeur Z_SDP
       psumfrac = 0.0;
       for (i = 2; i <= n+1; i++) 
       {
         if (xsdp[tab_indice[i-1]]!=1) continue;
         for (j = 2; j <= n+1; j++) 
         {
           if (xsdp[tab_indice[j-1]]==1)
           {
             psumfrac += p[i-2][j-2];
           }
         }
       }
       zred=psumfrac;
       for (i = 2; i <= n+1; i++)
       {
	       xred[tab_indice[i-1]]=0;
       }
      
       zsur = quadknap_EkQKP_Sur_bis(nred, cred, pred, wred, xred, kred, &opt_H_Sur_bis, Epsilon_bis);
            
      //Calculer la valeur globale qui est differente de l'addition!!!!
       l=0;
       for (i = 2; i <= n+1; i++) 
       {
         if (xsdp[tab_indice[i-1]]!=0)
	       {
	         if (xsdp[tab_indice[i-1]]<1)
	         {
	           xsdp[tab_indice[i-1]]=xred[l];
	           l++;
	         }
	       }
       }
       psum = 0;
       for (i = 2; i <= n+1; i++) 
       {
         if (!xsdp[tab_indice[i-1]]) continue;
	       for (j = 2; j <= n+1; j++) 
         {
           if (xsdp[tab_indice[j-1]])
	         {
	           psum += psdp[i-2][j-2];
	         }
         }
       }
       zred=psum;

       if(zred>zbest){
         zbest=zred;

         for (i = 0; i < n; i++) {
            xsol[i]=(int)xsdp[i];
         }
       }
       for (i = 0; i < n; i++) 
       {
	       x[i]=(int)xsdp[i];
       }
       zred = improve_EkQKP_hybrid(n, c, p, w, x, nbobj, zred);

       if(zred>zbest){
         zbest=zred;
         memcpy(xsol, x, MSIZE*sizeof(boolean));
        }
      }
    }
    if(zred>zbest){
      zbest=zred;
      for (i = 0; i < n; i++) 
        xsol[i]=x[i];
    }

  }
    
  z_H_Sur = quadknap_EkQKP_Sur(n, c, p, w, x, nbobj, &opt_H_Sur);
  zred=0.0;
    
  if (z_H_Sur != -1)
  {
    if (!opt_H_Sur)
    {
      cred=c; nred=n; kred=nbobj;
      for(i=0;i<n;i++)
      {
        for(j=0;j<n;j++)
	      {
          psdp[i][j]=pred[i][j]=p[i][j];
	      }
        wsdp[i]=wred[i]=w[i];
        xsdp[i]=2.0;
        tab_indice[i]=i;
        diagx[i+1]=2;
      }

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
      fix = k = 0;
      nbfix1=0;
//       while (j != EOF)
      for (l=0;;)
      {
        fscanf(sol,"%d",&i);
        fscanf(sol,"%d",&j);
        fscanf(sol,"%d",&i);
        fscanf(sol,"%d",&j);
        fscanf(sol,"%lf",&a);
        while (xsdp[k]==1.0 || xsdp[k]==0.0)
          k++;
        if (i==j)
        {
          if (a<Epsilon)
          {
            diagx[i-1]=0;
            xsdp[k]=0.0;
            tab_indice[i-1]=k;
	          fix++;
          }
          else
          {
            xsdp[k]=a;
            tab_indice[i-1]=k;
          }
          k++; l++;
        }
        if (l==nred) break;
      }

      fclose(sol);

      kred-=nbfix1;
      fix = k = 0;
      psum = 0.0;
      wsum = 0;
      for (i = 1; i < nred+1; i++)
      {
        if (diagx[i]==2)
        {
          wred[k]=wsdp[i-1];
//           printf("i= %d k=%d wred=%d wsdp=%d\n",i, k, wred[k],wsdp[i-1]);
          l=0;
          for (j = 1; j < nred+1; j++)
          {
            if (diagx[j]==2)
            {
              pred[k][l] = psdp[i-1][j-1];
            l++;
            }
          }
          k++;
          continue;
        }
        fix++;
//         printf("\nfix=%d et i=%d\n",fix,i);
        if (!diagx[i])
        {
          continue;
        }
        wsum += wsdp[i-1];
        for (j = 1; j < nred+1; j++)
        {
          if (diagx[j]==1)
          {
            psum += psdp[i-1][j-1];
          }
        }
      }
      cred-=wsum;
      zred+=psum;
      nred-=fix;
      
      //Calculer la valeur Z_SDP
      psumfrac = 0.0;
      for (i = 2; i <= n+1; i++) 
      {
        if (xsdp[tab_indice[i-1]]!=1) continue;
        for (j = 2; j <= n+1; j++) 
        {
          if (xsdp[tab_indice[j-1]]==1)
          {
            psumfrac += p[i-2][j-2];
          }
        }
      }
      zred=psumfrac;

      for (i = 2; i <= n+1; i++)
      {
	      xred[tab_indice[i-1]]=0;
      }
      
      //cpu_time(NULL);
      zsur = quadknap_EkQKP_Sur(nred, cred, pred, wred, xred, kred, &opt_H_Sur);
            
      //Calculer la valeur globale qui est differente de l'addition!!!!
      l=0;
      for (i = 2; i <= n+1; i++) 
      {
        if (xsdp[tab_indice[i-1]]!=0)
	      {
	        if (xsdp[tab_indice[i-1]]<1)
	        {
	          xsdp[tab_indice[i-1]]=xred[l];
	          l++;
	        }
	      }
      }
      psum = 0.0;
      for (i = 2; i <= n+1; i++) 
      {
//         printf("xsdp[%d]=%lf\n",tab_indice[i-1],xsdp[tab_indice[i-1]]);
        if (!xsdp[tab_indice[i-1]]) continue;
	      for (j = 2; j <= n+1; j++) 
        {
          if (xsdp[tab_indice[j-1]])
	        {
	          psum += psdp[i-2][j-2];
//             printf("xsdp[%d]=%lf xsdp[%d]=%lf p[%d][%d]=%d\n",tab_indice[i-1],xsdp[tab_indice[i-1]],tab_indice[j-1],xsdp[tab_indice[j-1]],i-2,j-2,psdp[i-2][j-2]);
	        }
        }
      }
      zred=psum;
    }
  }
  if (z_H_Sur>=zred)
  {
    zred=z_H_Sur;
  }
  else
  {
    for (i = 0; i < n; i++) 
    {
	     x[i]=(int)xsdp[i];
    }
  }
//    zred = improve_EkQKP_hybrid(n, c, p, w, x, nbobj, zred);
//    fflush(NULL);
  if(zred>zbest) {
    zbest=zred;
//    printf("zbest5=%f\n", zbest);
//    printf("Mise a jour solution heuristique 5\n");
    memcpy(xsol, x, MSIZE*sizeof(boolean));
  }
//    printf("%d: c %d nbobj %d time %.2lf valeur sol heur BestHybrid: %f\n\n", v, c, nbobj, time, zbest);
  printf("%d: c %d nbobj %d valeur sol heur BestHybrid: %f\n\n", v, c, nbobj, zbest);
  fflush(NULL);

// }
  FILE* f = fopen("heur_sol.data", "w") ; 
  fprintf(f, "val %f \n", zbest) ;

  for (i=0; i<n; i++)
    if(xsol[i]){
      // printf("xsol[%d] = %d\n", i, xsol[i]);
      fprintf(f, "%d %d \n", i+1, xsol[i]) ; 
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



/* =======================================================================
          error
   ======================================================================= */

void error(char *str)
{
//   printf("%s\n", str);
//   printf("PROGRAM IS TERMINATED !!!\n\n");
  exit(-1);
}


/* ======================================================================
         comparisons
   ====================================================================== */

int povercomp(ntype *a, ntype *b) {return DET(pover[*b],w[*b],pover[*a],w[*a]);}
int icomp(ntype *a, ntype *b) { return *a - *b; }
int ecomp(lpitem *a, lpitem *b) { 
  if (a->e == b->e) return 0; 
  return (a->e < b->e ? 1 : -1); 
}


/* ======================================================================
          palloc
   ====================================================================== */

void pfree(void *p)
{
  if (p == NULL) error("freeing null");
  free(p);
}


void *palloc(long size, long no)
{
  long *p;

  size *= no;
  if (size == 0) size = 1;
  if (size != (size_t) size) error("Alloc too big");
  p = malloc(size);
  if (p == NULL) error("No memory");
  return p;
}


/* ======================================================================
        partsort
   ====================================================================== */

void partsort(lpitem *f, lpitem *l, etype c)
{
  register lpitem *i, *j, *m;
  register etype me;
  register etype ws;
  register int d;

  for (;;) {
    compt++; /* modification LL: ajout de compt++ */
    d = l - f + 1;
    if (d > 1) {
      m = f + d / 2;
      if (f->e < m->e) SWAP(f, m);
      if (d > 2) {
        if (m->e < l->e) {
          SWAP(m, l); if (f->e < m->e) SWAP(f, m);
        }
      }
    }
    if (d <= 3) break;
    me = m->e; i = f; j = l; ws = 0;
    for (;;) {
      do { ws+=i->w; i++; } while (i->e > me);
      do {           j--; } while (j->e < me);
      if (i > j) break;
      SWAP(i, j);
    }
    if (ws <= c) { f = i; c -= ws; } else l = i-1;
  }
}


/* ======================================================================
                                  lpsolve
   ====================================================================== */

etype lpsolve(lpitem *a, etype c, int n)
{
  register lpitem *k, *m;
  register etype r;
  register etype ps;

  if (c <= 0) return 0.0;
  m = a+n;
  partsort(a, m-1, c);
  for (k = a, ps = 0.0, r = c; ; k++) {
    if (k == m) return ps;
    if (k->w > r) break;
    r -= k->w; ps += k->p;
    if (k->x != NULL) *(k->x) = 1.0;
  }
  if (k->x != NULL) *(k->x) = r / (etype) k->w;
  return ps + r * k->e;
}


/* ======================================================================
        findpover
   ====================================================================== */

void findpover(void)
{
  register int i, j;
  register lpitem *k;
  lpitem a[MSIZE];

  for (i = 0; i < n; i++) {
    /* copy profits with reference to the original item */
    for (j = 0, k = a; j < n; j++) {
      if (j == i) continue;
      k->p = p[j][i]; k->w = w[j];
      k->e = k->p / (etype) k->w; k->x = NULL;
      k++;
    }
    pover[i] = p[i][i] + lpsolve(a, c-w[i], n-1);
  }
}


/* ======================================================================
                                findorder
   ====================================================================== */

void findorder(ntype *t)
{
  int i;

  /* find upper planes */
  findpover();

  /* sort according to pover2 */
  for (i = 0; i < n; i++) t[i] = i;
  qsort(t, n, sizeof(etype), (funcptr) povercomp);
}


/* ======================================================================
          findminw
   ====================================================================== */

void findminw(void)
{
  int i;
  etype mw;

  mw = w[n-1];
  for (i = n-1; i >= 0; i--) {
    if (w[i] < mw) mw = w[i];
    minw[i] = mw;
  }
}


/* ======================================================================
                                  improve
   ====================================================================== */

void improve(int *xprime)
{
  register int i, j, gaini, gainj;
  register etype tot, gain, bgain;
  register etype res;
  etype q[MSIZE];

  res = c;
  /* printf("improve: n %d res %f c %d\n", n, res, c); */
  for (i = 0; i < n; i++) if (xprime[i]) res -= w[i];
  for (;;) {
    for (i = 0; i < n; i++) {
      for (tot = p[i][i], j = 0; j < n; j++) {
        if ((j != i) && (xprime[j] != 0)) tot += p[i][j] + p[j][i];
      }
      q[i] = tot;
    }
    bgain = 0.0;
    gaini = gainj = 0;
    for (i = 0; i < n; i++) {
      if (xprime[i] == 0) {
        if (w[i] <= res) {
          gain = q[i];
          if (gain > bgain) { bgain = gain; gaini = i; gainj = -1; }
        } else {
          for (j = 0; j < n; j++) {
            if (j == i) continue;
            if (xprime[j] == 0) continue;
            if (w[i] - w[j] <= res) {
              gain = q[i] - q[j] - (p[i][j] + p[j][i]);
              if (gain > bgain) { bgain = gain; gaini = i; gainj = j; }
            }
          }
        }
      }
    }
    /* printf("best gain %d i %d j %d res %f\n", bgain, gaini, gainj, res); */
    if (bgain == 0.0) break;
    xprime[gaini] = 1;
    if (gainj != -1) xprime[gainj] = 0; 
    if (gainj != -1) res += w[gainj] - w[gaini]; else res -= w[gaini];
    if (res < 0) error("Negative res");
  } /* end of while loop */
  for (gain = 0, res = c, i = 0; i < n; i++) {
    if (xprime[i] == 0) continue;
//     printf("1)res %f w[%d]=%f\n", res, i, w[i]);
    res -= w[i];
//     printf("2)res %f\n", res);
    for (j = 0; j < n; j++) {
      if (xprime[j]) gain += p[i][j];
    }
  }
  if (res+epsilon < 0) error("1)Excess capacity");
  /* printf("possible solution %d\n", gain); */
  if (gain > z) {
    /* printf("improved solution %d to %d\n", z+fixp, gain+fixp); */
    z = gain;
    for (i = 0; i < n; i++) xstar[i] = xprime[i];
  }
}


/* ======================================================================
          greedy
   ====================================================================== */

void greedy(void)
{
  int i, j, mini;
  etype psum, pi, minp;
  etype wsum;
  double eff, mineff;

  /* find global sums, and initialize solution vector */
  psum = 0.0; wsum = 0;
  for (i = 0; i < n; i++) {
    x[i] = 0; wsum += w[i];
    for (j = 0; j < n; j++) psum += p[i][j];
  }
  ptot = psum;
  wtot = wsum;

  /* run greedy heuristic */
  for (i = 0; i < n; i++) x[i] = 1;
  psum = ptot; wsum = wtot;
  for (;;) {
    mineff = ptot; mini = -1;
    for (i = 0; i < n; i++) {
      if (!x[i]) continue;
      pi = -p[i][i];
      for (j = 0; j < n; j++) if (x[j]) pi += p[j][i] + p[i][j];
      eff = pi / (double) w[i];
      if (eff < mineff) { mineff = eff; mini = i; minp = pi; }
    }
    if (mini == -1) error("No item found\n");
    i = mini; x[i] = 0;
    psum -= minp; wsum -= w[i];
    if (wsum <= c) break;
  }
  z = heur = psum;
  memcpy(xstar, x, MSIZE*sizeof(boolean));
  /* printf("greedy psum %d wsum %d c %d\n", psum, wsum, c); */
  improve(x);
}



/* ======================================================================
                                  reorder
   ====================================================================== */

void reorder(ntype *new, ntype *old)
{
  /* makes an ordering according to indices in new, returns a table */
  /* old, which when applied to reorder, brings back the table to origin */
  int i, j, i1, j1;
  etype p1[MSIZE][MSIZE];
  etype w1[MSIZE];
  itype x1[MSIZE];
  ntype o1[MSIZE];
  etype pover1[MSIZE];
  itype xstar1[MSIZE];

  /* copy to new table */
  for (i = 0; i < n; i++) {
    i1 = new[i];
    o1[i] = old[i1]; 
    /* o1[i1] = old[i]; */
    x1[i] = x[i1];
    w1[i] = (float)w[i1];
//    printf("Reorder: w[%d]=%d, w1[%d]=%f\n", i1, w[i1], i, w1[i]);
    pover1[i] = pover[i1];
    xstar1[i] = xstar[i1];
    for (j = 0; j < n; j++) {
      j1 = new[j];
      p1[i][j] = p[i1][j1];
    }
  }

  /* copy new table back */
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) p[i][j] = p1[i][j];
    xstar[i] = xstar1[i]; pover[i] = pover1[i]; 
    x[i] = x1[i]; w[i] = (int)w1[i];
//    printf("Reorder: w[%d]=%d\n", i, w[i]);
    old[i] = o1[i];
  }

}


/* ======================================================================
                                invertorder
   ====================================================================== */

static void invertorder(int *rev, int *org)
{
  int i;

  for (i = 0; i < n; i++) {
    rev[org[i]] = i;
  }
}


/* ======================================================================
                                  iterate
   ====================================================================== */

void iterate(int loopno)
{
  register int i, j;
  register float pij;
  register lpitem *k;
  register etype d, da, stp, dp, sumd, mu, *pi, *xi;
  register itype h, count;
  register etype zlpmin;
  etype zlp;
  lpitem a[MSIZE];
  etype pp[MSIZE][MSIZE];
  etype pbest[MSIZE][MSIZE];
  etype xover[MSIZE];
  etype x[MSIZE][MSIZE];
  etype delta[MSIZE][MSIZE];
  itype xprime[MSIZE];
  etype pover[MSIZE];

  /* copy to floats */
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) pp[i][j] = p[i][j];
  }

  /* iterate multipliers */
  zlpmin = INFTY; count = 0; mu = 1.0;
  for (h = 0; h < maxlagr; h++) {

    /* solve continuous knapsack for each row */
    for (i = 0; i < n; i++) {
      pi = pp[i]; xi = x[i]; xi[i] = 1.0; /* diagonal always chosen */
      for (j = 0, k = a; j < n; j++) {
        if (j == i) continue;
        k->p = pi[j]; k->w = w[j]; k->e = k->p / k->w;
        k->x = &xi[j]; xi[j] = 0.0;
        k++;
      }
      pover[i] = p[i][i] + lpsolve(a, c-w[i], n-1);
    }

    /* solve column knapsack problem */
    for (i = 0, k = a; i < n; i++, k++) {
      k->p = pover[i]; k->w = w[i]; k->e = k->p / k->w;
      k->x = &xover[i]; xover[i] = 0.0;
    }
    zlp = lpsolve(a, c, n);
    /* if (h % 50 == 0) printf("%d: bound %d\n", h, zlp+fixp);  */
    if ((h == 0) & (loopno == 1)) inibound = zlp+fixp; /* initial ub */

    /* find gradient */
    sumd = 0;
    for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) {
        d = x[i][j] * xover[i] - x[j][i] * xover[j];
        delta[i][j] = d; sumd += d * d;
      }
    }

    /* iterate pp[i][j] */
    for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) {
        d = delta[i][j];
        if (ABS(d) < EPSILON) continue;
        stp = mu * (zlp - z) / sumd;
        dp  = d * stp;
        if (dp > 0) {
          if (dp > pp[i][j]) dp = pp[i][j];
        } else {
          if (-dp > pp[j][i]) dp = -pp[j][i];
        }
        pp[i][j] -= dp; pp[j][i] += dp;
/*      if (h >= 199) // modification LL: ajout de condition
        printf(" pp[%d][%d]=%f et pp[%d][%d]=%f\n",i,j,pp[i][j],j,i,pp[j][i]);fflush(NULL); // modification LL: ajout printf*/
      }
    }
    count++;
    if (zlp < zlpmin) {
      zlpmin = zlp; count = 0;
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) pbest[i][j] = pp[i][j];
      }
    }
    if (count == 20) { mu /= 2; count = 0; }
    if (h % 2 == 0) {
      /* printf("heuristic for h %d\n", h); */
      for (i = 0; i < n; i++) xprime[i] = (int) xover[i]; 
      improve(xprime);
    }
  } /* end of one iteration */
  NKR += (200.0+(double) n)*((double) n+1.0); /* modification LL: ajout de l'affectation */
//   printf("compt=%d\n",compt); /* modification LL: ajout de printf */
//  printf("compt=%d n=%d et nb_iter=%f\n", compt,n,(double) compt/((200.0+(double) n)*((double) n+1.0))); /* modification LL: ajout de printf */
  /* copy to integer table */
  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      pij = p[i][j] + p[j][i];
      p[i][j] = pbest[i][j] + 0.5;
      p[j][i] = pij - p[i][j];
    }
  }
  lagbound = zlpmin+fixp; z0 = z+fixp;
  /* printf("best bound %f greedy z %f impr z %f\n", lagbound,heur,z+fixp); */
}


/* ======================================================================
                                  findbound
   ====================================================================== */

etype findbound(int f, float c, int h, int val)
{
  /* find continuous bound for items i >= t,   */
  /* capacity c, and item h fixed at value val */

  register int i, j;
  register lpitem *k;
  register etype *pi;
  register stype n1;
  register etype p1;
  register etype c1;
  lpitem a[MSIZE];
  etype pover[MSIZE];

  /* solve continuous knapsack for each row */
  for (i = f; i < n; i++) {
    pi = p[i]; c1 = c; p1 = 0.0; n1 = 0;
    for (j = f, k = a; j < n; j++) {
      if (j == h) { /* fixed variable */
        if (val == 1) { c1 -= w[j]; p1 += pi[j]; } continue;
      }
      if (j == i) { /* diagonal value always chosen */
        c1 -= w[j]; p1 += pi[j]; continue;
      }
      k->p = pi[j]; k->w = w[j]; k->e = k->p / (etype) k->w; k->x = NULL;
      k++; n1++;
    }
    pover[i] = p1 + lpsolve(a, c1, n1);
  }

  /* solve column knapsack problem */
  c1 = c; p1 = 0; n1 = 0;
  for (i = f, k = a; i < n; i++) {
    if (i == h) {
      if (val == 1) { c1 -= w[i]; p1 += pover[i]; } continue;
    }
    k->p = pover[i]; k->w = w[i]; k->e = k->p / (etype) k->w; k->x = NULL;
    k++; n1++;
  }
  return p1 + lpsolve(a, c1, n1);
}


/* ======================================================================
                                   reduce
   ====================================================================== */

ntype reduce(ntype *ord)
{
  register int i, j;
  stype u, r1, r0;

  r1 = r0 = 0;
  for (i = 0; i < n; i++) {
    ord[i] = i;
    if (xstar[i] == 1) {
      u = findbound(0, c, i, 0);
      if (u < z+1) {
        ord[i] = i+n; xstar[i] = 1; r1++;
        /* printf("item %d fixed at 1\n", i); */
      }
    } else {
      u = findbound(0, c, i, 1);
      if (u < z+1) {
        ord[i] = i+n; xstar[i] = 0; r0++;
        /* printf("item %d fixed at 0\n", i); */
      }
    }
  }
  qsort(ord, n, sizeof(ntype), (funcptr) icomp);
  for (i = 0; i < n; i++) if (ord[i] >= n) ord[i] -= n;
  /* printf("fixed %d at 1, %d at 0 out of %d\n", r1, r0, n); */
  return r0 + r1;
}


/* ======================================================================
                                   decrease
   ====================================================================== */

void decrease(ntype nf, etype *ps, etype *ws)
{
  register int i, j;
  etype psum;
  etype wsum;

  psum = wsum = 0.0;
  for (i = n-1; i >= n-nf; i--) {
    if (xstar[i] == 1) {
      for (j = 0; j < i; j++) p[j][j] += p[i][j] + p[j][i];
      psum += p[i][i]; wsum += w[i];
    }
  }
  *ps = psum; *ws = wsum;
}


/* ======================================================================
                                   increase
   ====================================================================== */

void increase(int nf)
{
  register int i, j;

  for (i = n; i < n+nf; i++) {
    if (xstar[i] == 1) {
      for (j = 0; j < i; j++) p[j][j] -= p[i][j] + p[j][i];
    }
  }
}


/* ======================================================================
          initsets
   ====================================================================== */

void initsets(void)
{
  register int i, j, l;
  register lpitem *k;
  register itemlist *u, *h, *g;
  register etype psum;
  register etype wsum, c1;
  lpitem a[MSIZE];
  rowknap *row;

  for (i = 0; i < n; i++) {
    row = rowk + i;
    for (j = 0, k = a, h = row->it; j < n; j++, h++) {
      if (j == i) continue;
      h->p = k->p = p[i][j];
      h->w = k->w = w[j];
      k->e = k->p / (etype) k->w;
      k->x = (etype *) j;
      k++;
    }
    qsort(a, (int) (k-a), sizeof(lpitem), (funcptr) ecomp);
    g = &(row->head);
    for (j = 0, l = 1, k = a; j < n; j++) {
      if (j == i) continue;
      h = row->it + (int) (k->x);
      h->no = l; g->next = h;
      g = h; k++; l++;
    }
    g->next = &(row->tail);

    for (h = &(row->head), g = &(row->tail); h != g; h = h->next) {
      h->next->prev = h;
    }

    /* now finally find break item */
    psum = wsum = 0; c1 = c - w[i];
    for (h = row->head.next, g = &(row->tail); h != g; h = h->next) {
      if (wsum + h->w > c1) break;
      psum += h->p; wsum += h->w;
    }
    row->psumb = psum; row->wsumb = wsum; row->b = h;
    row->head.no = 0; row->head.prev = NULL;
    row->tail.no = n; row->tail.next = NULL;
    row->head.p = 0; row->head.w = 0;
    row->tail.p = 0; row->tail.w = c;
    /* printf("row kp (%d,%d) c %d b (%d,%d)\n", psum,wsum,c1,h->p,h->w); */
  }
}


/* ======================================================================
         removeitem 
   ====================================================================== */

void removeitem(int t)
{
  register itemlist *h;
  register rowknap *row;
  register int i;

  for (i = t+1, row = rowk+t+1; i < n; i++, row++) {
    h = &(row->it[t]);
    /* printf("%d: removing (%d,%d) no %d\n", i, h->p, h->w, h->no); */
    h->next->prev = h->prev;
    h->prev->next = h->next;
    h->psumb = row->psumb;
    h->wsumb = row->wsumb;
    h->b     = row->b;
    if (row->b == h) { /* equals b, move to next */
      row->b = h->next; continue;
    }
    if (h->no < row->b->no) { /* before b, subtract */
      row->psumb -= h->p; row->wsumb -= h->w;
    }
  }
}


/* ======================================================================
         insertitem
   ====================================================================== */

void insertitem(int t)
{
  register itemlist *h;
  register rowknap *row;
  register int i;

  for (i = t+1, row = rowk+t+1; i < n; i++, row++) {
    h = &(row->it[t]);
    /* printf("inserting (%d,%d)\n", h->p, h->w); */
    h->next->prev = h;
    h->prev->next = h;
    row->b     = h->b; 
    row->psumb = h->psumb; 
    row->wsumb = h->wsumb;
  }
}


/* ======================================================================
                                 breakbound
   ====================================================================== */

etype breakbound(etype ps, etype ws, int t)
{
  register itemlist *h;
  register etype psumb;
  register etype wsumb, c1;
  register rowknap *row;
  register lpitem *k;
  register int i;
  lpitem d[MSIZE];
  etype u, psum;

  /* printf("breakbound ps %f, ws %d, c1 %d t %d\n", ps, ws, c-ws, t); */
  for (i = t, k = d, row = rowk+t; i < n; i++, row++) {
    c1 = c - ws - w[i]; if (c1 < 0) continue; /* no space for row */
    h = row->b; psumb = row->psumb; wsumb = row->wsumb;
    if (wsumb > (c1+0.00001)) { /* moving b bacward */
      do { 
        h = h->prev; psumb -= h->p; wsumb -= h->w;
     } while (wsumb > (c1+0.00001));
    } else { /* moving b forward */
      while (wsumb + h->w <= c1) {
        psumb += h->p; wsumb += h->w; h = h->next;
      }
    }
    row->b = h; row->psumb = psumb; row->wsumb = wsumb; 
    /* printf(" %d: found c1 %d b (%d,%d)\n", i, c1, h->p, h->w); */
    k->p = p[i][i] + psumb + (c1 - wsumb) * (double) h->p / h->w;
    k->w = w[i];
    k->e = k->p / (etype) k->w;
    k->x = NULL;
    k++;
  }
  psum = lpsolve(d, c-ws, (k-d));
  u = psum + ps;   
//   printf("derived u %d = %d z= %d\n", u, u+fixp, z);
  return u;
}


/* ======================================================================
                                quadbranch
   ====================================================================== */

void quadbranch(etype ps, etype ws, int t)
{
  register int j;
  etype u;
  double time;

  iterations++;
//   printf("ws=%f c=%f\n", ws, c);
  fflush(NULL);
  if (ws > c+epsilon) error("Excess weight");
/*  cpu_time(&time);
  if(time>30) {
//       printf("Time limit\n");
//       printf("PROGRAM IS TERMINATED !!!\n\n");
      if (ustar == 0 & t == 0) {
        ustar = breakbound (ps, ws, 1);
      }
      ERREUR=1;
      return;
    }*/
// error("Time limit");
//   printf("toto6.0\n");
  if (ps > z) {
//   printf("toto6.1\n");
    z = ps; 
    /* printf("improved solution to %d\n", z); */
    for (j = 0; j < n; j++)
      xstar[j] = x[j];
  }

  if ((t != n) && (ws + minw[t] <= c)) {
//   printf("t=%d\n",t);
    u = breakbound(ps, ws, t);
    if (u > z) {
//       ustar = u+fixp;
//   printf("toto6.3\n");
      removeitem(t);
      if (ws + w[t] <= c) {
        for (j = t+1; j < n; j++) p[j][j] += p[j][t] + p[t][j];
        x[t] = 1;
//   printf("toto6.4\n");
        quadbranch(ps + p[t][t], ws + w[t], t+1);
//       printf("derived u %d = %d z= %d ps = %d\n", u, u+fixp, z, ps);
//   printf("toto6.5\n");
        x[t] = 0;
        for (j = t+1; j < n; j++) p[j][j] -= p[j][t] + p[t][j]; 
      }
      if (t==indicemin && !ERREUR) {
//         printf("t = %d  indicemin = %d\n", t, indicemin);
        ustar = u+fixp;
/*        for (j = 0; j < n; j++)
          printf("x[%d]=%d\t", j, x[j]);
        printf("\n");*/
/*        for (j = 0; j < n; j++) xustar[j] = x[j];
        for (j = 0; j < n; j++)
          printf("xustar[%d]=%d\t", j, xustar[j]);
        printf("\n");*/
        indicemin++;
      }
      if (ustar == 0 && ERREUR && t == 0) {
        ustar = breakbound (ps, ws, 1);
      }
//   printf("toto6.6\n");
      quadbranch(ps, ws, t+1);
      insertitem(t);
//   printf("toto6.7\n");
    }
  }
}


/* ======================================================================
        quadknap
   ====================================================================== */

int quadknap_Sur(int no, float *cap, float *ptab, float *wtab, int *xtab, float *zqkp, int *stop, int nbobj, int capa, int *winit, float *zinit)
{
  int i, j;
  long time, ltime;
  etype ps;
  etype ws;
  ntype nf, fix;
  ntype ord[MSIZE], o1[MSIZE];
  double times; /* modification LL: ajout de times */
  stype objsum, wsum;
  etype psum;
 
  /* initialize global variables */
  n    = no;         /* number of items */
  c    = *cap;        /* capacity */
  z    = 0;          /* incumbent solution */
  nf   = 0;          /* number of variables fixed to some value */
  fixp = 0;          /* quadratic profit sum of variables fixed to 1 */
  fixw = 0;          /* weight sum of variables fixed to 1 */
  iterations = 0;    /* number of b&b nodes */
  maxlagr = 200 + n; /* max number of lagrange iterations */
  NKR=0; /* modification LL: ajout de NKR */
  compt=0; /* modification LL: ajout de compt */
  ERREUR=0;
  *stop = 0;
  ustar = 0;
  indicemin = 0;

  /* copy profits and weights to internal arrays */
  if (no > MSIZE) error("too many items");
  memcpy(p, ptab, sizeof(float) * MSIZE * MSIZE);
  memcpy(w, wtab, sizeof(float) * MSIZE);

//   printf("n=%d\t c=%f\n", n, c);fflush(NULL);
/*  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i==1)
        printf("p[%d][%d]=%d\t", i,j,p[i][j]);fflush(NULL);
    }
  }*/
//   printf("\n");(c1+0.00001)
//   for (i = 0; i < n; i++) {
//     printf("w[%d]=%f\t", i,w[i]);fflush(NULL);
//   }
//   printf("\n");fflush(NULL);

//   printf("toto0\n");
  /* find initial solution */
  greedy();
//   printf("toto1\n");
  /* repeat: tighten bound, reduce */
  for (j = 0; j < n; j++) o1[j] = j;
//  cpu_time(NULL); /* modification LL: ajout de cpu_time */
  for (i = 1;; i++) {
    /* iterate multipliers */
    iterate(i);

//   printf("toto2\n");
    /* perform reduction */
    fix = reduce(ord);
//     printf("%d: NOW REDUCED %d TOTAL REDUCED %d\n", i, fix, nf+fix);
//   printf("toto3\n");
    if (fix == 0) break;
    reorder(ord, o1);
    decrease(fix, &ps, &ws);
    nf += fix; n -= fix; z -= ps; c -= ws; fixp += ps; fixw += ws;
  }
//  cpu_time(&times); /* modification LL: ajout de cpu_time */
//  printf("time %.2lf\n",times); /* modification LL: ajout de printf */
//  printf("compt=%d et nb iter partsort moyen: %f\n",compt, (double) compt / NKR); /* modification LL: ajout de printf */
  /* printf("inibound %d, lagbound %d, greedy z %d\n", inibound,lagbound,z); */

//   printf("toto4\n");
  /* reorder according to descending values of pover */
  findorder(ord);
  reorder(ord, o1);

//   printf("toto5\n");
  /* construct table which makes it possible to derive bounds fast */
  initsets();
  findminw();

//   printf("toto6\n");
  /* branch-and-bound algorithm */
  for (j = 0; j < n; j++) x[j] = 0;
  quadbranch(0, 0, 0); 

//   printf("toto7\n");
  /* restore table */
  /* printf("increase n by %d, z by %d, c by %d\n", nf, fixp, fixw); */
  if (!ERREUR) {
    increase(nf);
    n += nf; z += fixp; c += fixw;
    invertorder(ord, o1);
    reorder(ord, o1);

//   printf("toto8\n");
  /* copy solution vector to xtab */
    memcpy(xtab, xstar, sizeof(int) * MSIZE);
    *zqkp = z;
    objsum = wsum = psum = 0;
    for (j = 0; j < n; j++) {
      if (xtab[j]) {
  objsum ++;
//  printf("w[%d]=%d\t", j, winitial[j]);
  wsum += winit[j];
        for (i = 0; i < n; i++) {
    if (xtab[i]){
      psum += p[j][i];
//             printf("p[%d][%d]=%d\n", j, i, p[j][i]);
    }
  }
      }
    }
    if ((wsum <= capa) & (objsum == nbobj)) *zinit = psum;
//     printf("\nobjsum= %d, nbobj= %d, wsum= %d, cap= %d, z= %d\n", objsum, nbobjets, wsum, capacity, *zinit);
  }
  else {
/*    for (i = 0; i < n; i++)
      printf("xustar[%d]=%d\t", i, xustar[i]);
    printf("\n");*/
//     memcpy(xtab, xustar, sizeof(int) * MSIZE);
/*    for (i = 0; i < n; i++)
      printf("xtab[%d]=%d\t", i, xtab[i]);
    printf("\n");
 */
   *stop=-1;
//     *zqkp = z;
    increase(nf);
    n += nf; z += fixp; c += fixw;
    invertorder(ord, o1);
    reorder(ord, o1);

//   printf("toto8\n");
  /* copy solution vector to xtab */
    memcpy(xtab, xstar, sizeof(int) * MSIZE);
    *zqkp = z;
/*    *zqkp = ustar;*/
//     printf("zqkp_best = %d\n", z);
//     printf("ustar = %d indicemin = %d\n", ustar, indicemin);
    objsum = wsum = psum = 0;
    for (j = 0; j < n; j++) {
      if (xtab[j]) {
  objsum ++;
//  printf("w[%d]=%d\t", j, winitial[j]);
  wsum += winit[j];
        for (i = 0; i < n; i++) {
    if (xtab[i]){
      psum += p[j][i];
//             printf("p[%d][%d]=%d\n", j, i, p[j][i]);
    }
  }
      }
    }
    if ((wsum <= capa) & (objsum == nbobj)) *zinit = psum;
//     printf("\nobjsum= %d, nbobj= %d, wsum= %d, cap= %d, z= %d\n", objsum, nbobjets, wsum, capacity, *zinit);
  }

  return 1;
}


/* ======================================================================
                                  improve
   ====================================================================== */

void improve_EkQKP(int *xprime)
{
  int i, j, gaini, gainj;
  etype tot, gain, bgain;
  stype res, kres;
  etype q[MSIZE];
  double eff, maxeff;
  int maxi;
  etype pi;

  res = c; kres = nbobjets;
  // printf("improve: n %d res %d c %d\n", n, res, c); 
  for (i = 0; i < n; i++) {
    if (xprime[i]) {
      res -= wbis[i];
      kres--;
    }
  }
  for (;;) {
    for (i = 0; i < n; i++) {
      for (tot = pbis[i][i], j = 0; j < n; j++) {
        if ((j != i) && (xprime[j] != 0)) tot += pbis[i][j] + pbis[j][i];
      }
      q[i] = tot;
    }
    bgain = gaini = gainj = 0;
    for (i = 0; i < n; i++) {
      if (xprime[i] == 0) {
        if (wbis[i] <= res && kres > 0) {
          gain = q[i];
          if (gain > bgain) { bgain = gain; gaini = i; gainj = -1; }
        } else {
          for (j = 0; j < n; j++) {
            if (j == i) continue;
            if (xprime[j] == 0) continue;
            if (wbis[i] - wbis[j] <= res) {
              gain = q[i] - q[j] - (pbis[i][j] + pbis[j][i]);
              if (gain > bgain) { bgain = gain; gaini = i; gainj = j; }
            }
          }
        }
      }
    }
//     printf("best gain %d i %d j %d res %d\n", bgain, gaini, gainj, res);
    if (bgain == 0) break;
    xprime[gaini] = 1; kres--;
    if (gainj != -1) {xprime[gainj] = 0; kres++;}
    if (gainj != -1) res += wbis[gainj] - wbis[gaini]; else res -= wbis[gaini];
    if (res < 0) {
//       printf("Negative res\n");
//       printf("PROGRAM IS TERMINATED !!!\n\n");
      ERREUR=1;
      return;
    }
//       error("Negative res");
  } // end of while loop 

  for (gain = 0, res = c, kres = nbobjets, i = 0; i < n; i++) {
    if (xprime[i] == 0) continue;
//     printf("1)wi=%d res=%d kres=%d\n",wbis[i], res, kres);
    res -= wbis[i]; kres--;
//     printf("2)wi=%d res=%d kres=%d pii=%d\n",wbis[i], res, kres, pbis[i][i]);
    for (j = 0; j < n; j++) {
      if (xprime[j]) gain += pbis[i][j];
    }
  }
  if (res < 0 || kres < 0) {
//       printf("2)Excess capacity\n");
//       printf("PROGRAM IS TERMINATED !!!\n\n");
      ERREUR=1;
      return;
    }
// error("Excess capacity");
  if (kres) {
    for (;;) {
      maxeff = 0; maxi = -1;
      for (i = 0; i < n; i++) {
        if (xbis[i]) continue;
        pi = -pbis[i][i];
        for (j = 0; j < n; j++) if (xbis[j]) pi += pbis[j][i] + pbis[i][j];
        eff = pi / (double) wbis[i];
        if (eff > maxeff && wbis[i] <= res) { maxeff = eff; maxi = i; }
      }
      if (maxi == -1) {
//       printf("No item found\n");
//       printf("PROGRAM IS TERMINATED !!!\n\n");
      ERREUR=1;
      return;
    }
// error("No item found\n");
      i = maxi; xbis[i] = 1;
      res -= wbis[i]; kres--;
      if (res >= 0 && kres == 0) break;
    }
  }
  for (gain = 0, res = c, kres = nbobjets, i = 0; i < n; i++) {
    if (xprime[i] == 0) continue;
    res -= wbis[i]; kres--;
    for (j = 0; j < n; j++) {
      if (xprime[j]) gain += pbis[i][j];
    }
  }
  if (res < 0 || kres < 0 || kres > 0) {
//       printf("3)Excess capacity\n");
//       printf("PROGRAM IS TERMINATED !!!\n\n");
      ERREUR=1;
      return;
    }
// error("Excess capacity");
//   printf("possible solution %d\n", gain);
  if (gain > z) {
//     printf("improved solution %d to %d\n", z+fixp, gain+fixp);
    z = gain;
    for (i = 0; i < n; i++) xstar[i] = xprime[i];
  }
}


/* ======================================================================
          greedy
   ====================================================================== */

void greedy_EkQKP(void)
{
  int i, j, mini, maxi;
  etype psum, pi, minp;
  stype wsum;
  double eff, mineff;

  // find global sums, and initialize solution vector 
  psum = 0; wsum = 0; maxi =0;
  for (i = 0; i < n; i++) {
    xbis[i] = 0; wsum += wbis[i];
    for (j = 0; j < n; j++) psum += pbis[i][j];
  }
  ptot = psum;
  wtot = wsum;
  maxi = n;

  // run greedy heuristic 
  for (i = 0; i < n; i++) xbis[i] = 1;
  psum = ptot; wsum = wtot;
  for (;;) {
    mineff = ptot; mini = -1;
    for (i = 0; i < n; i++) {
      if (!xbis[i]) continue;
      pi = -pbis[i][i];
      for (j = 0; j < n; j++) if (xbis[j]) pi += pbis[j][i] + pbis[i][j];
      eff = pi / (double) wbis[i];
      if (eff < mineff) { mineff = eff; mini = i; minp = pi; }
    }
    if (mini == -1) {
//       printf("4)Excess capacity\n");
//       printf("PROGRAM IS TERMINATED !!!\n\n");
      ERREUR=1;
      return;
    }
// error("No item found\n");
    i = mini; xbis[i] = 0;
    psum -= minp; wsum -= wbis[i]; maxi--;
    if (wsum <= c && maxi <= nbobjets) break;
  }
  z = heur = psum;
  memcpy(xstar, xbis, MSIZE*sizeof(boolean));
//   printf("greedy psum %d wsum %d c %d maxi %d nbobjets %d\n", psum, wsum, c, maxi, nbobjets);
  improve_EkQKP(xbis);
/*  for (i = 0; i < n; i++) 
  {
    if (x[i])
      printf("greedy: x[%d]=%d\n", i, x[i]);
    if (xstar[i])
      printf("greedy: xstar[%d]=%d\n", i, xstar[i]);
    if (xbis[i])
      printf("greedy: xbis[%d]=%d\n", i, xbis[i]);
  }*/

//   printf("greedy improved psum %d wsum %d c %d maxi %d nbobjets %d\n", psum, wsum, c, maxi, nbobjets);
}



/* ======================================================================
                                  sous_quadknap_E-kQKP_LR1gradient
   ====================================================================== */

void sous_gradient(int loopnb)
{
  int i, j, count;
  ptype mu, mul, mur, alpha, alphal, alphar;
  itype h, s;
  etype zqkpmin, zqkp, zqkp0, zqkp1, stop, zqkpbis;
  int status, end=0;
  int permute =0;
  etype zopt = 0;
  
  etype w1[MSIZE];
  etype cflo;
  stype objsum, wsum, objsum1, objsum0, objsuminf, wsum1, wsum0, wsuminf;

  itype xover[MSIZE];
  itype xprime[MSIZE];
  double time;
  int stoppe = 0;

  zqkpmin = INFTY; count = 0;

  /*solve Sur with mu=1*/
  for (i = 0; i < n; i++)
    xover[i] = 0;
  zqkp1=0.0;
  mu = 1.0;
  for (i = 0; i < n; i++)
    w1[i] = (etype) wbis[i] + mu;
  cflo = (etype) c + mu*nbobjets;

/*  for (i = 0; i < n; i++)
  {
    printf("i= %d w1=%f\n",i, w1[i]);  
    for (j=0;j<n;j++)
    {
      if (i==49)
        printf("i= %d j=%d p=%d\n",i, j, pbis[i][j]);  
    }
  }*/

  status = quadknap_Sur(n, &cflo, pbis, w1, xover, &zqkp1, &stop, nbobjets, c, wbis, &zopt);
  if (zopt > z0)
  {
    z0 = zopt;
    memcpy(xstar, xover, sizeof(int) * MSIZE);
  }
  if(zopt >= zqkp1) opt_H_Sur=1;
  count++;
  if (stop == -1) {
/*    wsum1 = 0;
    for (i = 0; i < n; i++)
    if (xover[i]) wsum1 +=wbis[i];
    objsum1 = 0;
    for (i = 0; i < n; i++)
      if (xover[i]) objsum1 ++;
    printf("\nzqkp1= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsum1 = %d et wsum1 = %d et end = %d\n", zqkp1, mu, c, cflo, nbobjets, objsum1, wsum1, end);*/
//     fflush(NULL);
//     printf("\nzqkp1= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et end = %d\n", zqkp1, zopt, mu, c, cflo, nbobjets, end);
//     fflush(NULL);
    return;
  }

  objsum1 = 0;
  for (i = 0; i < n; i++)
    if (xover[i]) objsum1 ++;
  wsum1 = 0;
  for (i = 0; i < n; i++)
    if (xover[i]) wsum1 +=wbis[i];
//   fflush(NULL);
//   printf("\nzqkp1= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsum1 = %d et wsum1 = %d et end = %d\n", zqkp1, zopt, mu, c, cflo, nbobjets, objsum1, wsum1, end);
//   fflush(NULL);
  if (zqkp1 < zqkpmin) zqkpmin = zqkp1;
  if ((wsum1 <= c) & (objsum1 <= nbobjets)) {    /*Cas a traiter: wsum <= c et objsum < nbobjets car cont2 en egalite*/
//     if (zqkp1 < zqkpmin) zqkpmin = zqkp1;
    end = 1;
  }
  else {
    alpha = (etype) (c - wsum1) / (etype) (objsum1 - nbobjets);
//     printf("alpha = %lf et objsum1 = %d et wsum1 = %d et end = %d\n", alpha, objsum1, wsum1, end);
    if (wsum1 > c) {
      for (i = 0; i < n; i++)
        xover[i] = 0;
      zqkp0=0.0;
      mu = 0.0;
      for (i = 0; i < n; i++)
        w1[i] = (etype) wbis[i];
      cflo = (etype) c;

//      cpu_time(&time);
//       printf(" time %.2lf\n", time);
//      cpu_time(NULL);

      status = quadknap_Sur(n, &cflo, pbis, w1, xover, &zqkp0, &stop, nbobjets, c, wbis, &zopt);
      if (zopt > z0)
      {
        z0 = zopt;
        memcpy(xstar, xover, sizeof(int) * MSIZE);
      }
      if(zopt >= zqkp0) opt_H_Sur=1;
      count++;
      if (stop == -1) {
        wsum0 = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) wsum0 +=wbis[i];
        objsum0 = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) objsum0 ++;
//         fflush(NULL);
//         printf("\nzqkp0= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsum0 = %d et wsum0 = %d et end = %d\n", zqkp0, zopt, mu, c, cflo, nbobjets, objsum0, wsum0, end);
//         printf("\nzqkp0= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et end = %d\n", zqkp0, zopt, mu, c, cflo, nbobjets, end);
//         fflush(NULL);
        return;
      }

      objsum0 = 0;
      for (i = 0; i < n; i++)
        if (xover[i]) objsum0 ++;
      wsum0 = 0;
      for (i = 0; i < n; i++)
        if (xover[i]) wsum0 +=wbis[i];
//       fflush(NULL);
//       printf("\nzqkp0= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsum0 = %d et wsum0 = %d et end = %d\n", zqkp0, zopt, mu, c, cflo, nbobjets, objsum0, wsum0, end);
//       fflush(NULL);
      if (zqkp0 < zqkpmin) zqkpmin = zqkp0;
      if (objsum0 <= nbobjets) {    /*Cas a traiter: wsum <= c et objsum < nbobjets car cont2 en egalite*/
//         if (zqkp0 < zqkpmin) zqkpmin = zqkp0;
        end = 1;
      }
      else {
        alphal = (etype) (c - wsum0) / (objsum0 - nbobjets);
        alphar = alpha;
      }
    }
    else {
      for (i = 0; i < n; i++)
        xover[i] = 0;
      zqkp=0.0;
      mu = 0.0;
      for (i = 0; i < n; i++)
        w1[i] = mu*wbis[i] + (1.0-mu);
      cflo = (etype) mu*c + (1.0-mu)*nbobjets;

//      cpu_time(&time);
//       fflush(NULL);
//       printf(" time %.2lf\n", time);
//       fflush(NULL);
//      cpu_time(NULL);

      status = quadknap_Sur(n, &cflo, pbis, w1, xover, &zqkp, &stop, nbobjets, c, wbis, &zopt);
      if (zopt > z0)
      {
        z0 = zopt;
        memcpy(xstar, xover, sizeof(int) * MSIZE);
      }
      if(zopt >= zqkp) opt_H_Sur=1;
/*      for (i = 0; i < n; i++)
        printf("xover[%d]=%d\t", i, xover[i]);
      printf("\n");*/
      count++;
      wsuminf = 0;
      for (i = 0; i < n; i++)
        if (xover[i]) wsuminf +=wbis[i];
      objsuminf = 0;
      for (i = 0; i < n; i++)
        if (xover[i]) objsuminf ++;
//       fflush(NULL);
//       printf("\nzqkpinf= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsuminf = %d et wsuminf = %d et end = %d\n", zqkp, zopt, mu, c, cflo, nbobjets, objsuminf, wsuminf, end);
//       fflush(NULL);
      if (zqkp < zqkpmin) zqkpmin = zqkp;
      if (wsuminf <= c & objsuminf <= nbobjets && !stop) {    /*Cas a traiter: wsum <= c et objsum < nbobjets car cont2 en egalite*/
        stop = 0;
      }
      if (zqkp <= z) {
        stop = 0;
        end = 1;
      }

      while (stop == -1)  {
        stoppe = 1;
        zqkpbis=zqkp;
        for (i = 0; i < n; i++)
          xover[i] = 0;
        zqkp=0;
        stop=0,
        mu += 0.01;
        for (i = 0; i < n; i++)
          w1[i] = mu*wbis[i] + (1.0-mu);
        cflo = (etype) mu*c + (1.0-mu)*nbobjets;

//        cpu_time(&time);
//         fflush(NULL);
//         printf(" time %.2lf\n", time);
//         fflush(NULL);
//        cpu_time(NULL);
        status = quadknap_Sur(n, &cflo, pbis, w1, xover, &zqkp, &stop, nbobjets, c, wbis, &zopt);
        if (zopt > z0)
        {
          z0 = zopt;
          memcpy(xstar, xover, sizeof(int) * MSIZE);
        }
        if(zopt >= zqkp) opt_H_Sur=1;
        count++;
        wsuminf = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) wsuminf +=wbis[i];
        objsuminf = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) objsuminf ++;
//         fflush(NULL);
//         printf("\nzqkpinf= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsuminf = %d et wsuminf = %d et end = %d\n", zqkp, zopt, mu, c, cflo, nbobjets, objsuminf, wsuminf, end);
//         fflush(NULL);
        if (zqkp < zqkpmin) zqkpmin = zqkp;
        if (wsuminf <= c & objsuminf <= nbobjets) {    //Cas a traiter: wsum <= c et objsum < nbobjets car cont2 en egalite
//           fflush(NULL);
//           printf("break\n");
//    fflush(NULL);
          break;
        }
        if (zqkpbis < zqkp && !stop) {
//           fflush(NULL);
//           printf("zqkpbis %d zqkp %d\n", zqkpbis, zqkp);
//           fflush(NULL);
    end = 1;
          break;
        }
        if (zqkp <= z) {
//           fflush(NULL);
//           printf("zqkp %d z %d\n", zqkp, z);
//           fflush(NULL);
    end = 1;
          break;
        }
        if (count >= 5 && stop) {
//           fflush(NULL);
//           printf("count %d\n", count);
//           fflush(NULL);
    end = 1;
          break;
        }
      }

      if (wsuminf <= c & objsuminf <= nbobjets & !stop) {    /*Cas a traiter: wsum <= c et objsum < nbobjets car cont2 en egalite*/
        wsuminf = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) wsuminf +=wbis[i];
        objsuminf = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) objsuminf ++;
//         fflush(NULL);
//         printf("\nzqkpinf= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsuminf = %d et wsuminf = %d et end = %d\n", zqkp, zopt, mu, c, cflo, nbobjets, objsuminf, wsuminf, end);
//         fflush(NULL);
        if (zqkp < zqkpmin) zqkpmin = zqkp;
        end = 1;
      }
      else {
        if (stop) alphal = 0.0;
        else alphal = (etype) (objsuminf - nbobjets) / (c - wsuminf);
        alphar = 1.0/alpha;
        permute = 1;
      }
    }
  }

//   printf("alphal = %lf et alphar = %lf et end = %d\n", alphal, alphar, end);
  mul = 0.0;
  mur = 1.0;
  stoppe = 0;
  /* iterate multipliers */

//   for (h = 0; h < maxlagr; h++) {
  while (!end) {
    if (alphal >= alphar) {
      /*for (i = 0; i < n; i++)
        xover[i] = 0;
      zqkp1=0.0;
      mu = mul;
      for (i = 0; i < n; i++)
        w1[i] = (etype) wbis[i] + mu;
      cflo = (etype) c + mu*nbobjets;

      cpu_time(&time);
      status = quadknap_Sur(n, &cflo, p, w1, xover, &zqkp1, nbobjets, c, w, &zopt);
      count++;
      if (zqkp1 == -1)  {
        wsum = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) wsum +=wbis[i];
        objsum = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) objsum ++;
        printf("\nzqkp1= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsum = %d et wsum = %d et end = %d\n", zqkp1, mu, c, cflo, nbobjets, objsum, wsum, end);
        return;
      }

      printf("\nzqkp1= %d, mu= %f, c= %d, cflo= %f nbobj= %d\n", zqkp1, mu, c, cflo, nbobjets);
      for (i = 0; i < n; i++)
        xover[i] = 0;
      zqkp0=0.0;
      mu = mur;
      for (i = 0; i < n; i++)
        w1[i] = (etype) wbis[i] + mu;
      cflo = (etype) c + mu*nbobjets;

      cpu_time(&time);
      status = quadknap_Sur(n, &cflo, p, w1, xover, &zqkp0, nbobjets, c, w, &zopt);
      count++;
      if (zqkp0 == -1) {
        wsum = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) wsum +=wbis[i];
        objsum = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) objsum ++;
        printf("\nzqkp0= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsum = %d et wsum = %d et end = %d\n", zqkp0, mu, c, cflo, nbobjets, objsum, wsum, end);
        return;
      }

      if (zqkp0 <= zqkp1) zqkpmin = zqkp0;
      else zqkpmin = zqkp1;
      */
      end = 1;
//       fflush(NULL);
//       printf("\nalphal >= alphar: zqkp= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et end = %d\n", zqkpmin, zopt, mu, c, cflo, nbobjets, end);
//       fflush(NULL);
    }
    else {
      mu = (alphal + alphar)/2;
      for (i = 0; i < n; i++)
        xover[i] = 0;
      zqkp=0.0;
     if (!permute) {
       for (i = 0; i < n; i++)
         w1[i] = (etype) wbis[i] + mu;
       cflo = (etype) c + mu*nbobjets;
     }
     else {
       for (i = 0; i < n; i++)
         w1[i] = (etype) mu*wbis[i] + 1.0;
       cflo = mu*c + (etype) nbobjets;
     }
//       printf("\nmu= %f, c= %d, cflo= %f nbobj= %d end = %d\n", mu, c, cflo, nbobjets, end);

//      cpu_time(&time);
//       fflush(NULL);
//       printf(" time %.2lf\n", time);
//       fflush(NULL);
//      cpu_time(NULL);
      status = quadknap_Sur(n, &cflo, pbis, w1, xover, &zqkp, &stop, nbobjets, c, wbis, &zopt);
      if (zopt > z0)
      {
        z0 = zopt;
        memcpy(xstar, xover, sizeof(int) * MSIZE);
      }

      if(zopt >= zqkp) opt_H_Sur=1;
      count++;
      if (stop == -1) {
        wsum = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) wsum +=wbis[i];
        objsum = 0;
        for (i = 0; i < n; i++)
          if (xover[i]) objsum ++;
//         fflush(NULL);
//         printf("\nzqkp= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsum = %d et wsum = %d et end = %d\n", zqkp, zopt, mu, c, cflo, nbobjets, objsum, wsum, end);
//         fflush(NULL);
        stoppe = 1;
//         fflush(NULL);
//         printf("\nStoppe: zqkp= %d, mu= %f, c= %d, cflo= %f nbobj= %d et end = %d\n", zqkp, mu, c, cflo, nbobjets, end);
//         fflush(NULL);
        break;
      }

      wsum = 0;
      for (i = 0; i < n; i++)
        if (xover[i]) wsum +=wbis[i];
      objsum = 0;
      for (i = 0; i < n; i++)
        if (xover[i]) objsum ++;
//       fflush(NULL);
//       printf("\nzqkp= %d, zopt= %d, mu= %f, c= %d, cflo= %f nbobj= %d et objsum = %d et wsum = %d et end = %d\n", zqkp, zopt, mu, c, cflo, nbobjets, objsum, wsum, end);
//       fflush(NULL);
      if (zqkp < zqkpmin) zqkpmin = zqkp;
      if (zqkp <= z) end = 1;
      if ((wsum <= c) & (objsum <= nbobjets) & !stoppe) {    /*Cas a traiter: wsum <= c et objsum < nbobjets car cont2 en egalite*/
        if (zqkp < zqkpmin) zqkpmin = zqkp;
        end = 1;
      }
      else {
        if (!permute) {
          alpha = (etype) (c - wsum) / (objsum - nbobjets);
          if (objsum > nbobjets) {
            alphal = alpha;
            mul = mu;
          }
          else {
            alphar = alpha;
            mur = mu;
          }
        }
        else {
          alpha = (etype) (nbobjets - objsum) / (wsum - c);
          if (wsum > c) {
            alphal = alpha;
            mul = mu;
          }
          else {
            alphar = alpha;
            mur = mu;
          }
        }
      }
//       printf("alphal = %lf et alphar = %lf et alpha =  %lf et end = %d\n", alphal, alphar, alpha, end);
    }
  } /* end of one iteration */
  if (zqkp < zqkpmin) zqkpmin = zqkp;
  lagbound = zqkpmin;

//   fflush(NULL);
//   printf("Surrogate upper bound %d lower bound %d greedy z %d subgradient nb iteration %d\n", lagbound, z0, z, count);
//   fflush(NULL);
}


/* ======================================================================
        quadknap_EkQKP_Sur
   ====================================================================== */

float quadknap_EkQKP_Sur(int no, int cap, float *ptab, int *wtab, int *xtab, int nbobj, int *H_Sur_opt)//, int z_glob)
{
  int i, j;
  long time, ltime;
  etype ps;
  stype ws;
  ntype nf, fix;
  double times; /* modification LL: ajout de times */
 
  /* initialize global variables */
  n    = no;         /* number of items */
  c    = cap;        /* capacity */
  z    = 0;          /* incumbent solution */
//   z    = z_glob;          /* incumbent solution */
  nf   = 0;          /* number of variables fixed to some value */
  fixp = 0;          /* quadratic profit sum of variables fixed to 1 */
  fixw = 0;          /* weight sum of variables fixed to 1 */
  iterations = 0;    /* number of b&b nodes */
  maxlagr = 10 + n; /* max number of lagrange iterations */
  nbobjets = nbobj; /*modification LL: number of items in the knapqack*/
  ERREUR = 0;
  opt_H_Sur = 0;

  /* copy profits and weights to internal arrays */
  if (no > MSIZE) error("too many items");
  memcpy(pbis, ptab, sizeof(float) * MSIZE * MSIZE);
  memcpy(wbis, wtab, sizeof(int) * MSIZE);

//   printf("n= %d c=%d nbobjets=%d\n",n, c, nbobjets);  
//   for (i = 0; i < n; i++)
//   {
//     printf("i= %d xtab=%d\n",i, xtab[i]);  
//     printf("i= %d w=%d\n",i, wbis[i]);  
//     for (j=0;j<n;j++)
//     {
//       pbis[i][j]=ptab[i][j];
//       if (i==1)
//         printf("i= %d j=%d p=%d\n",i, j, pbis[i][j]);  
//     }
//   }

/*  for (i = 0; i < n; i++) 
  {
    if (x[i])
      printf("quadknap_EkQKP_Sur1: x[%d]=%d\n", i, x[i]);
    if (xstar[i])
      printf("quadknap_EkQKP_Sur1: xstar[%d]=%d\n", i, xstar[i]);
  }*/

  /* find initial solution */
  greedy_EkQKP();
//   fflush(NULL);
//   printf(" greedy EkQKP z %d ", z);
//   fflush(NULL);
  z0 = z;

  if (!ERREUR) sous_gradient(1);
  else  z0 = -1;

//   printf("opt_H_Sur=%d\n",opt_H_Sur);
  *H_Sur_opt = opt_H_Sur;
//   printf("*H_Sur_opt=%d\n",*H_Sur_opt);

  /* copy solution vector to xtab */
  memcpy(xtab, xstar, sizeof(int) * MSIZE);
//   for (i = 0; i < n; i++)
//   {
//     if (xstar[i]) xtab[i]=1;
//     else xtab[i]=0;
// //     printf("xstab[%d] = %d\n",i,xtab[i]);
//   }

  return z0;
}


float quadknap_EkQKP_Sur_bis(int no, int cap, float *ptab, int *wtab, int *xtab, int nbobj, int *H_Sur_opt, double Epsilon)//, int z_glob)
{
  int i, j;
  long time, ltime;
  etype ps;
  stype ws;
  ntype nf, fix;
  double times; /* modification LL: ajout de times */
 
  /* initialize global variables */
  n    = no;         /* number of items */
  c    = cap;        /* capacity */
  z    = 0;          /* incumbent solution */
//   z    = z_glob;          /* incumbent solution */
  nf   = 0;          /* number of variables fixed to some value */
  fixp = 0;          /* quadratic profit sum of variables fixed to 1 */
  fixw = 0;          /* weight sum of variables fixed to 1 */
  iterations = 0;    /* number of b&b nodes */
  maxlagr = 10 + n; /* max number of lagrange iterations */
  nbobjets = nbobj; /*modification LL: number of items in the knapqack*/
  ERREUR = 0;

  /* copy profits and weights to internal arrays */
  if (no > MSIZE) error("too many items");
  memcpy(pbis, ptab, sizeof(float) * MSIZE * MSIZE);
  memcpy(wbis, wtab, sizeof(int) * MSIZE);

//   printf("n= %d c=%d nbobjets=%d\n",n, c, nbobjets);  
//   for (i = 0; i < n; i++)
//   {
//     printf("i= %d xtab=%d\n",i, xtab[i]);  
//     printf("i= %d w=%d\n",i, wbis[i]);  
//     for (j=0;j<n;j++)
//     {
//       pbis[i][j]=ptab[i][j];
//       if (i==1)
//         printf("i= %d j=%d p=%d\n",i, j, pbis[i][j]);  
//     }
//   }

  /* find initial solution */
  greedy_EkQKP();
//   fflush(NULL);
//   printf(" greedy EkQKP z %d ", z);
//   fflush(NULL);
  z0 = z;

//   opt_H_Sur_bis = *H_Sur_opt;
//   printf("*H_Sur_opt_bis=%d\n",*H_Sur_opt);
//   if (!ERREUR && (Epsilon==0.05 || Epsilon==0.00) && !*H_Sur_opt) sous_gradient(1);
//   else  z0 = -1;

//   printf("opt_H_Sur=%d\n",opt_H_Sur);
//   *H_Sur_opt = opt_H_Sur_bis;
//   *H_Sur_opt = 1;
//   printf("*H_Sur_opt=%d\n",*H_Sur_opt);

  /* copy solution vector to xtab */
  memcpy(xtab, xstar, sizeof(int) * MSIZE);
//   for (i = 0; i < n; i++)
//   {
//     if (xstar[i]) xtab[i]=1;
//     else xtab[i]=0;
// //     printf("xstab[%d] = %d\n",i,xtab[i]);
//   }

  return z0;
}


/* ======================================================================
                                  improve
   ====================================================================== */

void improve_EkQKP_end(int *xprime)
{
  int i, j, gaini, gainj;
  etype tot, gain, bgain;
  stype res, kres;
  etype q[MSIZE];
  double eff, maxeff;
  int maxi;
  etype pi;

  res = c; kres = nbobjets;
  // printf("improve: n %d res %d c %d\n", n, res, c); 
  for (i = 0; i < n; i++) {
    if (xprime[i]) {
      res -= w[i];
      kres--;
    }
  }
  for (;;) {
    for (i = 0; i < n; i++) {
      for (tot = p[i][i], j = 0; j < n; j++) {
        if ((j != i) && (xprime[j] != 0)) tot += p[i][j] + p[j][i];
      }
      q[i] = tot;
    }
    bgain = gaini = gainj = 0;
    for (i = 0; i < n; i++) {
      if (xprime[i] == 0) {
        if (w[i] <= res && kres > 0) {
          gain = q[i];
          if (gain > bgain) { bgain = gain; gaini = i; gainj = -1; }
        } else {
          for (j = 0; j < n; j++) {
            if (j == i) continue;
            if (xprime[j] == 0) continue;
            if (w[i] - w[j] <= res) {
              gain = q[i] - q[j] - (p[i][j] + p[j][i]);
              if (gain > bgain) { bgain = gain; gaini = i; gainj = j; }
            }
          }
        }
      }
    }
//     printf("best gain %d i %d j %d res %d\n", bgain, gaini, gainj, res);
    if (bgain == 0) break;
    xprime[gaini] = 1; kres--;
    if (gainj != -1) {xprime[gainj] = 0; kres++;}
    if (gainj != -1) res += w[gainj] - w[gaini]; else res -= w[gaini];
    if (res < 0) {
//       printf("Negative res\n");
//       printf("PROGRAM IS TERMINATED !!!\n\n");
      ERREUR=1;
      return;
    }
//       error("Negative res");
  } // end of while loop 

  for (gain = 0, res = c, kres = nbobjets, i = 0; i < n; i++) {
    if (xprime[i] == 0) continue;
//     printf("1)wi=%d res=%d kres=%d\n",w[i], res, kres);
    res -= w[i]; kres--;
//     printf("2)wi=%d res=%d kres=%d pii=%d\n",w[i], res, kres, p[i][i]);
    for (j = 0; j < n; j++) {
      if (xprime[j]) gain += p[i][j];
    }
  }
  if (res < 0 || kres < 0) {
//       printf("5)Excess capacity\n");
//       printf("PROGRAM IS TERMINATED !!!\n\n");
      ERREUR=1;
      return;
    }
// error("Excess capacity");
  if (kres) {
    for (;;) {
      maxeff = 0; maxi = -1;
      for (i = 0; i < n; i++) {
        if (x[i]) continue;
        pi = -p[i][i];
        for (j = 0; j < n; j++) if (x[j]) pi += p[j][i] + p[i][j];
        eff = pi / (double) w[i];
        if (eff > maxeff && w[i] <= res) { maxeff = eff; maxi = i; }
      }
      if (maxi == -1) {
/*      printf("No item found\n");
      printf("PROGRAM IS TERMINATED !!!\n\n");
 */     ERREUR=1;
      return;
    }
// error("No item found\n");
      i = maxi; x[i] = 1;
      res -= w[i]; kres--;
      if (res >= 0 && kres == 0) break;
    }
  }
  for (gain = 0, res = c, kres = nbobjets, i = 0; i < n; i++) {
    if (xprime[i] == 0) continue;
    res -= w[i]; kres--;
    for (j = 0; j < n; j++) {
      if (xprime[j]) gain += p[i][j];
    }
  }
  if (res < 0 || kres < 0 || kres > 0) {
//       printf("6)Excess capacity\n");
//       printf("PROGRAM IS TERMINATED !!!\n\n");
      ERREUR=1;
      return;
    }
// error("Excess capacity");
//   printf("possible solution %d\n", gain);
  if (gain > z) {
//     printf("improved solution %d to %d\n", z+fixp, gain+fixp);
    z = gain;
    for (i = 0; i < n; i++) xstar[i] = xprime[i];
  }
}


/* ======================================================================
        improve_EkQKP_hybrid
   ====================================================================== */

float improve_EkQKP_hybrid(int no, int cap, float *ptab, int *wtab, int *xtab, int nbobj, float z_glob)
{
  int i, j;
  long time, ltime;
  etype ps;
  stype ws;
  ntype nf, fix;
  double times; /* modification LL: ajout de times */
 
  /* initialize global variables */
  n    = no;         /* number of items */
  c    = cap;        /* capacity */
//   z    = 0;          /* incumbent solution */
  z    = z_glob;          /* incumbent solution */
  nbobjets = nbobj; /*modification LL: number of items in the knapqack*/
  ERREUR = 0;

  /* copy profits and weights to internal arrays */
  if (no > MSIZE) error("too many items");
  memcpy(p, ptab, sizeof(float) * MSIZE * MSIZE);
  memcpy(w, wtab, sizeof(int) * MSIZE);
  memcpy(x, xtab, MSIZE*sizeof(itype));
  memcpy(xstar, xtab, MSIZE*sizeof(boolean));

//   printf("n= %d c=%d nbobjets=%d\n",n, c, nbobjets);  
//   for (i = 0; i < n; i++)
//   {
//     printf("i= %d xtab=%d\n",i, xtab[i]);  
//     printf("i= %d w=%d\n",i, w[i]);  
//     for (j=0;j<n;j++)
//     {
//       p[i][j]=ptab[i][j];
//       if (i==1)
//         printf("i= %d j=%d p=%d\n",i, j, p[i][j]);  
//     }
//   }

  /* find initial solution */
  improve_EkQKP_end(x);
//   fflush(NULL);
//   printf("\nimprove EkQKP z %d\n", z);
//   fflush(NULL);
  z0 = z;

  /* copy solution vector to xtab */
  memcpy(xtab, xstar, sizeof(int) * MSIZE);
//   for (i = 0; i < n; i++)
//   {
//     if (xstar[i]) xtab[i]=1;
//     else xtab[i]=0;
// //     printf("xstab[%d] = %d\n",i,xtab[i]);
//   }

  return z0;
}



/* ======================================================================
                                   cpu_time
   ====================================================================== */

/* The following routine is used for measuring the time used
 *   cpu_time(NULL)    - start of time measuring
 *   cpu_time(&t)      - returns the elapsed cpu-time in variable t (double).
 */

/*void cpu_time(double *t)
{
  static char buffer[30];
  static struct timeval tv;

//   time_t curtime;

  if (t == NULL) {
    gettimeofday(&tv, NULL); 
    t1=tv.tv_sec+(tv.tv_usec/1000000.0);
  } else {
    gettimeofday(&tv, NULL);
    *t = tv.tv_sec+(tv.tv_usec/1000000.0) - t1;
//  static struct tms start, stop;

  if (t == NULL) {
    times(&start);
  } else {
    times(&stop);
    *t = (double) (stop.tms_utime - start.tms_utime) / sysconf(_SC_CLK_TCK);
  }
}*/


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
        printitems
   ====================================================================== */

void printitems(void)
{
  /* print instance to file for further handling */
  int i, j;
  FILE *out;

  out = fopen("save.out", "a");
  fprintf(out, "----------\nsize %d\n", n);
  for (j = 0; j < n; j++) fprintf(out," %3d", w[j]);
  fprintf(out, "\n\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) fprintf(out," %f", p[i][j]);
    fprintf(out, "\n");
  }
  for (j = 0; j < n; j++) fprintf(out," %3d", x[j]);
  fprintf(out, "\nc %d\n", c);
  fclose(out);
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
//   printf("\n\n nbobj: %d\n", nbobj);
//   c    = randm(wsum-50) + 50;
//   c    = randm((nbobj*100)-50) + 50;
  c    = randm((nbobj*30)-50) + 50;
}


/* ======================================================================
                                   sumdata
   ====================================================================== */

void sumdata(int n1, int r1, int pct1, float z, long c, double tottime)
{
  static long n;
  static long r;
  static long pct;
  static double time = 0.0;
  static float ztot   = 0;
  static float zsum   = 0;
  static long csum   = 0;

  if (n1 == 0) {
    fprintf(trace,"n          = %ld\n", n);
    fprintf(trace,"r          = %ld\n", r);
    fprintf(trace,"pct        = %ld\n", pct);
    fprintf(trace,"time       = %.2lf\n", time        / (double) TESTS);
    fprintf(trace,"ztot       = %.1lf\n", ztot        / (double) TESTS);
    fprintf(trace,"zsum       = %.0lf\n", zsum        / (double)     1);
    fprintf(trace,"csum       = %.0lf\n", csum        / (double)     1);
  } else {
    n       = n1;
    r       = r1;
    pct     = pct1;
    ztot   += z;                  /* sum of optimal solutions       */
    time   += tottime;            /* total computational time       */
    zsum    = (float) ((int) (zsum+z) % 1000);  /* controle sum of all solutions  */
    csum    = ((csum+c) % 1000);  /* controle sum of all capacities */
  }
}


/* ======================================================================
        checksolution
   ====================================================================== */

void checksolution(int c, float z)
{
  int i, j;
  float psum;
  int wsum, ksum;

  psum = 0.0;
  wsum = ksum = 0;
  for (i = 0; i < n; i++) {
    if (!xsol[i]) continue;
//    printf("w[%d]=%d \n", i, w[i]);
    wsum += w[i];
    ksum++;
//    printf("i=%d x=%d\n", i, x[i]);
    for (j = 0; j < n; j++) {
       if (xsol[j])
       {
         psum += p[i][j];
//         printf("psum=%f et p[%d][%d]=%f, \n", psum, i, j, p[i][j]);
       }
    }
  }
  printf("CHECK %d: psum %f z %f ksum %d nbobj %d\n", n, psum,z,ksum,nbobj);
// printf("CHECK %d: psum %f z %f wsum %d c %d ksum %d nbobj %d\n", n, psum,z,wsum,c,ksum,nbobj);
// if (wsum > c) terminate_qkp("excess weight");
 if (ksum != nbobj) terminate_qkp("bad nb item");
 if (psum != z) terminate_qkp("bad solution");
}

