 
/* ======================================================================
        QCR  Marie-Christine PLATEAU - Lucas LETOCART
   ====================================================================== */

/* This code is an exact solution method to solve linearly-constrained 
 * zero-one quadratic problems :
 *
 *  (QP) : Minimize   x^t Q x + c^t x + l
 *         s.t.       Ax=b 
 *                    A'x <= b' 
 *                    x \in {0, 1}^n 
 * where :
 *     - Q is a symmetric nxn real matrix, 
 *     - c is an n real vector, 
 *     - l is a real,
 *     - A is an mxn matrix,
 *     - A' is a pxn matrix. 
 *     - b is an m real vector, 
 *     - b' is a p real vector, 
 * 
 *
 * A description of the QCR method is found in: 
 *
 *   A. Billionnet, S. Elloumi, M.-C. Plateau,
 *   "Convex Quadratic Programming for Exact Solution of 0-1 Quadratic
 *    Programs"
 *   
 * The present code is written in ANSI-C, and has been compiled with
 * the GNU-C compiler using option "-ansi -pedantic". It may however
 * be necessary to customize the timing routines for other operating
 * systems. To obtain an exacutable code with the "quadknap" algorithm
 * compile with
 *
 *   gcc -O3 -o QCR_E-kQKP QCR_E-kQKP.c -lm
 *
 * The program reads three arguments
 *
 *   n     number of items 
 *   r     range of coefficients i.e. profits weights are in [1,r]
 *   pct   density of instance, i.e. expected frequency of nonzero
 *         profits
 *
 * having read the parameters, 10 instances are constructed and 
 * solved using the "quadknap" algorithm. Average values of the
 * obtained solutions, capacities, and solution times are printed
 * out to a file named "trace.c".
 *
 * This code can be used free of charge for research and academic purposes
 * only.*/

/*This algorithm solves linearly-constraintes zero-one quadratic problem*/

//#define TESTS       4   /* number of test to run in each series */
#define MSIZE      645   


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <limits.h> 
#include <string.h>
#include <math.h>
#include <sys/times.h>
#include <unistd.h>


/* ======================================================================
				   Macros
   ====================================================================== */

//#define srand(x)     srand48(x)
//#define randm(x)    (lrand48() % (long) (x))


/* ======================================================================
 			    global variables
   ====================================================================== */

static int n, h, o, m, p, nonulA, nonulAbis,hbis;
static float l;
//static  *trace;
static float matq[MSIZE][MSIZE];
static float matc[MSIZE];
static float matl;
static float matA[MSIZE][MSIZE];
static float matAbis[MSIZE][MSIZE]; 
static float matb[MSIZE];
static float matbbis[MSIZE];
FILE *fmatq;
FILE *fmatc;
FILE *fmatl;
FILE *fmatA;
FILE *fmatAbis; 
FILE *fmatb;
FILE *fmatbbis;
FILE *filesdp;

/* ======================================================================
                          reading global variables
   ====================================================================== */

/* The following function is used for reading  global variables where :
 *    - n : is the number of variables
 *    - h : the non-zero coefficients of Q
 *    - o : if o=1, at least one coefficient of vector c is nonnull
 *          else,   all coefficients of c are null.
 *    - l : value of l
 *    - m : number of equality constraints
 *    - p : number of inequality constraintes
 */

void read_data()
{
    int i,ind1,ind2,otemp;
    float val;
    otemp=0;
    fmatq=fopen("q.txt","r");
    fmatc=fopen("c.txt","r");
    fmatl=fopen("l.txt","r");
    fmatA=fopen("A.txt","r");
    fmatAbis=fopen("Abis.txt","r");
    fmatb=fopen("b.txt","r");
    fmatbbis=fopen("bbis.txt","r");  

    fscanf(fmatq,"%d %d",&n,&h);
    fscanf(fmatc,"%d",&o);
    fscanf(fmatl,"%f",&matl);
    fscanf(fmatA,"%d",&m);
//    fscanf(fmatA,"%d",&nonulA);
    fscanf(fmatAbis,"%d",&p);
//   fscanf(fmatAbis,"%d",&nonulAbis);
 


    for(i=1;i<=h;i++)	      //cr tion de la matrice Q
    {
        fscanf(fmatq,"%d %d %f",&ind1,&ind2,&val);
	matq[ind1-1][ind2-1]=val;
// 	printf("%f",matq[ind1-1][ind2-1]);
//	if (ind1-1==ind2-1 && val!= 0){
//	    matc[ind1-1]=val;
//	    otemp=1;
//	}
    }
    
    if(o==1)                  //cr tion de c
	for(i=1;i<=n;i++)
	{
	    fscanf(fmatc,"%f",&val);
	    matc[i-1]=val;
	}
//   o=otemp;
   hbis=h;
 
   if( o==1)                              //s'il existe un vecteur c
	hbis=hbis+n;
//    printf("hbis:%d\n",hbis);
   //      o=otemp;

   //      printf("\notemp:%d\n",otemp);
     if( m != 0)
    {
//	for(i=1;i<=nonulA;i++)	      
	for(i=1;i<=n*m;i++)	      //cr tion de la matrice A des contraintes d' alit 	
{
	    fscanf(fmatA,"%d %d %f",&ind1,&ind2,&val);
	    matA[ind1-1][ind2-1]=val;
	}
	for (i=0;i<m;i++)	      //cr tion de b des contraintes d' alit 	
{
	    fscanf(fmatb,"%f",&val);
	    matb[i]=val;
	}
    }

    if( p != 0)
    {
	for(i=1;i<=n;i++)	//cr tion de la matrice Abis des contraintes d'in.
	{
	    fscanf(fmatAbis,"%d %d %f",&ind1,&ind2,&val);
	    matAbis[ind1-1][ind2-1]=val;
	}
	for (i=0;i<p;i++)       	//cr tion de bbis des contraintes d'in.
	{
	    fscanf(fmatbbis,"%f",&val);
	    matbbis[i]=val;
	}
    }
    

    fclose(fmatq);
    fclose(fmatc);
    fclose(fmatl);
    fclose(fmatA);
    fclose(fmatAbis);
    fclose(fmatb);
    fclose(fmatbbis);

}


void create_file_sb(){
    
    int j,u,t,s,i,ru,r;
    //int hbis=h;
    filesdp=fopen("file.sb","w");
    int nb_cont=n+n*m+m+1+p;
    
    fprintf(filesdp,"%d\n",n+1);
    fprintf(filesdp,"SYMMETRIC_SPARSE\n");
    fprintf(filesdp,"%d %d\n",n+1,hbis);
    
//  for (i=0;i<n-1;i++)				//Ecriture de Q
//	for(j=i+1;j<n;j++)
    for (i=0;i<n;i++)				//Ecriture de Q
	for(j=i;j<n;j++)
	    fprintf(filesdp,"%d %d %f\n",i+1,j+1,-matq[i][j]);
    
    
 if(o==1)                                      //Ecriture de c
     for(i=0;i<n;i++)
     {
       fprintf(filesdp,"0 "); 
       fprintf(filesdp,"%d %f\n", i+1, -(float)matc[i]/2);
     }


fprintf(filesdp,"%d\n", nb_cont);


//CONTRAINTES -bkxi + sum akjXij
for (u=0;u<m;u++)   
   {
       for(t=0;t<n;t++)
       {
	fprintf(filesdp,"LOWRANK_SPARSE_DENSE\n");
	fprintf(filesdp,"%d", n+1 );
	fprintf(filesdp," 1 1\n");
	fprintf(filesdp,"%d", t+1);
	fprintf(filesdp," 0 1.\n");
	fprintf(filesdp,"%d", n+1);
	fprintf(filesdp," 1\n");  
	fprintf(filesdp,"%f\n", -matb[u]);
	for (ru=0;ru<n;ru++)
	    fprintf(filesdp,"%f\n",matA[u][ru]);	
       
	fprintf(filesdp,"= 0\n");
       }
   }


//CONTRAINTES Xii = xi
for (i=1;i<=n;i++)
 {
    fprintf(filesdp,"SYMMETRIC_SPARSE\n");
    fprintf(filesdp,"%d", n+1 );
    fprintf(filesdp," 2\n");
    fprintf(filesdp,"0 ");
    fprintf(filesdp,"%d", i);
    fprintf(filesdp," -1\n");
    fprintf(filesdp,"%d ", i);
    fprintf(filesdp,"%d ", i);
    fprintf(filesdp,"2.\n");  
    fprintf(filesdp,"= 0\n");

 }


//CONTRAINTES D EGALITE *2
 if(m != 0)
     for (s=0;s<m;s++)
     {
       fprintf(filesdp,"SYMMETRIC_SPARSE\n");
       fprintf(filesdp,"%d ", n+1);
       fprintf(filesdp,"%d\n", n);
       for (r=0;r<n;r++)
	{
	    fprintf(filesdp,"0 ");
	    fprintf(filesdp,"%d ",r+1 );
	    fprintf(filesdp,"%f\n",matA[s][r]);
	}
       fprintf(filesdp,"= ");
       fprintf(filesdp,"%f\n",2*matb[s]);
     }


//CONTRAINTES D'INEGALITE *2
 if( p != 0)
     for (s=0;s<p;s++)
     {
	 fprintf(filesdp,"SYMMETRIC_SPARSE\n");
	 fprintf(filesdp,"%d ", n+1);
	 fprintf(filesdp,"%d\n", n);
	 for(r=0;r<n;r++)
	 {
	     fprintf(filesdp,"0 ");
	     fprintf(filesdp,"%d ",r+1 );
	     fprintf(filesdp,"%f\n",matAbis[s][r]);
	 }
	 fprintf(filesdp,"< ");
	 fprintf(filesdp,"%f\n",2*matbbis[s]);
     }


//ECRITURE DE X11 = 1
 fprintf(filesdp,"SINGLETON\n");
 fprintf(filesdp,"%d ", n+1);
 fprintf(filesdp,"0 0 1\n");
 fprintf(filesdp,"= 1\n");
 
 fclose(filesdp);
 }


void create_file_sb_dimplus1(){
    
    int j,u,t,s,i,ru,r,htemp,rubis;
    //int hbis=h;
    htemp=0;
    if (o==1)
      htemp=n;

    filesdp=fopen("file_chr.sb","w");
//     int nb_cont=n+n*m+2*m+1+p;
//     int nb_cont=n+n*m+m+1+p;
    int nb_cont=n+2*m+1+p;
    htemp=htemp+n*(n-1)/2;
    fprintf(filesdp,"%d\n",n+1);
    fprintf(filesdp,"SYMMETRIC_SPARSE\n");
    fprintf(filesdp,"%d %d\n",n+2,htemp);
    
//    for (i=0;i<n-1;i++)				//Ecriture de Q
//	for(j=i+1;j<n;j++)
    for (i=0;i<n-1;i++)				//Ecriture de Q
      for(j=i+1;j<n;j++){
	//	if(-matq[i][j]!=0)
	fprintf(filesdp,"%d %d %f\n",i+1,j+1,-matq[i][j]);
      }
    
 if(o==1)                                      //Ecriture de c
     for(i=0;i<n;i++)
     {
       fprintf(filesdp,"0 "); 
       fprintf(filesdp,"%d %f\n", i+1, -(float)matc[i]/2);
     }


fprintf(filesdp,"%d\n", nb_cont);


//CONTRAINTES -bkxi + sum akjXij
// for (u=0;u<m;u++)   
//    {
//        for(t=0;t<n;t++)
//        {
// 	fprintf(filesdp,"LOWRANK_SPARSE_DENSE\n");
// 	fprintf(filesdp,"%d", n+2 );
// 	fprintf(filesdp," 1 1\n");
// 	fprintf(filesdp,"%d", t+1);
// 	fprintf(filesdp," 0 1.\n");
// 	fprintf(filesdp,"%d", n+2);
// 	fprintf(filesdp," 1\n");  
// 	fprintf(filesdp,"%f\n", -matb[u]);
// 	for (ru=0;ru<n;ru++)
// 	    fprintf(filesdp,"%f\n",matA[u][ru]);	
// 	fprintf(filesdp,"0\n");
// 	fprintf(filesdp,"= 0\n");
//        }
//    }

//CONTRAINTES D'EGALITE AU CARRE
for (u=0;u<m;u++)   
{
   rubis=0;
   for(ru=1;ru<=n;ru++)
     rubis+=ru;    

   fprintf(filesdp,"SYMMETRIC_SPARSE\n");
   fprintf(filesdp,"%d ", n+2);
   fprintf(filesdp,"%d\n", n+rubis);

   for (i=1;i<=n;i++)
   {
      fprintf(filesdp,"0 ");
      fprintf(filesdp,"%d ", i);
      fprintf(filesdp,"%f\n",(1-2*matb[0]));
   }
   for(ru=1;ru<=n;ru++)
   {
      for (rubis=ru;rubis<=n;rubis++)
      {
         if(rubis==ru)
	 {
            fprintf(filesdp,"%d ", ru);
            fprintf(filesdp,"%d ", rubis);
            fprintf(filesdp,"0\n");
//             fprintf(filesdp,"1\n");
	 }
	 else
	 {
            fprintf(filesdp,"%d ", ru);
            fprintf(filesdp,"%d ", rubis);
//             fprintf(filesdp,"1\n");
	    fprintf(filesdp,"2.\n");	
// 	    fprintf(filesdp,"%f\n",2*matA[u][rubis-1]);	
	 }
      }
   }
//    fprintf(filesdp,"= 0\n");
   fprintf(filesdp,"= %f\n", -2*matb[0]*matb[0]);
}

//CONTRAINTES Xii = xi
for (i=1;i<=n;i++)
 {
    fprintf(filesdp,"SYMMETRIC_SPARSE\n");
    fprintf(filesdp,"%d", n+2 );
    fprintf(filesdp," 2\n");
    fprintf(filesdp,"0 ");
    fprintf(filesdp,"%d", i);
    fprintf(filesdp," -1\n");
    fprintf(filesdp,"%d ", i);
    fprintf(filesdp,"%d ", i);
    fprintf(filesdp,"2.\n");  
    fprintf(filesdp,"= 0\n");

 }


//CONTRAINTES D EGALITE *2
 if(m != 0)
     for (s=0;s<m;s++)
     {
       fprintf(filesdp,"SYMMETRIC_SPARSE\n");
       fprintf(filesdp,"%d ", n+2);
       fprintf(filesdp,"%d\n", n);
       for (r=0;r<n;r++)
	{
	    fprintf(filesdp,"0 ");
	    fprintf(filesdp,"%d ",r+1 );
	    fprintf(filesdp,"%f\n",matA[s][r]);
	}
       fprintf(filesdp,"= ");
       fprintf(filesdp,"%f\n",2*matb[s]);
     }


//CONTRAINTES D'INEGALITE *2
 if( p != 0)
     for (s=0;s<p;s++)
     {
	 fprintf(filesdp,"SYMMETRIC_SPARSE\n");
	 fprintf(filesdp,"%d ", n+2);
	 fprintf(filesdp,"%d\n", n);
	 for(r=0;r<n;r++)
	 {
	     fprintf(filesdp,"0 ");
	     fprintf(filesdp,"%d ",r+1 );
	     fprintf(filesdp,"%f\n",matAbis[s][r]);
	 }
	 fprintf(filesdp,"< ");
	 fprintf(filesdp,"%f\n",2*matbbis[s]);
     }


//ECRITURE DE X11 = 1
 fprintf(filesdp,"SINGLETON\n");
 fprintf(filesdp,"%d ", n+2);
 fprintf(filesdp,"0 0 1\n");
 fprintf(filesdp,"= 1\n");
 
 fclose(filesdp);
 }

void transforme_prob_au()
{
    int i;
    float val;
    FILE *fprob,*fau;
    fprob=fopen("prob.sol","r");
    fau=fopen("au.txt","w");
    int nau=n*m+n+m+1+p;
    fprintf(fau,"%d %d\n",nau,1);
    for(i=1;i<=nau-1;i++)
    {
	fscanf(fprob,"%f",&val);
        fprintf(fau,"%f\n",val); 
    }
    fscanf(fprob,"%f",&val);
    fprintf(fau,"%f ",val); 

    fclose(fprob);
    fclose(fau);
}



void create_newQ_newc_newA (){

    int i,j,logemax, htemp;
    int emax;
    FILE *femax;
    femax=fopen("emax.txt","r");
    fscanf(femax,"%d",&emax);

    if (emax != 0)
	logemax=ceil(log((float)emax)/log(2)); //calcul du nbre de variables ajout s
    else
	logemax=0;
    if (emax==1)
      logemax=1;

    printf("emax:%d\n",emax);
   printf("logemax:%d\n",logemax);

    if (emax!=0)
    {
	htemp=((n+logemax)*(n+logemax-1))/2;

	fmatq=fopen("q.txt","w");
	fprintf(fmatq,"%d ", n+logemax);		// riture dans q.txt
	fprintf(fmatq,"%d\n", htemp);
	for (i=1;i<=(n+logemax-1);i++)
	    for (j=i+1;j<=n+logemax;j++)
	    {
		if (i==n+logemax-1 && j==n+logemax)
		    fprintf(fmatq,"%d %d %f ", i, j, matq[i-1][j-1]);
		else
		    fprintf(fmatq,"%d %d %f\n", i,j,matq[i-1][j-1]);
	    }
	
	
	if (o==1)
	{
	    fmatc=fopen("c.txt","w");
	    fprintf(fmatc,"%d\n", 1);			// riture dans c.txt
	    for (i=1;i<=(n+logemax);i++)
	    {
		if (i==(n+logemax))
		    fprintf(fmatc,"%f ", matc[i-1]);
		else
		    fprintf(fmatc,"%f\n", matc[i-1]);		
	    }
	}	     
	
	fmatA=fopen("A.txt","w");
	if (m==1)
	{
	    fmatA=fopen("A.txt","w");
	    fprintf(fmatA,"%i\n", 2);			// riture dans A.txt
	    for (i=1;i<=n;i++)
               fprintf(fmatA,"%d %d %f\n", 1, i, matA[0][i-1]);	
	    j=0;
	    	  for (i=(n+1);i<=(n+logemax);i++)
	    {
	    
	        fprintf(fmatA,"%d %d %d\n", 1, i, 0);
	      j++;	
	    }
	}

	
	    for (i=1;i<=n;i++)
               fprintf(fmatA,"%d %d %f\n", 2, i, matAbis[0][i-1]);	
	    j=0;
	    for (i=(n+1);i<=(n+logemax);i++)
	    {
		if (i==(n+logemax))
		    fprintf(fmatA,"%d %d %f ", 2, i, pow(2,j));		
		else
		    fprintf(fmatA,"%d %d %f\n", 2, i, pow(2,j));
	        j++;	
	    }


	fmatb=fopen("b.txt","w");
	fprintf(fmatb,"%f\n",matb[0]);
	fprintf(fmatb,"%f ",matbbis[0]);


	fclose(fmatq);
	fclose(fmatb);
	fclose(fmatc);
	fclose(fmatA);
	
    }
    else
    {
	fmatA=fopen("A.txt","w");
	if (m==1)
	{
	    fmatA=fopen("A.txt","w");
	    fprintf(fmatA,"%i\n", 2);			// riture dans A.txt
	    for (i=1;i<=n;i++)
               fprintf(fmatA,"%d %d %f\n", 1, i, matA[0][i-1]);	
	    j=0;
	    	  for (i=(n+1);i<=(n+logemax);i++)
	    {
	    
	        fprintf(fmatA,"%d %d %d\n", 1, i, 0);
	      j++;	
	    }
	}

	
	    for (i=1;i<=n;i++)
               fprintf(fmatA,"%d %d %f\n", 2, i, matAbis[0][i-1]);	
	    j=0;
	    for (i=(n+1);i<=(n+logemax);i++)
	    {
		if (i==(n+logemax))
		    fprintf(fmatA,"%d %d %f ", 2, i, pow(2,j));		
		else
		    fprintf(fmatA,"%d %d %f\n", 2, i, pow(2,j));
	        j++;	
	    }


	fmatb=fopen("b.txt","w");
	fprintf(fmatb,"%f\n",matb[0]);
	fprintf(fmatb,"%f ",matbbis[0]);


	fclose(fmatA);
	fclose(fmatb);
    }


}

int  main(){
    
     read_data();
     
//     create_newQ_newc_newA ();
//     system("cp Abis0.txt Abis.txt");
//     system("cp bbis0.txt bbis.txt");
//     read_data();
    create_file_sb_dimplus1();
    
     //create_file_sb();
     //system("sb -oy au.txt -f file.sb");	
    // printf("Avant SB2CSDP\n");
     system("./SB2CSDP file_chr.sb");
    // printf("Apres SB2CSDP\n");
//      system("../SB2CSDP file.sb");
//      system("../csdp file.sdpa prob.sol > sol.txt");
    // printf("Avant csdp\n");
    system("./csdp file_chr.sdpa prob.sol > sol_chr.txt");
    // printf("Apres csdp\n");
//      system("time (../csdp file.sdpa prob.sol > sol.txt)");
//      transforme_prob_au();
//      system("ampl -f ./pg.run");

return 0;
}
