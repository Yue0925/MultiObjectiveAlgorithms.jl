#include <stdio.h>
#include <stdlib.h>
#include <string.h>

FILE *SB;// Fichier SB
int n,m,nb,tmp; //n=dimension, m=nb de contraintes, nb=nombre d'elements non nul dans C
int nbEQ, nbINE; // Nombre de contraintes d'egalit  d'in alit double *c;
double *c;
float **C;
float ***A;
int *NBZA; // NBZA[i]= nombre d'elements non nuls dans A[i]
char *TYPE; // TYPE[i] est dans {'=','<','>'}

void lectureFichier()
{
  int i,j,l,k;
  int nlB,ncB,nlD,ncD;
  int N;
  int compteur;
  float val;
  int ind1,ind2;
  char * type;
  char tmpCHAR;
  int dummy;
  char * dumc;

  float ** B;
  float ** D;
  float ** BD;

  dumc = (char *)malloc(sizeof(char)*500);
  type = (char *)malloc(500*sizeof(char)); // pour reconnaitre le type de contrainte
  nbEQ = 0;
  nbINE= 0;

  // Trace constante: on oublie

  fscanf(SB,"%d",&dummy);
  fscanf(SB,"%s",dumc);

  // Dimension du probleme
  fscanf(SB,"%d %d",&n,&nb); //nb=nombre d'elements non nuls dans C

//   printf("nb non nuls dans C = %d\n",n);

  // Matrice definissant la fonction objectif : max C.X
  C = (float **)malloc(3*sizeof(float *));
  C[0] = (float *)malloc(nb*sizeof(float));  // premier indice
  C[1] = (float *)malloc(nb*sizeof(float));  // deuxieme indice
  C[2] = (float *)malloc(nb*sizeof(float));  // valeur

  // LECTURE DE LA MATRICE DEFINISSANT LA FONCTION: CX
  for(i=0;i<nb;i++)
    {
      fscanf(SB,"%f %f %f\n",&(C[0][i]),&(C[1][i]),&(C[2][i]));
    }

  // A contiendra les matrices definissant les contraintes Ai.X=ci

  // Nombre de contraintes
  fscanf(SB,"%d",&m);
//   printf("nb contraintes = %d\n",m);

  A =  (float ***)malloc(sizeof(float **)*m); //A[i]: ieme matrice
  NBZA =  (int *)malloc(sizeof(int)*m); // NBZA[i]: nombre d'elements non nuls dans A[i]
  TYPE = (char *)malloc(sizeof(char)*m);  // TYPE[i] est dans {'=','<','>'}
  c =  (double *)malloc(sizeof(double)*m); //c[i]: second membre
  B =  (float **)malloc(sizeof(float *)*n); // B et D sont utilisees
  for(i=0;i<n;i++) B[i] =  (float *)malloc(sizeof(float)*n);
  D =  (float **)malloc(sizeof(float *)*n); // pour LOWRANK_SPARSE_SPARSE
  for(i=0;i<n;i++) D[i] =  (float *)malloc(sizeof(float)*n);
  BD =  (float **)malloc(sizeof(float *)*n);
  for(i=0;i<n;i++) BD[i] =  (float *)malloc(sizeof(float)*n);  // BD=BD^T+DB^T

  // TRAITEMENT DES CONTRAINTES DANS SB
  for(i=0;i<m;i++)
    {
      // Type de contrainte: SINGLETON, LOWRANK_SPARSE_SPARSE, SYMMETRIC_SPARSE, LOWRANK_SPARSE_DENSE
      fscanf(SB,"%s",type);

      /****************************************
       *        SINGLETON                       *
       ****************************************/
      if (type[0]=='S' && type[1]=='I')
	{
	  A[i] =  (float **)malloc(sizeof(float *));
	  NBZA[i] = 1;
	  A[i][0] =  (float *)malloc(3*sizeof(float));
	  fscanf(SB,"%d",&dummy); // dimension: on oublie
	  fscanf(SB,"%f %f %f\n",&(A[i][0][0]),&(A[i][0][1]),&(A[i][0][2]));
	  // contrainte d'egalite ou d'inegalite ?
	  fscanf(SB,"%s",type);
	  tmpCHAR = type[0];
	  if (strlen(type) == 1)
	    {
	      fscanf(SB,"%s",type);
	      c[i] = atof(type);
	    }
	  else c[i] = atof(type+1);
	  if (tmpCHAR=='=')
	    {
	      nbEQ++;
	      TYPE[i] = '=';
	    }
	  if (tmpCHAR=='<')
	    {
	      nbINE++;
	      TYPE[i] = '<';
	    }
	  if (tmpCHAR=='>')
	    {
	      nbINE++;
	      TYPE[i] = '>';
	    }
	}
      /****************************************
       *        SYMMETRIC_SPARSE                *
       *****************************************/
      if (type[0]=='S' && type[1]=='Y')
	{
	  // ieme contrainte
	  fscanf(SB,"%d",&dummy); // dimension: on oublie
	  fscanf(SB,"%d",&N); // Nombre d'elements non nuls
	  A[i] =  (float **)malloc(sizeof(float *)*N);
	  NBZA[i] = N;
	  for(k=0;k<N;k++)
	    {
	      A[i][k] =  (float *)malloc(3*sizeof(float));
	      fscanf(SB,"%f %f %f",&(A[i][k][0]),&(A[i][k][1]),&(A[i][k][2]));
	    }
	  fscanf(SB,"%s",type);
	  tmpCHAR = type[0];
	  if (strlen(type) == 1)
	    {
	      fscanf(SB,"%s",type);
	      c[i] = atof(type);
	    }
	  else c[i] = atof(type+1);
	  // contrainte d'egalite ou d'inegalite ?
	  if (tmpCHAR=='=')
	    {
	      nbEQ++;
	      TYPE[i] = '=';
	    }
	  if (tmpCHAR=='<')
	    {
	      nbINE++;
	      TYPE[i] = '<';
	    }
	  if (tmpCHAR=='>')
	    {
	      nbINE++;
	      TYPE[i] = '>';
	    }
	}
      /****************************************
       *        LOWRANK_SPARSE_SPARSE           *
       *****************************************/
      if (type[0]=='L' && type[15] =='S')
	{
	  // RAZ des matrices B et D
	  for(j=0;j<n;j++)
	    {
	      for(k=0;k<n;k++)
		{
		  B[j][k] = 0;
		  D[j][k] = 0;
		  BD[j][k] = 0;
		}
	    }
	  // Lecture de la premiere matrice:B
	  fscanf(SB,"%d",&nlB); // nombre de lignes de B
	  fscanf(SB,"%d",&ncB); // nombre de colonnes de B
	  fscanf(SB,"%d",&N); // Nombre d'elements non nuls
	  for(k=0;k<N;k++)
	    {
	      fscanf(SB,"%d %d %f",&ind1,&ind2,&val);
	      B[ind1][ind2] = val;
	    }
	  // Lecture de la deuxieme matrice:D
	  fscanf(SB,"%d",&nlD); // nombre de lignes de D
	  fscanf(SB,"%d",&ncD); // nombre de colonnes de D
	  fscanf(SB,"%d",&N); // Nombre d'elements non nuls
	  for(k=0;k<N;k++)
	    {
	      fscanf(SB,"%d %d %f",&ind1,&ind2,&val);
	      D[ind1][ind2] = val;
	    }
	  // Calcul de BD^T+DB^T
	  NBZA[i] = 0;
	  for(j=0;j<nlB;j++)
	    {
	      for(k=0;k<=j;k++) // matrice symetrique
		{
		  for(l=0;l<ncD;l++) BD[j][k] += B[j][l]*D[k][l]+B[k][l]*D[j][l];
		  if (BD[j][k]) NBZA[i]++;
		}
	    }
	  // Stockage de la contrainte dans A[i] qui contient NBZA[i] elements non nuls
	  A[i] =  (float **)malloc(sizeof(float *)*NBZA[i]);
	  compteur = 0;
	  // Parcours de BD
	  for(j=0;j<nlB;j++)
	    {
	      for(k=0;k<=j;k++)
		{
		  if(BD[j][k])
		    {
		      A[i][compteur] =  (float *)malloc(3*sizeof(float));
		      A[i][compteur][0]=k;
		      A[i][compteur][1]=j;
		      A[i][compteur][2]=BD[j][k];
		      compteur++;
		    }
		}
	    }

	  // Contrainte d'egalite ou d'inegalite ?
	  fscanf(SB,"%s",type);
	  tmpCHAR = type[0];
	  if (strlen(type) == 1)
	    {
	      fscanf(SB,"%s",type);
	      c[i] = atof(type);
	    }
	  else c[i] = atof(type+1);
	  if (tmpCHAR=='=')
	    {
	      nbEQ++;
	      TYPE[i] = '=';
	    }
	  if (tmpCHAR=='<')
	    {
	      nbINE++;
	      TYPE[i] = '<';
	    }
	  if (tmpCHAR=='>')
	    {
	      nbINE++;
	      TYPE[i] = '>';
	    }
	} // Fin de LOWRANK_SPARSE_SPARSE

      /****************************************
       *        LOWRANK_SPARSE_DENSE            *
       *****************************************/
      if (type[0]=='L' && type[15] =='D')
	{
	  // RAZ des matrices B et D
	  for(j=0;j<n;j++)
	    {
	      for(k=0;k<n;k++)
		{
		  B[j][k] = 0;
		  D[j][k] = 0;
		  BD[j][k] = 0;
		}
	    }
	  // Lecture de la premiere matrice:B
	  fscanf(SB,"%d",&nlB); // nombre de lignes de B
	  fscanf(SB,"%d",&ncB); // nombre de colonnes de B
	  fscanf(SB,"%d",&N); // Nombre d'elements non nuls
	  for(k=0;k<N;k++)
	    {
	      fscanf(SB,"%d %d %f",&ind1,&ind2,&val);
	      B[ind1][ind2] = val;
	    }

	  // Lecture de la deuxieme matrice:D
	  fscanf(SB,"%d",&nlD); // nombre de lignes de D
	  fscanf(SB,"%d",&ncD); // nombre de colonnes de D
	  if(nlB != nlD || ncB != ncD)
	    {
	      printf("\nerreur dans les dimensions des matrices sparse_dense\n");
	      exit(1);
	    }
	  for(j = 0; j < nlD; j++)
	    {
	      for(k = 0; k < ncD; k++)
		{
		  fscanf(SB, "%f", &val);
		  D[j][k] = val;
		}
	    }

	  // Calcul de BD^T+DB^T
	  NBZA[i] = 0;
	  for(j=0;j<nlB;j++)
	    {
	      for(k=0;k<=j;k++) // matrice symetrique
		{
		  for(l=0;l<ncD;l++) BD[j][k] += B[j][l]*D[k][l]+B[k][l]*D[j][l];
		  if (BD[j][k]) NBZA[i]++;
		}
	    }
	  // Stockage de la contrainte dans A[i] qui contient NBZA[i] elements non nuls
	  A[i] =  (float **)malloc(sizeof(float *)*NBZA[i]);
	  compteur = 0;
	  // Parcours de BD
	  for(j=0;j<nlB;j++)
	    {
	      for(k=0;k<=j;k++)
		{
		  if(BD[j][k])
		    {
		      A[i][compteur] =  (float *)malloc(3*sizeof(float));
		      A[i][compteur][0]=k;
		      A[i][compteur][1]=j;
		      A[i][compteur][2]=BD[j][k];
		      compteur++;
		    }
		}
	    }

	  // Contrainte d'egalite ou d'inegalite ?
	  fscanf(SB,"%s",type);
	  tmpCHAR = type[0];
	  if (strlen(type) == 1)
	    {
	      fscanf(SB,"%s",type);
	      c[i] = atof(type);
	    }
	  else c[i] = atof(type+1);
	  if (tmpCHAR=='=')
	    {
	      nbEQ++;
	      TYPE[i] = '=';
	    }
	  if (tmpCHAR=='<')
	    {
	      nbINE++;
	      TYPE[i] = '<';
	    }
	  if (tmpCHAR=='>')
	    {
	      nbINE++;
	      TYPE[i] = '>';
	    }
	} // Fin de LOWRANK_SPARSE_DENSE


    } // Fin de la boucle en i
//   printf("nb EQ = %d, nb INEQ = %d\n",nbEQ,nbINE);
}

void ecriturePb(char *filename)
{
  FILE *fSDPA;
  int i,k;
  int NBecart; // indice de la variable d'ecart courante

  NBecart = 1;
  fSDPA = fopen(filename,"w");
  if(fSDPA == NULL)
    {
      printf("Error while opening %s.\n",filename);
      exit(-1);
    }

  fprintf(fSDPA,"%d =mdim\n",m);
  if (nbINE)
    {
      fprintf(fSDPA,"2 =nblocks\n");
      fprintf(fSDPA,"{%d, -%d}\n",n,nbINE);
    }
  else
    {
      fprintf(fSDPA,"1 =nblocks\n");
      fprintf(fSDPA,"{%d}\n",n);
    }
  // ECRITURE de c
  for(i=0;i<m;i++)
    {
      fprintf(fSDPA,"%f ",c[i]);
      //      printf("%d : %f\n",i,c[i]);
      //      if (!(i % 10)) fprintf(fSDPA,"\n");
    }
  fprintf(fSDPA,"\n");
  // ECRITURE DE C
  for(i=0;i<nb;i++)
    {
      // MATRICE O, BLOCK 1, i, j, val
      fprintf(fSDPA,"0 1 %d %d %f\n", (int)(C[0][i])+1,(int)(C[1][i])+1,C[2][i]);
    }

  // ECRITURE DES CONTRAINTES (STOCKEES DANS A)
  for(i=0;i<m;i++)
    {
      // CONTRAINTE D'EGALITE
      if(TYPE[i]=='=')
	{
	  // MATRICE i, BLOCK 1, ik, jk, valk
	  for(k=0;k<NBZA[i];k++)
	    fprintf(fSDPA,"%d 1 %d %d %f\n", i+1, (int)A[i][k][0]+1,(int)A[i][k][1]+1,A[i][k][2]);
	}
      // CONTRAINTE D'INEGALITE '<'
      if(TYPE[i]=='<')
	{
	  // MATRICE i, BLOCK 1, ik, jk, valk
	  for(k=0;k<NBZA[i];k++)
	    fprintf(fSDPA,"%d 1 %d %d %f\n", i+1, (int)A[i][k][0]+1,(int)A[i][k][1]+1,A[i][k][2]);
	  // VARIABLE D'ECART: MATRICE i, BLOCK 2, NBecart, NBecart, 1.0
	  fprintf(fSDPA,"%d 2 %d %d 1.0\n", i+1, NBecart, NBecart);
	  NBecart++;
	}
      // CONTRAINTE D'INEGALITE '>'
      if(TYPE[i]=='>')
	{
	  // MATRICE i, BLOCK 1, ik, jk, valk
	  for(k=0;k<NBZA[i];k++)
	    fprintf(fSDPA,"%d 1 %d %d %f\n", i+1, (int)A[i][k][0]+1,(int)A[i][k][1]+1,A[i][k][2]);
	  // VARIABLE D'ECART: MATRICE i, BLOCK 2, NBecart, NBecart, -1.0
	  fprintf(fSDPA,"%d 2 %d %d -1.0\n", i+1, NBecart, NBecart);
	  NBecart++;
	}
    }
}

int main(int argc,char **argv)
{
  char buf[100];

  if(argc != 2)
    {
      printf("\n\t Usage : SB2CSDP <filename.sb>\n");
      exit(-1);
    }

  SB = fopen(argv[1],"r");
  if(SB == NULL)
    {
      printf("Error while opening %s\n",argv[1]);
      exit(-2);
    }

  lectureFichier();
  sprintf(buf,"%s",argv[1]);
  sprintf(buf+(strlen(argv[1])-3),".sdpa");
  ecriturePb(buf);
//   printf("File .sdpa written.\n\n");
  return 0;
}

