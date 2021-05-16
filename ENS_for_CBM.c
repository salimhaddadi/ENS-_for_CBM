/******************************************************************************

This computer code accompanies the paper:

   " Exponential neighborhood search for Consécutive Block Minimization "

by: Salim Haddadi

Submitted to:  International Transactions in Operational Research

May 15, 2021

Note: This code uses the solver "Linkern" which is freely downloadable from the
website:

      http://www.math.uwaterloo.ca/tsp/concorde.html
      
The executable of the solver should be put in a place where it can be accessed

Computing times are provided by the linux command:

$ time -p ./exec

where exec the executable code of this computer code

The running time of Linkern is given by an internal procedure

Results:

- The best configuration is in binary matrix a_best[][]
- The number of 1-blocks in a_best[][] is numbcos_best

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define INFINI  999999999
#define NB_ITER 500
#define RATIO   5
#define M	-1000

unsigned int m,n,numbcos,numbcos_best;
char **a,**a_best;
int iter,nbcol,i,j;
char **b;

void read_data(void);
void fill_file_in_TSPLIB_format(void);
void ens(void);
void Compute_numbcos_best(void);

main()
{
  srand(time(0));
  read_data();
//************** Compute initial tour by using Linkern ************************
  fill_file_in_TSPLIB_format();
  system("linkern -Q -o sol tsp");
  Compute_numbcos_best();
//******************** Initial configuration **********************************
  numbcos = numbcos_best;
  for ( i=0; i<m; i++ )
  for ( j=0; j<=n; j++ ) a[i][j] = a_best[i][j];
//****************** Proceed to the expoential neighborhood search ************
  for ( iter=1; iter<=NB_ITER; iter++ ) ens();
  printf("\nBest number of bco's %d\n\n",numbcos_best);
}

void ens(void)
{
  FILE *ff,*g;
  int i,j,k,e,f,p,kk;
  double r;
  long *numbers,max,nnbcol,nbsub;
  int *col,*perm,**small,*limits,**lenght;
  int *chosen;

  nbcol = n / RATIO;
//*****************************************************************************
  numbcos = numbcos_best;
  for ( i=0; i<m; i++ )
  for ( j=0; j<=n; j++ ) a[i][j] = a_best[i][j];
  nnbcol = 2 * nbcol + 1;
//************ generate n random numbers and find the 2* nbcol smallest *******
  numbers = (long *) malloc ( ( n + 1 ) * sizeof(long) );
  col = (int*) malloc ( ( nnbcol + 1 ) * sizeof(int) );
//*****************************************************************************
  for ( j=1; j<=n; j++ ) numbers[j] = rand();
  for ( i=1; i<=nnbcol; i++ )
  {
    max = 0;
    for ( j=1; j<=n; j++ ) if ( max < numbers[j])
    {
      max = numbers[j];
      f = j;
    }
    numbers[f] = 0;
    col[i] = f;
  }
//*********************** Sort the table col **********************************
  perm = (int*) malloc ( ( nnbcol + 1 ) * sizeof(int) );
//*****************************************************************************
  for(i=1; i<=nnbcol; i++) perm[i] = i;
  for(i=1;i<=nnbcol-1;i++)
  for(j=1;j<=nnbcol-i;j++) if ( col[j] >= col[j+1] )
  {
    f = col[j+1];
    e=perm[j+1];
    col[j+1]=col[j];
    perm[j+1]=perm[j];
    col[j]=f;
    perm[j]=e;
  }
  col[1] = 1;
  col[nnbcol] = n;
//************** Do this to avoid successive values of col[i] *****************
  for ( i=1; i<=nnbcol; i++ )
  {
    f = i/2;
    if ( i - 2 * f == 0 ) col[i] = 0;
  }
  nbcol = 0;
  for ( i=1; i<=nnbcol; i++ )
  {
    if ( col[i] > 0 )
    {
      nbcol++;
      col[nbcol] = col[i];
    }
  }
  nbsub = nbcol-1;
  p = 2 * nbsub;
//********************* Define matrix small of chosen columns *****************
  limits = (int *) malloc ( ( p + 1 ) * sizeof(int) );
  small = (int **) malloc ( m * sizeof(int*) );
  for ( i=0; i<m; i++ ) small[i] = (int *) malloc ( p * sizeof(int) );
//*****************************************************************************
  for ( j=1; j<=nbsub; j++ )
  {
    limits[2*j-1]= col[j] + 1;;
    limits[2*j]= col[j+1];
  }
  limits[1] = 1;
//************************* compute small matrix ******************************
  for ( j=1; j<=p; j++ )
  for ( i=0; i<m; i++ )
  {
    small[i][j]= a[i][limits[j]];
  }
  for ( i=0; i<m; i++ ) small[i][0]= 0;
//************************* compute the lenght matrix *************************
  lenght = (int **) malloc ( ( p + 1 ) * sizeof(int*) );
  for ( i=0; i<=p; i++ ) lenght[i] = (int *) malloc ( ( p + 1 ) * sizeof(int) );
//*****************************************************************************
  for ( i=0; i<=p-1; i++ )
  for ( j=i+1; j<=p; j++ )
  {
    e = 0;
    for ( k=0; k<m; k++ ) e += (1-small[k][i])*small[k][j]+(1-small[k][j])*small[k][i];
    lenght[i][j]=e;
    lenght[j][i] = e;
  }
  for ( i=0; i<=nbsub; i++ ) lenght[i][i] = 0;
  for ( j=1; j<=nbsub; j++ )
  {
    lenght[2*j-1][2*j] = M;
    lenght[2*j][2*j-1] = M;
  }
//***************** Fill file "tsp" for linkern *******************************/
  if ( ( ff = fopen("tsp","w") ) == NULL )
  {
    puts("erreur d'ouverture de fichier");
    exit(1);
  }
  fprintf(ff,"NAME: tsp\n");
  fprintf(ff,"TYPE: TSP\n");
  fprintf(ff,"DIMENSION: %d\n",p + 1);
  fprintf(ff,"EDGE_WEIGHT_TYPE: EXPLICIT\n");
  fprintf(ff,"EDGE_WEIGHT_FORMAT: UPPER_ROW\n");
  fprintf(ff,"EDGE_WEIGHT_SECTION\n");
  for ( i=0; i<=p - 1; i++ )
  {
    for ( j=i+1; j<=p; j++ ) fprintf(ff,"%d ",lenght[i][j]);
    fprintf(ff,"\n");
  }
  fprintf(ff,"EOF");
  fclose(ff);
  system("linkern -Q -otsp.sol tsp >> aux");
//******************** compute the lenght e of the optimal tour ***************
  chosen = (int *) malloc ( ( p + 1 ) * sizeof(int) );
//*****************************************************************************
  if ( ( g = fopen("tsp.sol","r") ) == NULL )
  {
    puts("erreur d'ouverture de fichier");
    exit(1);
  }
  fscanf(g,"%d",&e);
  fscanf(g,"%d",&e);
  e = -nbsub * M;
  fscanf(g,"%d",&kk);
  fscanf(g,"%d",&kk);
  fscanf(g,"%d",&kk);
  e += kk;
  for ( j=1; j<= p; j++ )
  {
    fscanf(g,"%d",&k);
    chosen[j] = k;
    fscanf(g,"%d",&k);
    fscanf(g,"%d",&k);
    e += k;
  }
  fclose(g);
//************************************ Réeorganize columns ********************
  e = 0;
  for ( j=1; j<=nbsub; j++ )
  {
    k = chosen[2*j-1];
    kk = chosen[2*j];
    if ( k < kk )
    {
      for ( f=limits[k]; f<=limits[kk]; f++ )
      {
	e++;
	for ( i=0; i<m; i++ ) b[i][e] = a[i][f];
      }
    }
    else if ( k > kk )
    {
      for ( f=limits[k]; f>=limits[kk]; f-- )
      {
	e++;
	for ( i=0; i<m; i++ ) b[i][e] = a[i][f];
      }
    }
  }
  for ( i=0; i<m; i++ )
  for ( j=1; j<=n; j++ ) a[i][j] = b[i][j];
//*************************** Compute new number of 1-blocks ******************
  numbcos = 0;
  for ( i=0; i<m; i++ )
  {
    for ( j=1; j<=n-1; j++ ) if ( ( a[i][j]==1 ) && ( a[i][j+1] == 0 ) ) numbcos++;
    if ( a[i][n] == 1 ) numbcos++;
  }
  if ( numbcos_best > numbcos )
  {
//    printf("Number of bcos %d    Iter  %d\n",numbcos,iter);
    numbcos_best = numbcos;
    for ( i=0; i<m; i++ )
    for ( j=1; j<=n; j++ ) a_best[i][j] = a[i][j];
  }
}

void read_data(void)
{
  FILE *f;
  char file[20];
  static int e;
  unsigned int k;
  register unsigned int i,j;
  int *th,*h;

  puts("Filename ?");
  scanf("%20s",file);
  if((f=fopen(file,"r"))==NULL)
  {
    puts("erreur d'ouverture du fichier");
    exit(1);
  }
  fscanf(f,"%d%d",&m,&n);
//******************************* Memory allocation ***************************
  a_best = (char **) malloc ( m * sizeof(char*) );
  for ( i=0; i<m; i++ ) a_best[i] = (char *) malloc ( ( n + 1 ) * sizeof(char) );
  a = (char **) malloc ( m * sizeof(char*) );
  for ( i=0; i<m; i++ ) a[i] = (char *) malloc ( ( n + 1 ) * sizeof(char) );
  b = (char **) malloc ( m * sizeof(char*) );
  for ( i=0; i<m; i++ ) b[i] = (char *) malloc ( ( n + 1 ) * sizeof(char) );
  th = (int *) malloc ( (m+1) * sizeof(int) );
//***************** Estimate the number of non null entries *******************
  k = m * n * 0.2;
  h = (int *) malloc ( k * sizeof(int) );
//*****************************************************************************
  k=0;
  for(i=1;i<=m;i++)
  {
    fscanf(f,"%d",&e);
    for(j=k+1;j<=k+e;j++)fscanf(f,"%d",&h[j]);
    k+=e;
    th[i]=k;
  }
  fclose(f);
  th[0]=0;
  for(i=1;i<=m;i++)if( th[i]-th[i-1] == 0 )
  {
    puts("null row");
    printf("i= %d ",i);
    exit(1);
  }
//********************  Recover the binary constraint matrix a ****************
  for(i=0;i<m;i++)
  for(j=0;j<=n;j++) a[i][j]=0;
  for(i=0;i<m;i++)
  for(j=th[i]+1;j<=th[i+1];j++) a[i][h[j]]=1;
  free(h);
  free(th);
}

void fill_file_in_TSPLIB_format(void)
{
  register unsigned int i,j;
  FILE *f;
  char file[20];
  int e,k;
  int **length;

  length = (int **) malloc ( ( n + 2 ) * sizeof(int*) );
  for ( i=0; i<=n+1; i++ ) length[i] = (int *) malloc ( ( n + 2 ) * sizeof(int) );
//****************** Compute the symmetric length matrix **********************
  for(i=0;i<=n;i++)
  for(j=i+1;j<=n;j++)
  {
    e=0;
    for(k=0; k<m; k++) e += (1-a[k][i]) * a[k][j] + a[k][i] * (1-a[k][j]);
    length[i][j] = e;
    length[j][i] = e;
  }
  for ( i=0; i<=n; i++ ) length[i][i] = 0;
//******************** Open and fill file *************************************
  if((f=fopen("tsp","w"))==NULL)
  {
    puts("erreur d'ouverture du fichier");
    exit(1);
  }
  fprintf(f,"NAME: tsp\n");
  fprintf(f,"TYPE: TSP\n");
  fprintf(f,"COMMENT: \n");
  fprintf(f,"DIMENSION: %d\n",n+1);
  fprintf(f,"EDGE_WEIGHT_TYPE: EXPLICIT\n");
  fprintf(f,"EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
  fprintf(f,"EDGE_WEIGHT_SECTION\n");
  for(i=0; i<=n; i++)
  {
    for(j=0; j<=n; j++)fprintf(f," %5d",length[i][j]);
    fprintf(f,"\n");
  }
  fprintf(f,"EOF\n");
  fclose(f);
//************************ Free memory ****************************************
  for ( i=0; i<=n+1; i++ ) free(length[i]);
  free(length);
}

void Compute_numbcos_best(void)
{
  FILE *f;
  static int e,g,h;
  register unsigned int i;
  int *pi;

  pi = (int *) malloc ( ( n + 1 ) * sizeof(int) );
  if((f=fopen("sol","r"))==NULL)
  {
    puts("erreur d'ouverture du fichier");
    exit(1);
  }
  fscanf(f,"%d%d",&e,&g);
  for ( i=0; i<=n; i++ )
  {
    fscanf(f,"%d%d%d",&e,&g,&h);
    pi[i] = e;
  }
  fclose(f);

  for ( i=0; i<m; i++ )
  for ( j=0; j<=n; j++ ) a_best[i][j] = a[i][pi[j]];
  numbcos_best = 0;
  for ( i=0; i<m; i++ )
  {
    for ( j=1; j<=n-1; j++ ) if ( a_best[i][j] == 1 && a_best[i][j+1] == 0 ) numbcos_best++;
    if ( a_best[i][n] == 1 ) numbcos_best++;
  }
  printf("\nInitial number of bcos : %d",numbcos_best);
  free(pi);
}
