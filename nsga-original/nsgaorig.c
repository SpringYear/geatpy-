/********************************************************************\
***                      NON-DOMINATED SORTING                     ***
***                             using                              ***
***                       GENETIC ALGORITHM                        ***
***                              for                               ***
***                  MULTI OBJECTIVE OPTIMIZATION                  ***
***                                                                ***
***        Developed by : Kalyanmoy Deb and Mayank Goyal           ***
***           The Department of Mechanical Engineering             ***
***            Indian Institute of Technology, Kanpur              ***
***                  Kanpur, PIN 208 016, INDIA                    ***
***................................................................***
*** This is a GA implementation for multi-objective optimization.  ***
*** For multi-objective optimization, non-domonated sorting has    ***
*** been used.  The design variables can be Binary, Integer, Real  ***
*** or Enumerated data type.  Moreover a design vector can have    ***
*** design variables where each variable can be any of the above   ***
*** mentioned types.  The order in which the design variables      ***
*** are needed is also maintained in the code.  For multi-objective***
*** optimization, sharing is done using either of two ways :       ***
*** Sharing on Fitness Space or Sharing on Parameter Space. After  ***
*** completion, results are stored in a 'result.out' and for       ***
*** detailed inspection in a file 'report'. For a detail discussion***
*** please refer to Journal paper:                                 ***
*** Srinivas, N. and Deb, K. (1994). Multiobjective optimization   ***
***  using nondominated sorting genetic algorithms. Evolutionary   ***
***  Computation. Vol. 2, No. 3. Pages 221--248.                   ***
***                                                                ***
*** Please send your comments or bug information at                ***
*** deb@ls11.informatik.uni-dortmund.de or deb@iitk.ernet.in       ***
***                                                                ***
*** All rights reserved. Not to be used for commercial purposes.   ***
*** In case of a publication resulted by use of this code,         ***
*** an acknowledgement will be appreciated.                        ***
***                                                                ***
*** Last update 15 June 1998       Kalyanmoy Deb                   ***
\********************************************************************/

#include<stdio.h>
#include<math.h>

#define BITS_PER_BYTE 8
#define UINTSIZE      (BITS_PER_BYTE*sizeof(unsigned))
#define INFINITY      1e7
#define EPSILON       1e-6
#define PI            3.1415927
#define MAXVECSIZE    30 
#define MAXPOPSIZE    500
#define MAXENUMDATA   100
#define MAXOBJ        5
#define MINSCHEME     1   /* do not change */
#define TRUE          1   /*     ,,        */
#define FALSE         0   /*     ,,        */
#define ONESITE       1   /*     ,,        */
#define UNIF          2   /*     ,,        */
#define BIN           1   /*     ,,        */
#define INT           2   /*     ,,        */
#define ENUM          3   /*     ,,        */
#define CONT          4   /*     ,,        */
#define square(x)  ((x)*(x)) 
#define cube(x)  ((x)*(x)*(x)) 
#define PENALTY_COEFF 1.0e3

/*=============================
  Choose your problem here (put your problem at the end of code) 
  ===========================*/
#define book

/*=================
  TYPE DEFINTIONS :
  =================*/
struct indiv      
{  
  unsigned *chrom;		/* chrosome string 	*/
  float fitness[MAXOBJ];	/* fitness functions    */
  float x[MAXVECSIZE];	        /* variables  		*/
  float dumfitness;	        /* modified objective   */ 
  int   flag;                   /* flag as follows
			            0 => untouched
			            1 => dominated
		                    2 => non-dominated
		        	    3 => exhausted   */
  int   front;                  /* Front of the indiv.  */
  int   parent1;                /* two parents */
  int   parent2;	
};
typedef struct indiv INDIVIDUAL ;
typedef INDIVIDUAL *POPULATION ;	/* array of individuals */

/*====================
  FUNCTION PROTOTYPES :
  ====================*/
float 	randomperc(); 
float 	get_beta();
float	get_delta();
float   get_beta_rigid();
double  noise();
float 	rndreal();
int   	small(float);

/*==================
  GLOBAL VARIABLES  :
  ==================*/
int   	pop_size,   		/* Population Size  		 */
  gen_no,			/* Current generation number	 */
  max_gen, 		/* Maximum no. of generations    */
  no_xover, 		/* No. of cross overs done 	 */
  no_mutation,		/* No. of mutations done 	 */
  num_var,		/* Number of total design variables */
  num_bin_var,            /* Number of binary variables    */
  num_int_var,            /* Number of integer variables   */
  num_enum_var,           /* Number of enumerated variables */
  num_cont_var,           /* Number of continuous variables */
  num_obj,                /* Number of objectives           */
  lchrom[MAXVECSIZE], /* Length of chromosome		  */
  chromsize,          /* Number of bytes needed to store 
			 lchrom strings */
  x_strategy, 	    /* Crossover strategy UNIF,ONESITE etc. */
  REPORT,		    /* Flag for Full reports (True/False) */
  var_RIGID[MAXVECSIZE],  /* Flag for rigid boundaries (T/F) */
  BINGA,			/* Flag for binary GA (T/F)	*/
  REALGA,			/* Flag for real-GA (T/F)	  */
  FITNESS,                /* Flag for sharing strategy (T/F)  */
  PARAM,
  minmax[MAXOBJ],          /* min or max flag */
  tourneylist[MAXPOPSIZE],/* List of indices of individuals for
			     tournament selection routine */
  tourneypos, 		/* Current position of tournament   */
  tourneysize,		/* Tournament size ( = 2 for binary )*/
  var_type[MAXVECSIZE],   /* Temp. variable type            */
  TotalChromLength,       /* Total Chromosome Length       */
  num_enum_data[MAXVECSIZE];/* Number of Enumerated Data for
			       each enumerated variable    */

float 	seed,			/* Random seed number	   	*/
  n_distribution_c, n_distribution_m, /* Distribution
					 index for SBX 	*/
  p_xover, 		/* Cross over probability	*/
  p_mutation, 		/* Mutation probability		*/
  closeness,		/* Closeness epsilon		*/
  minx[MAXVECSIZE],	/* Minimum value of design variables */
  maxx[MAXVECSIZE],   	/* Maximum value of design variables */
  x_lower[MAXVECSIZE], 	/* Lower and Upper bounds on each  */
  x_upper[MAXVECSIZE], 	/*        design variable	   */
  afmin[MAXOBJ],          /* approx min value of obj funs. */
  afmax[MAXOBJ],          /* approx max value of obj funs  */
  dshare,	                /* Sharing distance		   */
  max_spread, 		/* Maximum spread */
  maxf[MAXOBJ],           /* Maximum fitness values         */
  minf[MAXOBJ],           /* Minimum fitness values         */
  avgfitness[MAXOBJ],     /* Average fitness values         */
  dum_avg, min_dum,       /* Avg. and min. dummy fitness      */
  init_dum,delta_dum,     /* Initial and delta dummy fitness  */
  c1=0.0,c2=0.0,          /* Children                         */
  weightage[MAXOBJ],      /* Weightage assigned to fitness    */
  EnumData[MAXVECSIZE][MAXENUMDATA]; /* Enumerated Data       */

POPULATION oldpop, newpop;	/* Old and New populations   */

int *choices, nremain;
float *fraction;

/*====================================================================
  SUBROUTINE FOR INPUTTING GLOBAL PARAMETERS :
  ====================================================================*/
input_parameters()
{
  int   k,temp2,count;
  char  ans;
  float temp1, sumweight;
  
  printf("\nNumber of objective functions (2) --> ");
  scanf("%d",&num_obj);
  printf("\nSpecify minimization (1) or maximization (-1) for each function ");
 for (count=0; count<num_obj; count++)
    {
      printf("\n  Objective function #%2d (1 or -1)      --> ",count+1); 
      scanf("%d",&minmax[count]);
    }
  
  printf("\nNumber of variables (Maximum %2d) --- : ",MAXVECSIZE);
  scanf("%d",&num_var);
  
  BINGA = FALSE;
  REALGA = FALSE;
  TotalChromLength=0;
  
  num_bin_var=0; num_int_var=0; num_enum_var=0; num_cont_var=0;
  for (k=0; k<= num_var-1; k++)
    {
      printf("\nVariable Type for variable #%2d : \n",k+1);
      printf("\t\t 1 for Binary\n");
      printf("\t\t 2 for Integer\n");
      printf("\t\t 3 for Enumerated\n");
      printf("\t\t 4 for Continuous/Real");
      printf("\n  Give your choice (1/2/3/4) ------- : ");
      scanf("%d",&var_type[k]);
      
      if (var_type[k]!=ENUM)  /* If not ENUM, ask for bounds */
	{
	  printf("\nLower and Upper bounds of x[%d] ----- : ",k+1);
	  scanf("%f %f",&x_lower[k],&x_upper[k]);
	}
      
      if (var_type[k]!=CONT)
	var_RIGID[k] = TRUE;
      
      if (var_type[k]==BIN) /* If BIN Type ask for chromosome length */
	{
	  num_bin_var++;
	  BINGA = TRUE;
	  printf("\nChromosome Length ------------------ : ");
	  scanf("%d",&lchrom[k]);
	  TotalChromLength += lchrom[k];
	}
      else
	lchrom[k] = 0;       /* End of BIN */
      
      if (var_type[k]==ENUM)  /* If ENUM Type read the enumerated data */
	{
	  num_enum_var++;
	  REALGA=TRUE;
	  printf("\nEnter the data for variable [%d] in ASCENDING order.\n",k+1);
	  printf("End of data recognized by 999999: ");
	  temp2=0;
	  scanf("%f",&temp1);
	  EnumData[k][0] = temp1;
	  temp2++;
	  while (!small(temp1-999999.0))
	    {
	      scanf("%f",&temp1);
	      EnumData[k][temp2] = temp1;
	      temp2++;
	    }
	  num_enum_data[k]=temp2-1;
	  x_lower[k]=EnumData[k][0];
	  x_upper[k]=EnumData[k][num_enum_data[k]-1];
	}
      else
	num_enum_data[k]=0;  /* End of ENUM */
      
      if (var_type[k]==INT) /* If INT Type form the enumerated data */
	{
	  num_int_var++;
	  REALGA=TRUE;
	  num_enum_data[k] = x_upper[k]-x_lower[k]+1;
	  for (count=0;count<num_enum_data[k];count++)
	    EnumData[k][count]=x_lower[k]+count;
	}                       /* End of INT */
      
      if (var_type[k]==CONT)             /* If CONT Type */
	{
	  num_cont_var++;
	  REALGA=TRUE;
	  printf(" Are these bounds rigid ---- ? (y/n) : ");
	  do { ans = getchar(); } while (ans!= 'y' && ans !='n');
	  if (ans == 'y')  var_RIGID[k] = TRUE;
	  else             var_RIGID[k] = FALSE;
	}
    }  /* End of for loop for variable info. */
  printf("\n");
  
  for (count=0;count<num_var;count++)
    {
      printf("Variable (#%2d) Type : ",count+1);
      switch (var_type[count])
	{
	case BIN  : printf("BINARY\n");break;
	case INT  : printf("INTEGER\n");break;
	case ENUM : printf("ENUMERATED\n");break;
	case CONT : printf("CONTINUOUS\n");break;
	}
      printf("Lower Bound   : %f\n",x_lower[count]);
      printf("Upper Bound   : %f\n",x_upper[count]);
    }
  
  FITNESS=PARAM=FALSE;
  printf("\n Sharing Strategy :"); 
  printf("\n Sharing on Fitness Space   (f)");
  printf("\n Sharing on Parameter Space (p)");
  printf("\n Give your choice (f/p) p preferable : ");
  do { ans = getchar(); } while (ans!= 'f' && ans !='p');
  if (ans == 'f') 	FITNESS = TRUE;
  else                  PARAM   = TRUE;
  
  if (FITNESS)
    {
      printf("\nEqual weightage for all objectives (y or n)?: ");
      do { ans = getchar(); } while (ans!= 'y' && ans != 'n');
      if (ans == 'n') 
	{
	  for (count=0, sumweight=0.0; count<num_obj; count++)
	    {
	      printf("\n  Weight for objective #%2d (0.5) -> ",count+1);
	      scanf("%f", &weightage[count]);
	      sumweight += weightage[count];
	      printf("\n  Lower and Upper values of function #%2d (approx.) -> ",count+1);
	      scanf("%f %f", &afmin[count], &afmax[count]);
	    } 
	  for (count=0; count<num_obj; count++)
	    weightage[count] /= sumweight; 
	}
      else 
	for (count=0; count<num_obj; count++)
	  weightage[count] = 1.0/(double)num_obj;
   }
  /* give a suggestion of dshare */
  printf("\nSigma share value accordingly (%8.3f) -- : ",0.5*pow(0.1,1.0/num_var));
  scanf("%f",&dshare);
  
  printf("\nReports to be printed ----- ? (y/n)  : ");
  do { ans = getchar(); } while (ans!= 'y' && ans !='n');
  if (ans == 'y') REPORT = TRUE;
  else            REPORT = FALSE;
  
  printf("\nPopulation Size ? (~ 20 times N) ---- : ");
  scanf("%d", &pop_size ); 
  if (pop_size > MAXPOPSIZE) 
    {
      printf("\n Increase the value of MAXPOPSIZE in program");
      printf("  and re-run the program");
      exit(-1);
    }
  
  printf("\nCross Over Probability ? ( 0.7 to 1 )  : ");
  scanf("%f",&p_xover);
  
  printf("\nMutation Probability ? ( 0 to 0.2 ) -- : ");
  scanf("%f",&p_mutation);
  
  printf("\n Give the strategy for X-over");
  printf("\n 1 : Single-point crossover");
  printf("\n 2 : Uniform crossover"); 
  printf("\n Give your choice -------------(1/2) : ");
  scanf("%d",&x_strategy);
  
  if (REALGA)
    {
      printf("\nDistribution Index for crossover and mutation (30 50) : ");
      scanf("%f %f",&n_distribution_c, &n_distribution_m);
    }     
  
  printf("\nHow many generations (100) ? --------- : ");
  scanf("%d",&max_gen);
  
  printf("\nGive random seed (0 to 1.0)          : ");
  scanf("%f",&seed);
  
  input_app_parameters();
}

/*====================================================================
  Initialses zero'th generation and global parameters
  ====================================================================*/
initialize()
{
  float u;
  int tmp,k,k1,i,j,j1,stop;
  double temp[MAXVECSIZE],coef;
  unsigned mask=1,nbytes;
  
  randomize();
  app_initialize();
  oldpop = (INDIVIDUAL *)malloc(pop_size*sizeof(INDIVIDUAL));
  newpop = (INDIVIDUAL *)malloc(pop_size*sizeof(INDIVIDUAL));
  if (oldpop == NULL) nomemory("oldpop in initialize()");
  if (newpop == NULL) nomemory("newpop in initialize()");
  
  if (BINGA)    /* If Binary Variables, allocate memory for chromosomes */
    {
      chromsize = (TotalChromLength/UINTSIZE); 
      if (TotalChromLength%UINTSIZE) chromsize++; 
      nbytes = chromsize*sizeof(unsigned); 
      
      for(j = 0; j < pop_size; j++) 
	{ 
	  if((oldpop[j].chrom = (unsigned *) malloc(nbytes)) == NULL) 
            nomemory("oldpop chromosomes"); 
	  
	  if((newpop[j].chrom = (unsigned *) malloc(nbytes)) == NULL) 
            nomemory("newpop chromosomes"); 
	} 
    }             /* End of If (BINGA) */
  
  /* Initializing continuous, integer 
     and enumerated variables    */
  for (k=0; k<= pop_size-1; k++) 
    {  
      oldpop[k].flag = 0;
      oldpop[k].parent1 = oldpop[k].parent2 = 0;
      oldpop[k].dumfitness = 0.0;
      for (j=0; j<=num_var-1; j++)
	{ 
	  if (var_type[j]!=BIN)
	    {
	      u = randomperc();
	      oldpop[k].x[j] = x_lower[j] * (1-u) + x_upper[j] * u;   
	      
	      if ((var_type[j]==INT)||(var_type[j]==ENUM))
		{
		  for(k1=0;k1<num_enum_data[j]-1;k1++)
		    {
		      if (oldpop[k].x[j]<EnumData[j][0]) 
			{
			  oldpop[k].x[j]=EnumData[j][0]; 
			  continue;
			}
		      else if (oldpop[k].x[j]>EnumData[j][num_enum_data[j]-1]) 
			{
			  oldpop[k].x[j]=EnumData[j][num_enum_data[j]-1]; 
			  continue;
			}
		      else if ((EnumData[j][k1]<oldpop[k].x[j])&&
			       (EnumData[j][k1+1]>oldpop[k].x[j]))
                        oldpop[k].x[j]=( (EnumData[j][k1+1]-oldpop[k].x[j])>
					 (oldpop[k].x[j]-EnumData[j][k1]) ) ? 
			  EnumData[j][k1] : EnumData[j][k1+1];
		    }    /* End of for k1 */
		}       /* End of if (INT or ENUM) */
	    }             /* End of if (not BIN)  */
	}
      
      /* Initializing chromosomes
	 for binary variables     */
      for(k1 = 0; k1 < chromsize; k1++) 
	{
	  oldpop[k].chrom[k1] = 0; 
	  if(k1 == (chromsize-1)) 
            stop = TotalChromLength - (k1*UINTSIZE); 
	  else 
            stop = UINTSIZE; 
	  /* A fair coin toss */ 
	  for(j1 = 1; j1 <= stop; j1++) 
	    { 
	      if (flip(0.5)) 
		oldpop[k].chrom[k1] = oldpop[k].chrom[k1]|mask; 
	      if (j1 != stop) 
               oldpop[k].chrom[k1] = oldpop[k].chrom[k1]<<1; 
	    }                 /* End of for (j1) */ 
	}                    /* End of for (k1) */
    }                       /* End of for (k) */
  
  no_xover = no_mutation = 0;
  
  init_dum  = pop_size;             /* For Sharing */
  min_dum   = 0.0;
  delta_dum = 0.1*init_dum;
  
  /* Decoding binary strings and initializing 
     fitness values for each individual */
  for (k=0; k<= pop_size-1; k++)    
    {                                    
      decode_string(&(oldpop[k]));      
      objective(&(oldpop[k]));          
    } 
}

/*====================================================================
  Decodes the string of the individual (if any) and puts the values in
  the array of floats.
  ====================================================================*/
decode_string(ptr_indiv)
     INDIVIDUAL *ptr_indiv;
{
  double *temp,coef[MAXVECSIZE]; 
  int j;
  
  if (ptr_indiv == NULL) error_ptr_null("ptr_indiv in decode_string");
  if (BINGA)
    {
      temp = (double *) malloc(num_var * sizeof(double));
      
      for(j=0; j<num_var; j++) 
	temp[j] = 0.0;
      
      decodevalue(ptr_indiv->chrom,temp);
      
      for(j=0; j<num_var; j++) 
	if (var_type[j]==BIN)
	  {
	    coef[j] = pow(2.0,(double)(lchrom[j])) - 1.0;
	    temp[j] = temp[j]/coef[j];
	    ptr_indiv->x[j] = temp[j]*x_upper[j] + (1.0 - temp[j])*x_lower[j];
	  }
      free(temp);
    }
}

/*====================================================================
  Prints an error message and terminates the program
  ====================================================================*/
nomemory(string) 
     char *string; 
{ 
  printf("\nmalloc: out of memory making %s!!\n",string); 
  printf("\n Program is halting .....");
  exit(-1); 
} 

/*==============================================================
  Gives error message of null pointer  and terminates the program. 
  ==============================================================*/
error_ptr_null(string)
     char *string;
{
  printf("\n Error !! Pointer %s found Null !",string);
  printf("\n Program is halting .....");
  exit(-1);
}

/*====================================================================
  Calculates statistics of current generation :
  ====================================================================*/
statistics(pop,gen)
     POPULATION pop;
     int gen;
{
  int j,iobj;
  float dumsumfit = 0.0;
  float sumfitness[MAXOBJ];
  
  for (iobj=0; iobj<num_obj; iobj++)
    {
      minf[iobj] = maxf[iobj] = 1.0e8;
      sumfitness[iobj] = 0.0;
    }
  
  for (j=0; j<=pop_size-1; j++)
    {
      dumsumfit   += pop[j].dumfitness;
      for (iobj=0; iobj<num_obj; iobj++)
	{
	  sumfitness[iobj] += pop[j].fitness[iobj];
	  if (pop[j].fitness[iobj] > maxf[iobj]) maxf[iobj] = pop[j].fitness[iobj];
	  if (pop[j].fitness[iobj] < minf[iobj]) minf[iobj] = pop[j].fitness[iobj];
	}
    }
  
  dum_avg = dumsumfit/(double)pop_size;
  for (iobj=0; iobj<num_obj; iobj++)
    avgfitness[iobj] = sumfitness[iobj]/(double)pop_size;
  
  app_statistics();
}

/*====================================================================
  Decodes the value of a group of binary strings and puts the decoded
  values into an array 'value'.
  ====================================================================*/
decodevalue(chrom,value) 
     unsigned *chrom; 
     double value[];
{ 
  int temp,k,count,tp,mask=1,sumlengthstring,stop,j,PrevTrack;
  int bitpos,flag1,flag2;
  int CountTrack;
  double pow(), bitpow; 
  
  if (BINGA == FALSE) return;
  if (chrom == NULL) error_ptr_null("chrom in decodevalue");
  
  for(count=0;count<num_var;count++)
    value[count]=0; 
  
  for(k = 0; k < chromsize ; k++) 
    { 
      if (k == (chromsize-1))
	stop = TotalChromLength-(k*UINTSIZE);
      else
	stop = UINTSIZE;
      /* loop thru bits in the current byte */ 
      tp = chrom[k]; 
      for (j=0;j<stop;j++)
	{
	  bitpos=j+k*UINTSIZE;
          /* test for current bit 0 or 1 */ 
          if((tp&mask) == 1) 
	    {
	      sumlengthstring=0;
	      flag1=FALSE;
	      flag2=FALSE;   /* To keep track of the 
			        position of 'j' in the chromosome */
              for(count=0;count<num_var;count++)        
		{                                           
		  if ((var_type[count]==BIN)&&(!flag2))  
		    {                                     
		      if (bitpos>=sumlengthstring)     
			flag1=TRUE;
		      else
			flag1=FALSE;
		      sumlengthstring+=lchrom[count];
		      if ((bitpos<sumlengthstring)&&(flag1)) 
			{  
			  flag2=TRUE;
			  CountTrack=count;
			  PrevTrack=sumlengthstring-lchrom[count];
			}                     /* End of if (bitpos...)  */ 
		    }                        /* End of if (vartype...) */
		}                           /* End of for (count)     */
              
              bitpow = pow(2.0,(double)(bitpos-PrevTrack));
              value[CountTrack] += bitpow;
	    }                   /* End of if (tp&mask) */   
          tp = tp>>1; 
	}                      /* End of for (j) */ 
    }                         /* End of for (k) */ 
} 

/*====================================================================
  GENERATION OF NEW POPULATION through SELECTION, XOVER & MUTATION :
  ====================================================================*/
generate_new_pop()
{
  int j,k,k1,mate1,mate2;
  
  app_computation();
  preselect();
  for (k=0; k<= pop_size-1; k += 2)
    {
      mate1 = select();       /* Stoc. Rem. Roulette Wheel Selection */
      mate2 = select();
      
      switch( x_strategy )          /* Crossover */
	{
	case ONESITE :  cross_over_1_site(mate1,mate2,k,k+1); break;
	case UNIF    :  cross_over_unif(mate1,mate2,k,k+1)  ; break;
	}
      
      mutation(&newpop[k]);        /* Mutation */
      decode_string(&(newpop[k]));
      update_x_BIN_ENUM(&(newpop[k]));
      objective(&(newpop[k]));
      newpop[k].parent1 = mate1+1;
      newpop[k].parent2 = mate2+1;
      
      mutation(&newpop[k+1]);
      decode_string(&(newpop[k+1])); 
      update_x_BIN_ENUM(&(newpop[k+1]));
      objective(&(newpop[k+1]));
      newpop[k+1].parent1 = mate1+1;
      newpop[k+1].parent2 = mate2+1;
    }                          /* For whole population */
}

/***************************************************************
  For integer and enumerated data, a permissible solution is
  calculated here 
  **************************************************************/
update_x_BIN_ENUM(indv)
     INDIVIDUAL *indv;
{
  int j, k1;
  
  for(j=0;j<num_var;j++)
    if ((var_type[j]==INT)||(var_type[j]==ENUM))
      {
	for(k1=0;k1<num_enum_data[j]-1;k1++)
	  {
	    if (indv->x[j]<EnumData[j][0]) 
	      {
		indv->x[j]=EnumData[j][0]; 
		continue;
	      }
	    else if (indv->x[j]>EnumData[j][num_enum_data[j]-1]) 
	      {
		indv->x[j]=EnumData[j][num_enum_data[j]-1]; 
		continue;
	      }
	    else if ((EnumData[j][k1]< indv->x[j])&&
		     (EnumData[j][k1+1]> indv->x[j]))
	      indv->x[j]=( (EnumData[j][k1+1] - indv->x[j]) >
			   (indv->x[j]-EnumData[j][k1]) ) ? 
		EnumData[j][k1] : EnumData[j][k1+1];
	  }                  /* End of for (k1) */
      }
}

/*====================================================================
  Binary cross over routine.
  ====================================================================*/
binary_xover (parent1, parent2, child1, child2) 
     unsigned *parent1, *parent2, *child1, *child2; 
     /* Cross 2 parent strings, place in 2 child strings */ 
{ 
  int j, jcross, k; 
  unsigned mask, temp; 
  
  if (BINGA == FALSE) return;
  if (parent1 == NULL) error_ptr_null("parent1 in binary_xover");
  if (parent2 == NULL) error_ptr_null("parent2 in binary_xover");
  if (child1== NULL) error_ptr_null("child1 in binary_xover");
  if (child2== NULL) error_ptr_null("child2 in binary_xover");
  
  jcross = rnd(1 ,(TotalChromLength - 1));/* Cross between 1 and l-1 */ 
  for(k = 1; k <= chromsize; k++) 
    { 
      if(jcross >= (k*UINTSIZE)) 
	{ 
                child1[k-1] = parent1[k-1]; 
                child2[k-1] = parent2[k-1]; 
	} 
            else if((jcross < (k*UINTSIZE)) && (jcross > ((k-1)*UINTSIZE))) 
            { 
	      mask = 1; 
	      for(j = 1; j <= (jcross-1-((k-1)*UINTSIZE)); j++) 
                { 
		  temp = 1; 
		  mask = mask<<1; 
		  mask = mask|temp; 
                } 
	      child1[k-1] = (parent1[k-1]&mask)|(parent2[k-1]&(~mask)); 
	      child2[k-1] = (parent1[k-1]&(~mask))|(parent2[k-1]&mask); 
            } 
      else 
	{ 
	  child1[k-1] = parent2[k-1]; 
	  child2[k-1] = parent1[k-1]; 
	} 
    } 
} 

/*====================================================================
  Creates two children from parents p1 and p2, stores them in addresses
  pointed by c1 and c2.  low and high are the limits for x values and
  rand_var is the random variable used to create children points.
  ====================================================================*/
create_children(p1,p2,c1,c2,low,high,RIGID)
     float p1,p2,*c1,*c2,low,high;
     int RIGID;
{
  float difference, x_mean, beta, temp;
  float u, u_l, u_u, beta_l=0.0, beta_u=0.0, beta1=0.0, beta2=0.0;
  int flag;
 
 
  if (c1 == NULL) error_ptr_null("c1 in create_children");
  if (c2 == NULL) error_ptr_null("c2 in create_children");
  flag = 0;
  if ( p1 > p2) { temp = p1; p1 = p2; p2 = temp; flag = 1; }
  x_mean = ( p1 + p2) * 0.5;
  difference = p2 - p1;
  if(difference <= 1.0e-9)
    {
      *c1 = p1;
      *c2 = p2;
    }
  else
    {
      if (RIGID)
        {
          beta_l = 1.0 + 2.0*(p1-low)/(p2-p1);
          beta_u = 1.0 + 2.0*(high-p2)/(p2-p1);
          if (MINSCHEME)
            {
              if(beta_l <= beta_u) beta_u = beta_l;
              else beta_l = beta_u;
              u = rndreal(0.0,1.0);
              beta1 = beta2 = get_beta_rigid(u,beta_l);
            }
          else
            {
              u_l = rndreal(0.0,1.0);
              beta1 = get_beta_rigid(u_l,beta_l);
 
              u_u = rndreal(0.0,1.0);
              beta2 = get_beta_rigid(u_u,beta_u);
            }
        }
      else
        {
          u = rndreal(0.0,1.0);
          beta1 = beta2 = get_beta(u);
        }
      *c1 = x_mean - beta1 * 0.5 * difference;
      *c2 = x_mean + beta2 * 0.5 * difference;
    }
  if (flag == 1)
    {
      temp = *c1; *c1 = *c2; *c2 = temp;
    }
  if (*c1 < low)  *c1 = low;
  if (*c2 > high) *c2 = high;
 
  if (*c1 < low || *c2 > high)
    {
      printf("Error!! p1 = %f, p2 = %f, c1 = %f, c2 = %f\n",p1,p2,*c1,*c2);
      exit(-1);
    }
}

/*====================================================================
cross over using strategy of 1 cross site with swapping .
  A random variable is chosen and crossed over. The variables on left
  side are passed as it is and those on right are swapped between the
  two parents.
====================================================================*/
cross_over_1_site(first,second,childno1,childno2)
int first,second,childno1,childno2;
{
    int j,k,site;
   
    if (flip(p_xover))   /* Cross over has to be done */
    {
       no_xover++;
       if (BINGA)
       binary_xover(oldpop[first].chrom   ,oldpop[second].chrom,
                    newpop[childno1].chrom,newpop[childno2].chrom); 

       site = rnd(0,num_var-1);
  
       for (k=0; k<=site-1; k++)   /* Passing the values straight
				      before the cross-site */ 
         if (var_type[site]!=BIN)  /* Only if variable type is not BINARY */
	 {
	    newpop[childno1].x[k] = oldpop[first].x[k]; 
            newpop[childno2].x[k] = oldpop[second].x[k]; 
 	 }

       for (k=site+1; k<=num_var-1; k++)  /* Swapping the values 
				             after the cross-site */ 
         if (var_type[site]!=BIN)
	 {
	    newpop[childno2].x[k] = oldpop[first].x[k]; 
            newpop[childno1].x[k] = oldpop[second].x[k]; 
	 }

                                     /* If variable != BINARY create children 
					at the cross site variable */ 
	if (var_type[site]!=BIN)
	      create_children(oldpop[first].x[site], oldpop[second].x[site], 
			    &(newpop[childno1].x[site]),
			    &(newpop[childno2].x[site]), 
			    x_lower[site],x_upper[site],var_RIGID[site]);
    }                 /* Cross over done */

    else              /* Passing x-values straight */
    {   
      if (BINGA)
      for (k=0; k<=chromsize; k++)
      {
	newpop[childno1].chrom[k] = oldpop[first].chrom[k];
	newpop[childno2].chrom[k] = oldpop[second].chrom[k];
      }
      for (k=0; k<=num_var-1; k++)
      {
	newpop[childno1].x[k] = oldpop[first].x[k];
        newpop[childno2].x[k] = oldpop[second].x[k];
      }
    }
}

/*====================================================================
CROSS - OVER  USING strategy of uniform 50% variables 
  For one variable problem, it is crossed over as usual.
  For multivariables, each variable is crossed over with a probability
  of 50 % , each time generating a new random beta.
====================================================================*/
cross_over_unif(first,second,childno1,childno2)
int first,second,childno1,childno2;
{
    float difference,x_mean,beta;
    int site,k;
   
    if (flip(p_xover))   /* Cross over has to be done */
    {
     no_xover++;
     if (BINGA)
     binary_xover(oldpop[first].chrom,oldpop[second].chrom,
		  newpop[childno1].chrom,newpop[childno2].chrom); 
      for (site = 0; site<=num_var-1; site++)
      {
	if (var_type[site]!=BIN)
          if (flip(0.5) || (num_var==1))
          {
            create_children(oldpop[first].x[site],oldpop[second].x[site],
		      &(newpop[childno1].x[site]),&(newpop[childno2].x[site]),
		      x_lower[site],x_upper[site],var_RIGID[site]);
          }
          else 
          {
	    newpop[childno1].x[site] = oldpop[first].x[site];
            newpop[childno2].x[site] = oldpop[second].x[site];
          }
      }               /* for loop */
    }                 /* Cross over done */

    else              /* Passing x-values straight */
    {   
      if (BINGA)
      for (k=0; k<=chromsize; k++)
      {
	newpop[childno1].chrom[k] = oldpop[first].chrom[k];
	newpop[childno2].chrom[k] = oldpop[second].chrom[k];
      }
      for (site=0; site<=num_var-1; site++)
      {
	newpop[childno1].x[site] = oldpop[first].x[site];
        newpop[childno2].x[site] = oldpop[second].x[site];
      }
    }
}

/*===================================================================
Calculates beta value for given random number u (from 0 to 1)
If input random numbers (u) are uniformly distributed for a set of
inputs, Binary Probability distribution simulation
====================================================================*/
float get_beta(u)
     float u;
{
  float beta;
 
  if (1.0-u < EPSILON ) u = 1.0 - EPSILON;
  if ( u < 0.0) u = 0.0;
  if (u < 0.5) beta = pow(2.0*u,(1.0/(n_distribution_c+1.0)));
  else beta = pow( (0.5/(1.0-u)),(1.0/(n_distribution_c+1.0)));
  return (beta);
}

/*****************************************************************
  Calculates beta for rigid boundaries 
  ****************************************************************/
float get_beta_rigid(u,betaa)
     float u, betaa;
{
  float bet, beta;
 
  if (u <= 0.0+1.0e-9)      beta = 0.0;
  else if (u >= 1.0-1.0e-9) beta = betaa;
  else
    {
      bet = 1.0/betaa;
      bet = 2.0 - pow(bet, (n_distribution_c + 1.0));
      if (u <= 1.0/bet)
        beta = pow(bet * u, (1.0 / (n_distribution_c + 1.0)));
      else
        beta = pow(1.0/(2.0-u*bet), (1.0 / (n_distribution_c + 1.0)));
    }
  return beta;
}

/*==================================================================
For given u value such that   -1 <= u <= 1, this routine returns a
value of delta from -1 to 1. Exact value of delta depends on specified
n_distribution. This is called by mutation().
====================================================================*/
float get_delta(u, delta_l, delta_u)
     float u, delta_l, delta_u;
{
  float delta, aa;
  
  if (u >= 1.0-1.0e-9)      delta = delta_u;
  else if (u <= 0.0+1.0e-9) delta = delta_l;
  else
    {
      if (u <= 0.5)
        {
          aa = 2.0*u + (1.0-2.0*u)*pow((1+delta_l),
				       (n_distribution_m + 1.0));
          delta = pow(aa, (1.0 / (n_distribution_m + 1.0))) - 1.0;
        }
      else
        {
          aa = 2.0*(1-u) + 2.0*(u-0.5)*pow((1-delta_u),
					   (n_distribution_m + 1.0));
          delta = 1.0 - pow(aa, (1.0 / (n_distribution_m + 1.0)));
        }
    }
  if (delta < -1.0 || delta > 1.0)
    {
      printf("Error in mutation!! delta = %f\n",delta);
      exit(-1);
    }
  return (delta);
}

/*==================================================================
Binary mutation routine ( borrowed from sga.c )
====================================================================*/
binmutation(child) 
unsigned *child; 
/* Mutate an allele w/ pmutation, count # of mutations */ 
{ 
    int j, k, stop; 
    unsigned mask, temp = 1; 
 
    if (BINGA == FALSE) return;
    if (child== NULL) error_ptr_null(" child in binmutation");
    for(k = 0; k < chromsize; k++) 
    { 
        mask = 0; 
        if(k == (chromsize-1)) 
            stop = TotalChromLength - ((k-1)*UINTSIZE); 
        else 
            stop = UINTSIZE; 
        for(j = 0; j < stop; j++) 
        { 
            if(flip(p_mutation)) 
            { 
                mask = mask|(temp<<j); 
            } 
        } 
        child[k] = child[k]^mask; 
    } 
} 

/*===================================================================
Mutation Using polynomial probability distribution. Picks up a random 
site and generates a random number u between -1 to 1, ( or between 
minu to maxu in case of rigid boudaries) and calls the routine
get_delta() to calculate the actual shift of the value.
====================================================================*/
mutation(indiv)
INDIVIDUAL  *indiv;
{
   float distance,x,delta,u,delta_l,delta_u;
   int k, site;

   if (indiv == NULL) error_ptr_null("indiv in mutation");

   if (flip(p_mutation))  no_mutation++;
   if(flip (p_mutation) && REALGA)
     {
       site = rnd(0,num_var - 1);
       
       if (var_type[site]!=BIN)
	 {
	   if(fabs(x_upper[site] -x_lower[site]) < EPSILON) return;
	   
	   /* calculation of bounds on delta */  
 	   if(var_RIGID[site])
	     { 
	       x = indiv->x[site];
	       distance = x_lower[site] - x;
	       delta_l = distance / (x_upper[site] - x_lower[site]);
	       if (delta_l < -1.0)  delta_l = -1.0;
	       
	       distance = x_upper[site] - x;
	       delta_u = distance / (x_upper[site] - x_lower[site]);
	       if (delta_u > 1.0)   delta_u = 1.0;
	       
	       if (MINSCHEME) 
		 {    
		   if (-1.0*delta_l < delta_u) delta_u = -1.0 * delta_l;
		   else delta_l = -1.0 * delta_u;
		 }
	     }
	   else  /* flexible */
	     {
	       delta_l = -1.0;
	       delta_u =  1.0;
	     }
	   u = rndreal(0.0, 1.0);
	   delta = get_delta(u, delta_l, delta_u)
	     * (x_upper[site] - x_lower[site]);
	   indiv->x[site] += delta;
	 }
       
     }    /* if flip() */
   
   if (BINGA) binmutation(indiv->chrom);
}

/*====================================================================
  Reporting the user-specified parameters :
fp is the file pointer to output file.
====================================================================*/
initreport(fp)
FILE *fp;
{
   int k, iobj;

   if (fp == NULL) error_ptr_null(" File fp in initreport");
   fprintf(fp,"\n=============================================");
   fprintf(fp,"\n             INITIAL REPORT                  ");
   fprintf(fp,"\n=============================================");
   fprintf(fp,"\n");
   fprintf(fp,"\n Number of objective functions : %2d",num_obj);
   for (iobj=0; iobj<num_obj; iobj++) 
     {
       fprintf(fp,"\n   Objective function #%2d : ",iobj+1); 
       if (minmax[iobj] == 1) fprintf(fp,"Minimize");
       else if (minmax[iobj] == -1) fprintf(fp,"Maximize");
     }
   fprintf(fp,"\n\n CROSSOVER TYPE             : ");
   if (BINGA) fprintf(fp,"Binary GA (Single-pt)");
   fprintf(fp,"\n                              ");
   switch (x_strategy)
   {
    case ONESITE : fprintf(fp,"\n STRATEGY                   : 1 cross - site with swapping"); break;
    case UNIF    : fprintf(fp,"\n STRATEGY                   : Uniformly all variables 50 % "); break;
    default      : fprintf(fp,"\n CROSS OVER NOT SET CORRECTLY "); break;
   }
   fprintf(fp,"\n Population size            : %d",pop_size);
   fprintf(fp,"\n Total no. of generations   : %d",max_gen);
   fprintf(fp,"\n Cross over probability     : %6.4f",p_xover);
   fprintf(fp,"\n Mutation probability       : %6.4f",p_mutation);
   if (BINGA)
      fprintf(fp,"\n String length              : %d",TotalChromLength);
   fprintf(fp,"\n Number of variables");
   fprintf(fp,"\n                Binary      : %d",num_bin_var);
   fprintf(fp,"\n                Integer     : %d",num_int_var);
   fprintf(fp,"\n                Enumerated  : %d",num_enum_var);
   fprintf(fp,"\n                Continuous  : %d",num_cont_var);
   fprintf(fp,"\n                  TOTAL     : %d",num_var);
   fprintf(fp,"\n Epsilon for closeness      : %g",closeness);
   if (REALGA)
     {
       fprintf(fp,"\n Distribution Index (cross) : %g",n_distribution_c);
       fprintf(fp,"\n Distribution Index (mut)   : %g",n_distribution_m);
     }
   fprintf(fp,"\n Sigma-share value          : %6.4f",dshare);
   fprintf(fp,"\n Sharing Strategy           : ");
   if (PARAM) 
      fprintf(fp,"sharing on Parameter Space");
   else if (FITNESS) 
     {
       fprintf(fp,"sharing on Fitness Space");
       for (iobj=0; iobj<num_obj; iobj++)
	 fprintf(fp,"\n Weightage for Obj. Fun. #%2d : %6.4f",iobj+1,weightage[iobj]);
     }
   fprintf(fp,"\n Lower and Upper bounds     :");
   for (k=0; k<=num_var-1; k++)
     fprintf(fp,"\n   %8.4f   <=   x%d   <= %8.4f",x_lower[k],k+1,x_upper[k]);
   fprintf(fp,"\n=================================================\n");
   app_initreport();
}

/*====================================================================
Writes a given string of 0's and 1's
puts a `-` between each substring (one substring for one variable)
Leftmost bit is most significant bit.
====================================================================*/
writechrom(chrom,fp) 
unsigned *chrom; 
FILE *fp;
{ 
    int j, k, stop,bits_per_var,count=0; 
    unsigned mask = 1, tmp; 
    int temp[100];

   for (j=0;j<100;j++)
    temp[j]=0;
   if (fp == NULL) error_ptr_null(" File fp in initreport");
   if (BINGA == FALSE) return;
   if (chrom == NULL) error_ptr_null("chrom in writechrom");
    for(k = 0; k < chromsize; k++) 
    { 
        tmp = chrom[k]; 
        if(k == (chromsize-1)) 
            stop = TotalChromLength - (k*UINTSIZE); 
        else 
            stop = UINTSIZE; 
        for(j = 0; j < stop; j++) 
        { 
            if(tmp&mask) 
	    {
                /*fprintf(fp,"1"); */
		temp[j+count*UINTSIZE]=1;
	    }
            else 
	    {
                /*fprintf(fp,"0"); */
		temp[j+count*UINTSIZE]=0;
	    }
            tmp = tmp>>1; 
        } 
	count++;
    } 

    count=0;
    for(j=0;j<num_var;j++)
    {
       for(k=count+lchrom[j]-1;k>=count;k--)
       {
	  fprintf(fp,"%d",temp[k]);
       }
       count+=lchrom[j];
       if (j!=num_var-1) fprintf(fp,"_");
    }
} 

/*====================================================================
Reporting the statistics of current population ( gen. no. 'num'):
  fp is file pointer to output file.
====================================================================*/
report(fp,num)
FILE *fp;
int num;
{
  int k,j,iobj; 
  char string[30];

  if (fp == NULL) error_ptr_null(" file fp in report()");
   /* ----------------------------------------- */
   /* WRITING IN THE OUTPUT FILE FOR INSPECTION */
   /* ----------------------------------------- */
   fprintf(fp,"\n======================== Generation # : %3d =================================",num);
   if (BINGA) fprintf(fp,"\n  No.   Rank         x   Obj. Fun. Values (f1,f2, etc.) Parents   String");
   else fprintf(fp,"\n  No.   Rank         x   Obj. Fun. Values (f1,f2, etc.   Parents  ");
   fprintf(fp,"\n-----------------------------------------------------------------------------");
   for (k=0; k<= pop_size-1; k++)
   {
     fprintf(fp,"\n %3d.   %3d   [%7.3f ] ",k+1,oldpop[k].front,oldpop[k].x[0]);
     for (j= 1; j<=num_var-1; j++)
       fprintf(fp,"\n              [%7.3f ] ",oldpop[k].x[j]);
     for (j=0;j<num_obj;j++) fprintf(fp," %9.3f ", oldpop[k].fitness[j]);
     fprintf(fp," (%3d %3d)", oldpop[k].parent1, oldpop[k].parent2);
     if (BINGA) 
     {  
	fprintf(fp,"     ");
        writechrom(oldpop[k].chrom,fp); 
     }
   }
   fprintf(fp,"\n-----------------------------------------------------------------------------");
   fprintf(fp,"\nAverage Dummy Fitness = %8.3f", dum_avg);
   for(iobj=0; iobj < num_obj; iobj++)
     fprintf(fp,"\nObj. Function #%2d: Max. Fitness = %8.3f  Min. Fitness = %8.3f  Avg. Fitness = %8.3f",iobj+1,maxf[iobj],minf[iobj],avgfitness[iobj]);
   fprintf(fp,"\nNo. of mutations = %d ;  No. of x-overs = %d",
				no_mutation,no_xover);

   fprintf(fp,"\n=============================================================================");
   fprintf(fp,"\n\n");
   
  app_report();
}

/*====================================================================
Reporting the statistics of current population ( gen. no. 'num'):
  fp is file pointer to output file.
====================================================================*/
result(fp,num)
FILE *fp;
int num;
{
  int k,j; 
  char string[30];


  if (fp == NULL) error_ptr_null(" file fp in report()");

   /* --------------------------------------*/
   /* WRITING IN THE OUTPUT FILE FOR RESULTS*/
   /* --------------------------------------*/
   fprintf(fp,"\n#============== Generation # : %3d ===========================================",num);
   if (BINGA) 
      fprintf(fp,"\n#  No.          x         Obj. Fun. Values (f1,f2,etc.)     String");
   else 
      fprintf(fp,"\n#  No.          x         Obj. Fun. Values (f1,f2,etc.)  ");
   fprintf(fp,"\n#=============================================================================");
   for (k=0; k<= pop_size-1; k++)
   {
     /* for now deb 14/8/98
    fprintf(fp,"\n %3d. x%d = [%10.3f ]  ",k+1,1,oldpop[k].x[0]);
     for (j= 1; j<=num_var-1; j++)
       fprintf(fp,"\n      x%d = [%10.3f ]   ",j+1,oldpop[k].x[j]);
       */
     for (j=0; j<num_obj; j++)
       fprintf(fp," %10.3f ", oldpop[k].fitness[j]);
     /*  if (BINGA) 
       {
	 fprintf(fp,"        ");
	 writechrom(oldpop[k].chrom,fp);
       }
     */
     fprintf(fp,"\n");
   }
   fprintf(fp,"\n#=============================================================================");
   fprintf(fp,"\n\n");
   
   app_report();
}

/*====================================================================
Releases the memory for all mallocs
====================================================================*/
free_all()
{
  int i; 
   
  for(i = 0; i < pop_size; i++) 
  {   
      free(oldpop[i].chrom); 
      free(newpop[i].chrom); 
  } 
  free(oldpop);
  free(newpop);
  select_free();
  app_free();
}

/*====================================================================
Population Fronts Created Here 

     Classifying entire population into different fronts 
     and assigning dummy fitness(dumfitness) accordingly.     

     flags as follows  :  flag = 0  for untouched
                                 1      dominated
                                 2      non-dominated
                                 3      exhausted             
====================================================================*/
MakeFronts(void)
{
  int i,j,front_index,pop_count,ijk, iobj, flagobj;
  
  for(i=0; i<pop_size; i++)                /* initializing */
    {
      oldpop[i].flag = 0;   /* making all indivs. untouched */
      oldpop[i].dumfitness = 0.0;
    }
  
  pop_count = 0;   /* for checking if all individuals are 
		      assigned a front */
  
  front_index = 1; /* first front */
  while(pop_count < pop_size)
    {
      for(j=0; j<pop_size; j++)
	{
	  if(oldpop[j].flag == 3) continue;/* already assigned
					      a front, do not consider*/
	  
	  for(i=0; i<pop_size; i++)
	    {
	      if(i == j) continue; /* one is not compared with the same */
	      
	      else if(oldpop[i].flag == 3) continue; /* already assigned*/
	      
	      else if(oldpop[i].flag == 1) continue; /* marked dominated*/
	      
	      else   /* check for domination */
		{
		  flagobj = 0;
		  for (iobj=0; iobj<num_obj && !flagobj; iobj++) 
		    if(minmax[iobj] * oldpop[j].fitness[iobj] <= minmax[iobj] * oldpop[i].fitness[iobj])
		      flagobj = 1;
		  if (flagobj == 0)
		    { oldpop[j].flag = 1; break; }
		}		  
	    }                       /*End of For loop --i--*/
	  
	  if(oldpop[j].flag == 0)  /* non-dominated solutions */
	    {
	      oldpop[j].flag = 2 ;
	      pop_count++ ;
	    }
	}                           /*End of For loop --j--*/
      
      if(front_index == 1)
	{                      /* all in first front are assigned a
                                  dummy fitness = popsize */
	  for(i=0; i<pop_size; i++)
	    if(oldpop[i].flag == 2) 
	      {
		oldpop[i].dumfitness = init_dum;
	        oldpop[i].front = front_index;
	      }
	  
	  phenoshare();  /*Phenotypic sharing of non-dominated strings*/
	  
	  minimum_dum();        /* Finding minimum dummy-fitness 
				   among shared strings*/
	  
	  front_index++ ;   /* incremented for next front */
	} 
      
      else {      /* for all other fronts 2, 3, ... */
	for(i=0; i<pop_size; i++)
	  if(oldpop[i].flag == 2) /* member of current front */
	    {  
	      oldpop[i].front = front_index ;
	      if(min_dum > delta_dum) 
		oldpop[i].dumfitness = min_dum-delta_dum;
	      /* smaller than the smallest dummy fitness in 
		 previous front */
	      else adjust(front_index);
	    }
	
	phenoshare();
	
	minimum_dum(); 
	
	front_index++ ;
	
      }                                /* --if else-- */
      
      for(i=0; i<pop_size; i++)
	{
	  /* call members of current front exhausted */ 
	  if(oldpop[i].flag == 2) oldpop[i].flag =  3;
	  
	  /* and ummark all current dominated solutions */
	  else if(oldpop[i].flag == 1) oldpop[i].flag = 0;
	}
      
    }                                    /*End of while loop */
  
}

adjust(index)
     int index;
{      /* jack up the fitness of all solutions assigned in a front
          to accomodate remaining ones */
  int i;
  double diff;
  
  diff = 2.0 * delta_dum - min_dum ;
  for(i=0; i<pop_size; i++)
    if(oldpop[i].flag == 1 || oldpop[i].flag == 0) continue;
    else oldpop[i].dumfitness += diff ;
  
  minimum_dum();
}

/*====================================================================
  Sharing in a front.  oldpop.dumfitness is divided by nichecount.
  =====================================================================*/
phenoshare()
{
  
  int i,j;
  float  dvalue,d,nichecount;
  double pow();
  float  distance();
  
  for(j=0; j<pop_size; j++)
    {
      nichecount = 1.0;
      if(oldpop[j].flag == 2)
	{
	  for(i=0; i<pop_size; i++)
	    {
	      if(i == j) continue;
	      
	      if(oldpop[i].flag == 2)
		{
		  /** distance() returns the phenotypic distance between
		    two individuals **/
		  d = distance(&(oldpop[j]),&(oldpop[i]));
		  if (d < 0.0) d = (-1.0)*d ;
		  if (d <= 0.000001) nichecount++ ;
		  else if(d < dshare) 
		    nichecount += (1.0-(d/dshare))*(1.0-(d/dshare));
	        }
	    }   /* for i loop */
	}    /* for oldpop[j].flag==2 loop */
      oldpop[j].dumfitness /= nichecount ;
      
    } /* j loop */
}

/*======================================================
  distance() returns the phenotypic distance
  between two individuals :  
  o in n-dimensional space. (number of variables = nx)
  o in fitness space.       (fitness1 - fitness2 space)   
  =======================================================*/
float  distance(critter1,critter2)
     INDIVIDUAL *critter1;
     INDIVIDUAL *critter2;
{
  int i, iobj;
  double dst,sqrt();
  
  dst = 0.0;
  if (FITNESS)             /* Sharing on fitness space */
    {
      for(iobj=0; iobj<num_obj; iobj++)
        dst += weightage[iobj]*square(critter1->fitness[iobj] - critter2->fitness[iobj])/square(afmax[iobj]-afmin[iobj]);
    }
  else if (PARAM)         /* Sharing on parameter space */
    {
      for (i=0;i<num_var;i++)
	dst += square(critter1->x[i]-critter2->x[i])/square(x_upper[i]-x_lower[i]);
    }
  dst = sqrt(dst);
  return(dst);
}

minimum_dum()
{    /* finding the minimum dummy fitness in the current front */
  int i;
  
  min_dum = 1000000000.0 ;
  
  for(i=0; i<pop_size; i++)
    {
      if(oldpop[i].flag == 2)
	{
	  if(oldpop[i].dumfitness < min_dum) min_dum = oldpop[i].dumfitness;
	}
    }
}

/*====================================================================
  MAIN PROGRAM ;
  ====================================================================*/
main()
{
  FILE *fp_out;      	/* File pointer for output file  	*/ 
  FILE *fp_report;     /* File pointer for report */ 
  int j,k,k1;
  POPULATION 	temp;	/* A temporary pointer of population 	*/
  
  /*---------------------------*/
  /* Program starts here :     */ 
  /*---------------------------*/
  printf("**********************************************************************\n");
  printf("***             MULTI-OBJECTIVE OPTIMIZATION SOFTWARE              ***\n");
  printf("***                                                                ***\n");
  printf("***        Developed by : Kalyanmoy Deb and Mayank Goyal           ***\n");
  printf("***           The Department of Mechanical Engineering             ***\n");
  printf("***        Indian Institute of Technology, Kanpur - INDIA          ***\n");
  printf("***................................................................***\n");
  printf("***  This software gives a set of solutions to a multi-objective   ***\n");
  printf("***  optimization problem.  The parameters can be represented as   ***\n");
  printf("***  any of the following following type (or their combination) :  ***\n");
  printf("***              o BINARY STRING TYPE                              ***\n");
  printf("***              o INTEGER TYPE                                    ***\n");
  printf("***              o ENUMERATED DATA TYPE                            ***\n");
  printf("***              o REAL/CONTINUOUS TYPE                            ***\n");
  printf("***                                                                ***\n");
  printf("***  For Multiobjective optimization, two strateies are there :    ***\n");
  printf("***              o Share on FITNESS SPACE                          ***\n");
  printf("***              o Share on PARAMETER SPACE                        ***\n");
  printf("***  Sigma-share value should be chosen accordingly.               ***\n");
  printf("***                                                                ***\n");
  printf("***  You can also assign importance of one fitness function over   ***\n");
  printf("***  other by giving them weightages.                              ***\n");
  printf("\n***  All rights reserved. Not to be used for commercial purposes.  ***\n");
  printf("\n***  Please send bug information and your comments to              *** ");
  printf("\n***    deb@ls11.informatik.uni-dortmund.de or deb@iitk.ernet.in    ***");
  printf("\n**********************************************************************");
  
  input_parameters();
  
  fp_out = fopen("result.out","w+");
  fp_report = fopen("report","w+");
  
  select_memory(); 
  //  initreport(fp_out);   
  initreport(fp_report);   
  
  gen_no=0;
  initialize();
  MakeFronts();
  statistics(oldpop,gen_no);
  
  // result(fp_out,gen_no);  
  if (REPORT) report(fp_report,gen_no);  
  
  printf("\n =====================================");
  printf("======================================== ");
  printf("\n Please Wait ");
  
  for(gen_no = 1; gen_no<=max_gen; gen_no++)
    {
      printf("."); if (gen_no%60==0) printf("\n");
      fflush(stdout);
      
      generate_new_pop();
      
      temp   = oldpop; 
      oldpop = newpop; 
      newpop = temp  ; 
      
      MakeFronts();
      statistics(oldpop,gen_no);
      
      if (gen_no==max_gen) result(fp_out,gen_no);  
      if (REPORT) report(fp_report,gen_no);  
    };                        /* One GA run is over  */
  printf("\n ============================================================================= ");
  
  free_all();
  
  fclose(fp_out);
  fclose(fp_report);
  app_closure();
  printf("\n Results are stored in file 'result.out' ");
  puts("\n O.K Good bye !!!"); 
}

/**************** End of Main Program ***************************/

/*-------------------------------------------------------  */ 
/* random.c - contains random number generator and related */
/* utilities,                                              */
/* Source : sga.c  (c) E.Goldberg 1986
   /*-------------------------------------------------------  */ 

/* variables are declared static so that they cannot       */
/* conflict with names of other global variables in other  */
/* files.  See K&R, p 80, for scope of static              */ 

static double oldrand[55];   /* Array of 55 random numbers */ 
static int jrand;                 /* current random number */ 
static double rndx2;    /* used with random normal deviate */ 
static int rndcalcflag; /* used with random normal deviate */ 

int small(float number)
{
  if (fabs(number)<=EPSILON) 
    return TRUE;
  else 
    return FALSE;
}

initrandomnormaldeviate() 
     /* initialization routine for randomnormaldeviate */ 
{ 
  rndcalcflag = 1; 
} 

double noise(mu ,sigma) 
     /* normal noise with specified mean & std dev: mu & sigma */ 
     double mu, sigma; 
{ 
  double randomnormaldeviate(); 
  
  return((randomnormaldeviate()*sigma) + mu); 
} 

double randomnormaldeviate() 
     /* random normal deviate after ACM algorithm 267 / Box-Muller Method */ 
{ 
  double sqrt(), log(), sin(), cos(); 
  float randomperc();  
  double t, rndx1; 
  
  if(rndcalcflag) 
    { 
      rndx1 = sqrt(- 2.0*log((double) randomperc())); 
      t = 6.2831853072 * (double) randomperc(); 
      rndx2 = sin(t); 
      rndcalcflag = 0; 
      return(rndx1 * cos(t)); 
    } 
  else 
    { 
      rndcalcflag = 1; 
      return(rndx2); 
    } 
} 

advance_random() 
     /* Create next batch of 55 random numbers */ 
{ 
  int j1; 
  double new_random; 
  
  for(j1 = 0; j1 < 24; j1++) 
    { 
      new_random = oldrand[j1] - oldrand[j1+31]; 
      if(new_random < 0.0) new_random = new_random + 1.0; 
      oldrand[j1] = new_random; 
    } 
  for(j1 = 24; j1 < 55; j1++) 
    { 
      new_random = oldrand [j1] - oldrand [j1-24]; 
      if(new_random < 0.0) new_random = new_random + 1.0; 
      oldrand[j1] = new_random; 
    } 
} 

int flip(prob) 
     /* Flip a biased coin - true if heads */ 
     float prob; 
{ 
  float randomperc(); 
  
  if(randomperc() <= prob) 
    return(1); 
  else 
    return(0); 
} 

randomize() 
     /* Get seed number for random and start it up */ 
{ 
  int j1; 
  
  for(j1=0; j1<=54; j1++) oldrand[j1] = 0.0; 
  jrand=0; 
  
  warmup_random(seed); 
  initrandomnormaldeviate(); 
} 

float randomperc() 
     /* Fetch a single random number between 0.0 and 1.0 -  */
     /* Subtractive Method . See Knuth, D. (1969), v. 2 for */ 
     /* details.Name changed from random() to avoid library */
     /* conflicts on some machines                          */ 
{ 
  jrand++; 
  if(jrand >= 55) 
    { 
      jrand = 1; 
      advance_random(); 
    } 
  return((float) oldrand[jrand]); 
} 

int rnd(low, high) 
     /* Pick a random integer between low and high */ 
     int low,high; 
{ 
  int i; 
  float randomperc(); 
  
  if(low >= high) 
    i = low; 
  else 
    { 
      i = (randomperc() * (high - low + 1)) + low; 
      if(i > high) i = high; 
    } 
  return(i); 
} 

float rndreal(lo ,hi) 
     /* real random number between specified limits */ 
     float lo, hi; 
{ 
  return((randomperc() * (hi - lo)) + lo); 
} 

warmup_random(random_seed) 
     /* Get random off and running */ 
     float random_seed; 
{ 
  int j1, ii; 
  double new_random, prev_random; 
  
  oldrand[54] = random_seed; 
  new_random = 0.000000001; 
  prev_random = random_seed; 
  for(j1 = 1 ; j1 <= 54; j1++) 
    { 
      ii = (21*j1)%54; 
      oldrand[ii] = new_random; 
      new_random = prev_random-new_random; 
      if(new_random<0.0) new_random = new_random + 1.0; 
      prev_random = oldrand[ii]; 
    } 
  
  advance_random(); 
  advance_random(); 
  advance_random(); 
  
  jrand = 0; 
} 

/*----------------------------------------------------------*/
/* Files for tournament selection :                         */
/* Source : sga.c (c) E.Goldberg                            */
/*----------------------------------------------------------*/

select_memory() 
{  /* allocates auxiliary memory for stochastic remainder selection*/
  
  unsigned nbytes;
  int j;
  
  choices = NULL;
  fraction = NULL;
  nbytes = pop_size*sizeof(int);
  if((choices = (int *) malloc(nbytes)) == NULL)
    nomemory(stderr,"choices");
  nbytes = pop_size*sizeof(float);
  if((fraction = (float *) malloc(nbytes)) == NULL)
    nomemory(stderr,"fraction");
}

select_free()
{
  /* frees auxiliary memory for stochastic remainder selection */
  choices = NULL;
  fraction = NULL;
  free(choices);
  free(fraction);
}

preselect()
     /* preselection for stochastic remainder method */
{
  int j, jassign, k;
  float expected;
  int flip();
  
  if(dum_avg == 0)
    {
      for(j = 0; j < pop_size; j++) choices[j] = j;
    }
  else
    {
      j = 0;
      k = 0;
      /* Assign whole numbers */
      do
        {
	  expected = (float)((oldpop[j].dumfitness)/dum_avg);
	  jassign = (int)expected;
	  /* note that expected is automatically truncated */
	  fraction[j] = expected - (float)jassign;
	  while(jassign > 0)
	    {
	      jassign--;
	      choices[k] = j;
	      k++;
	    }
	  j++;
        }
      while(j < pop_size);
      
      j = 0;
      /* Assign fractional parts */
      while(k < pop_size)
	{
	  if(j >= pop_size) j = 0;
	  if(fraction[j] > 0.0)
	    {
	      /* A winner if true */
	      if(flip(fraction[j]))
		{
		  choices[k] = j;
		  fraction[j] = fraction[j] - 1.0;
		  k++;
		}
	    }
	  j++;
	}
    }
  nremain = pop_size - 1;
}

int select()
     /* selection using remainder method */
{
  int jpick, slect;
  int rnd();
  
  jpick = rnd(0, nremain);
  slect = choices[jpick];
  choices[jpick] = choices[nremain];
  nremain--;
  return(slect);
}


reset1()    
     /* Name changed from reset because of clash with lib. function - RBA */
     /* Shuffles the tourneylist at random */ 
{ 
  int i, rand1, rand2, temp_site; 
  
  for(i=0; i<pop_size; i++) tourneylist[i] = i; 
  
  for(i=0; i < pop_size; i++) 
    { 
      rand1= rnd(0,pop_size-1);
      rand2=  rnd(0,pop_size-1);
      temp_site = tourneylist[rand1]; 
      tourneylist[rand1]=tourneylist[rand2]; 
      tourneylist[rand2]=temp_site; 
    } 
}   
/******************* APPLICATION ORIENTED ROUTINES ***************/
/**** Change these routines for your particular application ******/

input_app_parameters()
     /* Input your application dependent parameters here and put the
	output in global variables */
{
}

app_computation() 
     /* this routine should contain any application-dependent computations */ 
     /* that should be performed before each GA cycle.
	called by generate_new_pop    */ 
{ 
} 

app_free() 
     /* application dependent free() calls, called by free_all() */ 
{ 
} 

app_initialize() 
     /* application dependent initialization routine called by intialize() */ 
{ 
} 


app_initreport() 
     /* Application-dependent initial report called by initreport() */ 
{ 
} 

app_report() 
     /* Application-dependent report, called by report() */ 
{ 
} 


app_statistics() 
     /* Application-dependent statistics calculations called by statistics() */ 
{ 
} 

app_closure()
     /* Application-dependent work which should be done before closure of
	the main program. called by main() */
{
}

/*====================================================================
  OBJECTIVE FUNCTION : Change it for different applications
  
  Given a sample application where the two objective functions are :
  f1(x) = x*x;
  f2(x) = (x-2)*(x-2)  and the aim is to minimize both of them.
  
  NORMALIZE your fitness functions before use. (Divide them by their
  maximum possible values)
  
  The code is designed for minimizing the fitness functions.
  To maximize a function, either use
  o 1/(1+f(x)) or
  o -f(x)  instead of f(x).
  ===================================================================*/
objective(person)
     INDIVIDUAL *person;
{
  float a, a1, a2;
  float penalty, g[10];
  int i, nc;
   
  if (person == NULL) error_ptr_null("person in objective()");
  
  /* First problem */
#ifdef f1
  a=person->x[0];
  person->fitness[0] = square(a); 
  person->fitness[1] = square((2.0-a));
  nc = 0;
#endif
  
  /* Second problem */
#ifdef f2
  a=person->x[0];
  person->fitness[0]=(a<=1.0) ? -a : ((a<=3) ? -2+a : ((a<=4) ? 4-a : -4+a));
  person->fitness[1]=(a-5)*(a-5); 
  nc = 0;
#endif
  
  /* Third problem */
#ifdef f3
  a1 = person->x[0];
  a2 = person->x[1];
  person->fitness[0]=(a1-2)*(a1-2)+(a2-1)*(a2-1)+2;
  person->fitness[1]=9*a1-(a2-1)*(a2-1);
  
  nc = 2;
  g[0] = -1.0*(a1*a1 + a2*a2 - 225.0);
  g[1] = -1.0*(a1 - 3.0*a2 + 10.0);
#endif
	
#ifdef book
  a1 = person->x[0];
  a2 = person->x[1];
  person->fitness[0] = a1;
  person->fitness[1] = (1.0+a2)/a1;
  nc = 0;
#endif

  penalty = 0.0;
  for (i=0; i<nc; i++)
    if (g[i] < 0.0) penalty += PENALTY_COEFF * g[i] * g[i];
  
  for (i=0; i<num_obj; i++)
    person->fitness[i] = person->fitness[i] + minmax[i] * penalty;
}
/********************** END  OF  FILE **************************/




