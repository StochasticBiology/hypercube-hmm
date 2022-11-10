// this is the workhorse code for the HyperTraPS algorithm
// it takes command line arguments that dictate the data file(s), the structure of the inference run, and various parameters
// the output file is a set of samples from the posterior distribution inferred over hypercube parameters

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define RND drand48()

// maximum number of datapoints (just for memory allocation)
#define _MAXN 2000

// number of trajectories N_h, and frequencies of sampling for posteriors and for output
#define BANK 2000
#define TMODULE 100

#define _EVERYITERATION 0

// control output
#define VERBOSE 0
int SPECTRUM_VERBOSE = 0;

// produce gaussian random number
double gsl_ran_gaussian(const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * RND;
      y = -1 + 2 * RND;

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

// pick a new locus to change in state "state"; return it in "locus" and keep track of the on-course probability in "prob". "ntrans" is the transition matrix
void PickLocus(int *state, double *ntrans, int *targ, int *locus, double *prob, double *beta, int LEN)
{
  int i, j;
  double *rate;
  double totrate, nobiastotrate;
  double *cumsum;
  double r;

  rate = (double*)malloc(sizeof(double)*LEN);
  cumsum = (double*)malloc(sizeof(double)*LEN);
  
  nobiastotrate = 0;

  /* compute the rate of loss of gene i given the current genome */
  for(i = 0; i < LEN; i++)
    {
      /* ntrans must be the transition matrix. ntrans[0]-ntrans[LEN] are the bare rates. then ntrans[LEN+j*LEN+i] is the modifier for i from j*/
      if(state[i] == 0)
	{
	  rate[i] = ntrans[i];
	  for(j = 0; j < LEN; j++)
	    rate[i] += state[j]*ntrans[LEN+j*LEN+i];
	  rate[i] = exp(rate[i]);
	}
      else /* we've already lost this gene */
	rate[i] = 0;

      /* roulette wheel calculations as normal */
      cumsum[i] = (i == 0 ? 0 : rate[i-1]+cumsum[i-1]);
      nobiastotrate += rate[i];
    }

  totrate = 0;

  /* compute the rate of loss of gene i given the current genome */
  for(i = 0; i < LEN; i++)
    {
      /* ntrans must be the transition matrix. ntrans[0]-ntrans[LEN] are the bare rates. then ntrans[LEN+j*LEN+i] is the modifier for i from j*/
      if(state[i] == 0 && targ[i] != 0)
	{
	  rate[i] = ntrans[i];
	  for(j = 0; j < LEN; j++)
	    rate[i] += state[j]*ntrans[LEN+j*LEN+i];
	  rate[i] = exp(rate[i]);
	}
      else /* we've already lost this gene OR WE DON'T WANT IT*/
	rate[i] = 0;

      /* roulette wheel calculations as normal */
      cumsum[i] = (i == 0 ? 0 : rate[i-1]+cumsum[i-1]);
      totrate += rate[i];
    }

  /* normalised, additive rates -- is this sensible? */
  for(i = 0; i < LEN; i++)
    cumsum[i] /= totrate;

  r = RND;
  for(i = 0; i < LEN-1; i++)
    {
      if(cumsum[i] < r && cumsum[i+1] > r) { break; }
    }

  *locus = i;

  *prob = totrate/nobiastotrate;
  *beta = nobiastotrate;

  free(rate);
  free(cumsum);
}


// compute HyperTraPS probability of a transition from "startpos" to "targ" given transition matrix "P"
double LikelihoodMultiple(int *targ, double *P, int LEN, int *startpos, double tau1, double tau2)
{
  int *bank;
  int n0, n1;
  double *reject;
  int i, j, r;
  int locus;
  int *attempt;
  double min;
  double mean;
  double *prodreject;
  double *summand;
  int fail, score;
  int *hits;
  double totalsum;

  // new variables
  double u, prob_path, vi, betaci, nobiastotrate;
  double analyticI1, analyticI2;
  double sumI1, sumI2;
  int n;
  double tmprate;
  double *recbeta;
  // nobiastotrate is retain to match role in PickLocus but basically corresponds to -u
  int exitcount = 0;
  
  // allocate memory for BANK (N_h) trajectories
  bank = (int*)malloc(sizeof(int)*LEN*BANK);
  reject = (double*)malloc(sizeof(double)*BANK);
  hits = (int*)malloc(sizeof(int)*BANK);
  prodreject = (double*)malloc(sizeof(double)*BANK);
  recbeta = (double*)malloc(sizeof(double)*LEN*BANK);

  attempt = (int*)malloc(sizeof(int)*LEN);
  summand = (double*)malloc(sizeof(double)*LEN);
  
  // initialise each trajectory at the start state; count 0s and 1s
  for(i = 0; i < LEN*BANK; i++)
    bank[i] = startpos[i%LEN]; 
  n0 = 0;
  for(i = 0; i < LEN; i++)
    n0 += startpos[i];

  n1 = 0;
  for(i = 0; i < LEN; i++)
    n1 += (targ[i] == 1 || targ[i] == 2);

  if(n0 > n1)
    {
      // the target comes before the source
      printf("Wrong ordering, or some other problem with input file. Remember pairs of rows should be ancestor then descendant.\n");
      exit(0);
    }

  mean = 1;
  totalsum = 0;

  // check we're not already there
  fail = 0;
  for(i = 0; i < LEN; i++)
    fail += (targ[i] != startpos[i]);
  if(fail == 0) totalsum = 1;

  for(i = 0; i < BANK; i++)
    prodreject[i] = 1;

  // loop through the number of evolutionary steps we need to make
  for(r = n0; r < n1; r++)
    {
      for(i = 0; i < BANK; i++)
	hits[i] = 0;

      // loop through each trajectory
      for(i = 0; i < BANK; i++)
	{
	  for(j = 0; j < LEN; j++)
	    attempt[j] = bank[LEN*i+j];
	  // pick the locus to change at this step, and record the probability that we stay on track to the target
	  PickLocus(&bank[LEN*i], P, targ, &locus, &reject[i], &recbeta[LEN*i + (r-n0)], LEN);
	  bank[LEN*i+locus] = 1;

	  fail = 0;
	  // count whether we're there or not
	  for(j = 0; j < LEN; j++)
	    {
	      if(bank[LEN*i+j] != targ[j] && targ[j] != 2) fail = 1;
	    }
	  hits[i] += (1-fail);

	}

      // keep track of total probability so far, and record if we're there
      summand[r] = 0;
      for(i = 0; i < BANK; i++)
	{
	  prodreject[i] *= reject[i];
	  summand[r] += prodreject[i]*hits[i];
	}
      summand[r] /= BANK;

      totalsum += summand[r];
    }

  if(n0 == n1)
    {
      prob_path = 1;
    }
  else
    {
      prob_path = 0;
      for(n = 0; n < BANK; n++)
	{
	  prob_path += prodreject[n]*1./BANK;
	}
    }
          
  free(bank);
  free(reject);
  free(hits);
  free(prodreject);
  free(recbeta);

  free(attempt);
  free(summand);

  return prob_path;
}

// get total likelihood for a set of changes
double GetLikelihoodCoalescentChange(int *matrix, int len, int ntarg, double *ntrans, int *parents, double *tau1s, double *tau2s)
{
  double loglik, tloglik, tlik;
  int i, j;
  int multiple;
  int *startpos;

  startpos = (int*)malloc(sizeof(int)*len);
  
  // initialise and start at one corner of the hypercube
  loglik = 0;
  for(i = 0; i < len; i++)
    startpos[i] = 0;

  // loop through each ancestor/descendant pair
  for(i = 0; i < ntarg/2; i++)
    {
      // output if desired
      if(VERBOSE)
	{
	  printf("Target %i: ", i);
	  for(j = 0; j < len; j++) printf("%i", matrix[2*i*len+len+j]);
	  printf(" parent is: " );
	  for(j = 0; j < len; j++) printf("%i", matrix[2*i*len+j]); 
	  printf("\n");
	}
      // get log-likelihood contribution from this pair (transition) using HyperTraPS
       tlik = LikelihoodMultiple(&(matrix[2*i*len+len]), ntrans, len, &(matrix[2*i*len]), tau1s[i], tau2s[i]);
      tloglik = log(tlik);
      if(tlik < 0)
	{
	  printf("Somehow I have a negative likelihood, suggesting a lack of numerical convergence. Terminating to avoid unreliable posteriors.\n");
	  exit(0);
	}

      // output if required
      if(VERBOSE)
	printf("--- %i %f %f\n", i, exp(tloglik), tloglik);
      loglik += tloglik;
    }

  free(startpos);
  
  // return total log likelihood
  return loglik;
}
 
// main function processes command-line arguments and run the inference loop
int main(int argc, char *argv[])
{
  int parents[_MAXN];
  FILE *fp;
  int *matrix;
  int len, ntarg;
  double *trans, *ntrans;
  int t;
  int i, j;
  char ch;
  double lik, nlik;
  int *rec, *tmprec;
  int maxt, allruns;
  int seed;
  char str[200];
  char fstr[200];
  char shotstr[200], bestshotstr[200];
  double DELTA, MU;
  int NVAL;
  int expt;
  double acc, rej, lacc, lrej;
  int chain1, chain2;
  double prob;
  double *tmpmat;
  double r;
  char fstr1[100], fstr2[100];
  time_t timer;
  char buffer[25];
  struct tm* tm_info;
  double taus[_MAXN], tau1s[_MAXN], tau2s[_MAXN];
  int ntau;
  int nancount = 0;
  int spectrumtype;
  double bestlik = 0;
  int lengthindex, kernelindex;
  int SAMPLE;
  int losses;

  printf("\nHyperTraPS\n\n");
  
  // process command-line arguments
  if(argc != 6)
    {
      printf("Usage:\n   hypertraps-dt.ce [observations-file] [random number seed] [length index] [kernel index] [considering losses (1) or gains (0)]\n\n");
      return 0;
    }
  seed = atoi(argv[2]);
  lengthindex = atoi(argv[3]);
  kernelindex = atoi(argv[4]);
  losses = atoi(argv[5]);

  printf("Running with:\n[observations-file]: %s\n[random number seed]: %i\n[length index]: %i\n[kernel index]: %i\n[losses (1) or gains (0)]: %i\n", argv[1], seed, lengthindex, kernelindex, losses);

  // initialise and allocate
  maxt = (lengthindex == 1 ? 1000 : (lengthindex == 2 ? 10000 : lengthindex == 3 ? 100000 : (lengthindex == 4 ? 1000000 : 100)));
  if(maxt <= 10000) SAMPLE = 100; else SAMPLE = 1000;

  if(_EVERYITERATION)
    SAMPLE = 1;

  srand48(121+seed);
  matrix = (int*)malloc(sizeof(int)*100000);

  // choose parameterisation based on command line
  expt = kernelindex;
  switch(expt)
    {
    case 0: DELTA = 0; break;
    case 1: DELTA = 0.005; MU = 0.1; break;
    case 2: DELTA = 0.05; MU = 1.; break;
    case 3: DELTA = 0.05; MU = 1.; break;
    case 4: DELTA = 0.1; MU = 1.; break;
    case 5: DELTA = 0.25; MU = 1.; break;
    case 6: DELTA = 0.5; MU = 1.; break;
    default: DELTA = 0.75; MU = 1.; break;
    }
  
  // read data on changes from input file
  // if we're thinking about losses, we're regarding gene losses as feature acquisitions; and thus inverting the data
  fp = fopen(argv[1], "r");
  if(fp == NULL)
    {
      printf("Couldn't find observations file %s\n", argv[1]);
      return 0;
    }
  i = 0; len = 0;
  do{
    ch = fgetc(fp);
    switch(ch)
      {
      case '0': matrix[i++] = (losses == 1 ? 1 : 0); break;
      case '1': matrix[i++] = (losses == 1 ? 0 : 1); break;
      case '2': matrix[i++] = 2; break;
      case '\n': if(len == 0) len = i; break;
      }
  }while(!feof(fp));
  ntarg = i/len;
  NVAL = len*(len+1);
  fclose(fp);

  ntau = ntarg/2;
  
  for(i = 0; i < ntau; i++)
    {
      tau1s[i] = 0;
      tau2s[i] = INFINITY;
    }

  printf("Observed transitions:\n");
  for(i = 0; i < ntarg/2; i++)
    {
      for(j = 0; j < len; j++) printf("%i", matrix[2*len*i+j]);
      printf(" -> ");
      for(j = 0; j < len; j++) printf("%i", matrix[2*len*i+len+j]);
      printf("\n");
    }
  if(losses == 1) printf("(where 1 is absence)\n\n");
  if(losses == 0) printf("(where 1 is presence)\n\n");
  
  printf("Number of features is %i, I found %i observation pairs and made up %i (0-inf) ghost timings\n", len, ntarg/2, ntau);
  printf("\n");
  
  // allocate memory and initialise output file
  trans = (double*)malloc(sizeof(double)*NVAL); 
  ntrans = (double*)malloc(sizeof(double)*NVAL);
  tmpmat = (double*)malloc(sizeof(double)*NVAL);

  // prepare output files
  sprintf(shotstr, "%s-posterior-%i-%i-%i-%i.txt", argv[1], spectrumtype, seed, lengthindex, kernelindex);
  fp = fopen(shotstr, "w"); fclose(fp);
  sprintf(bestshotstr, "%s-best-%i-%i-%i-%i.txt", argv[1], spectrumtype, seed, lengthindex, kernelindex);
  fp = fopen(bestshotstr, "w"); fclose(fp);

  // initialise with an agnostic transition matrix
  for(i = 0; i < len; i++)
    trans[i] = 1;
  for(i = len; i < len*(len+1); i++)
    trans[i] = 0;

  // compute initial likelihood given this matrix
  lik = GetLikelihoodCoalescentChange(matrix, len, ntarg, trans, parents, tau1s, tau2s);

  // initialise counters for acceptance ratio
  acc = rej = 0;
  lacc = lrej = 0;

  // run the chain
  for(t = 0; t < maxt; t++)
    {
      if(t % SAMPLE == 0)
	{
	  // periodically output progress to a tracker file
	  time(&timer);
	  tm_info = localtime(&timer);

	  strftime(buffer, 25, "%Y:%m:%d %H:%M:%S", tm_info);

	  fp = fopen("alltrackernew.txt", "a");
	  fprintf(fp, "%s %i %i %i %s\n", shotstr, t, maxt/5, maxt, buffer);
	  fclose(fp);
	}
      if(lik > bestlik || t == 0)
	{
	  bestlik = lik;
	  fp = fopen(bestshotstr, "w");
	  for(i = 0; i < len; i++)
	    {
	      for(j = 0; j < len; j++)
		{
		  //  ntrans[LEN+j*LEN+i] is the modifier for i from j*/
		  if(i == j) { fprintf(fp, "%f ", trans[i]); }
		  else { fprintf(fp, "%f ", trans[len+j*len+i]); }
		}
	    }
	  fprintf(fp, "\n");
	  fclose(fp);
	}

      // output some info periodically
      if(t % SAMPLE == 0)
	printf("%i - ", t);

      if(t > maxt/5 && t % SAMPLE == 0)
	{
	  // if we're burnt in, periodically sample the current parameterisation to an output file
	  fp = fopen(shotstr, "a");
	  for(i = 0; i < len; i++)
	    {
	      for(j = 0; j < len; j++)
		{
		  //  ntrans[LEN+j*LEN+i] is the modifier for i from j*/
		  if(i == j) { fprintf(fp, "%f ", trans[i]); }
		  else { fprintf(fp, "%f ", trans[len+j*len+i]); }
		}
	    }
	  fprintf(fp, "\n");
	  fclose(fp);
	}

      // apply a perturbation to the existing parameterisation
      // non-uniform priors can be employed here if desired 
      for(i = 0; i < NVAL; i++)
	{
	  ntrans[i] = trans[i];
	  r = RND;
	  if(r < MU)
	    {
	      ntrans[i] += gsl_ran_gaussian(DELTA);
	    }
	  if(ntrans[i] < -10) ntrans[i] = -10;
	  if(ntrans[i] > 10) ntrans[i] = 10;
	}

      // compute likelihood for the new parameterisation
      nlik = GetLikelihoodCoalescentChange(matrix, len, ntarg, ntrans, parents, tau1s, tau2s);

      // keep track of NaNs in calculations
      if(isnan(nlik))
	{
	  nancount++;
	}

      // compare likelihood to previous
      if(nlik >= lik || -(lik-nlik) > log(RND))
	{
	  // accept this new parameterisation
	  acc++; lacc++;
	  //	      if(t % SAMPLE == 0)
	  //	printf("acc ");
	  lik = nlik;
	  for(i = 0; i < NVAL; i++)
	    trans[i] = ntrans[i];
	}
      else 
	{
	  // reject the change
	  rej++; lrej++;
	  // if(t % SAMPLE == 0)
	  //	printf("rej ");
	}

      //      if(t % SAMPLE == 0) printf("NaN count %i of %i\n", nancount, t);

      // output information periodically
      if(t % TMODULE == 0)
	{
	  printf("Iteration %i likelihood %f total-acceptance %f recent-acceptance %f trial-likelihood %f\n", t, lik, acc/(acc+rej), lacc/(lacc+lrej), nlik);
	  lacc = lrej = 0;
	}
    }

  return 0;
}
 

