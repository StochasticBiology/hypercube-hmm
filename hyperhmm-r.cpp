#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Use the Rcpp namespace to access Rcpp functions
using namespace Rcpp;

// Use the arma namespace to access Armadillo functions
using namespace arma;

#define _USE_CODE_FOR_R 1
#include "hyperhmm.cpp"

List HyperHMM(NumericMatrix obs,
	      NumericVector seed,
	      NumericVector nboot,
	      NumericVector fullsample,
	      NumericVector longitudinal);

// [[Rcpp::export]]
List HyperHMM(NumericMatrix obs,
	      NumericVector seed = 1,
	      NumericVector nboot = 100,
	      NumericVector fullsample = 1,
	      NumericVector longitudinal = 0)
{
  int _longitudinal = longitudinal[0];
  int _fullsample = fullsample[0];
  int _nboot = nboot[0];
  int _seed = seed[0];
  double time;
  List infout;
  int L;

  L = obs.ncol();
  if(_longitudinal == 1)
    {
      vector<string> data;
      vector<int> data_count;
      for(int i = 0; i < obs.nrow(); i++)
	{
	  NumericVector this_row = obs(i, _);
	  CharacterVector rowAsString = as<CharacterVector>(wrap(this_row));
	  data.push_back(as<string>(rowAsString));
	  cout << rowAsString;
	  data_count.push_back(1);
	}

      infout = run_inference_longitudinal(data, data_count, L, _nboot, "", time, _fullsample);
    }
  else
    {
      vector<string> data;
      for(int i = 0; i < obs.nrow(); i++)
	{
	  NumericVector this_row = obs(i, _);
	  CharacterVector rowAsString = as<CharacterVector>(wrap(this_row));
 	  data.push_back(as<string>(rowAsString));
	  cout << rowAsString;
	}
      infout = run_inference(data, L, _nboot, "", time, _fullsample);
    }
  return infout;
}


 
