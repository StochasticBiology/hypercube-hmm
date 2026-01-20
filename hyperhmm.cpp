#include <armadillo>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>
#include <cmath>
#include <string>
#include <iterator>
#include <random>

std::mt19937_64 rng;

using namespace std;


void myexit(int code)
{
#ifndef _USE_CODE_FOR_R
  exit(0);
#else
  Rcpp::stop("exiting");
#endif
}

/*
  A function for finding 2^r
  Input variables:
  - int r: the power of 2
  Output:
  - 2^r
*/
long int mypow2(int r)
{
  long int v = 1; //Store the output
  int i;
  for(i = 0; i < r; i++)
    v *= 2; //Multiply v with 2 r times
  return v;
}


/*
  A function for converting a number into a binary string.
  Input variables:
  - int n: The number you want to convert from int to binary string
  - int L: The total length of the binary string. This number needs to be greater then log_2(n)
  Output:
  - A string of length L representing the integer n as a binary string
*/
string number2binary(int n, int L){
  string binary = ""; //Store the output value

  //Creat thebinary number (without added 0's to the left)
  for(int v = mypow2(L-1); v >= 1; v /= 2){
    if(n >= v){
      n -= v;
      binary = binary + "1";
    }else{
      binary = binary + "0";
    }
  }
  //Add 0's to the left such that it has length L
  while(binary.length() < L){
    binary = "0" + binary;
  }
  return binary;
}



/*
  A function for converting a binary string into a number.
  Input variables:
  - string bin: The binary string you want to convert to a number
  - int L: The total length of the binary string.
  Output:
  - An integer form of the binary string
*/
int binary2int(string bin, int length){
  int number = 0; //Store the output value

  //Loop through the binary number
  for(int i=0; i<length; i++){

    //Check if the i'th place is a 1
    if(bin[i] == '1'){
      number += mypow2(length-i-1); //Add the corresponding value to the output
    }
  }
  return number;
}




/*
  Function for finding the r'th row of pascals triangle.
  Input variables:
  - int r: The row you want to find the values of
  Output:
  vector<int>: Containing all the elements of the r'th row
*/
vector<int> r_row_pascal(int r){
  vector<int> row;
  row.push_back(1);

  if (r == 0){
    return row;
  }

  vector<int> prev = r_row_pascal(r - 1);

  for(int i = 1; i < prev.size(); i++){
    int curr = prev[i - 1] + prev[i];
    row.push_back(curr);
  }
  row.push_back(1);
  return row;
}




int row_col_to_idx(int r, int c, arma::vec row_ptr, arma::vec col_idx){
  int idx = row_ptr(r);
  for(int i=row_ptr(r); i<row_ptr(r+1); i++){
    if(col_idx(idx) == c){
      break;
    }else{
      idx ++;
    }
  }
  return idx;
}



/*
  A function for creating numerous lists with information about the possible states at any given time t.
  Input variables (All of the vectors needs to be empty as input):
  - vector<int> n_states: The number of possible states at time t, where t is given by the index to that number. (This will always be the L'th row of the Pascal triangle)
  - vector<int> cumulative_states: A vector where the element at index i represents the index where t=i starts in the vector possible
  - vector<string> states: A vector containing all the indices for states in correct order.
  - int L: The total length of the binary string.
  Output:
  - Updated versions of n_states, cumulative_states, and states.
*/
void states_at_time_t(vector<int>& n_states, vector<int>& cumulative_states, vector<int>& states, int L){
  int index = 0;
  vector<int> count;
  string binary;

  //Find the r'th row of Pascals triangle
  n_states = r_row_pascal(L);
  //Add teh first element
  cumulative_states.push_back(1);

  //Redefine count such that it is of the same lenght as n_states
  count.resize(n_states.size(), 0);

  //Find the correct values of cumulative_states
  for(int i = 1; i < n_states.size(); i++){
    cumulative_states.push_back(cumulative_states[i-1] + n_states[i]);
  }
  //Add the first element of states
  states[0] = 0;

  //Loop through the remaining states
  for(int j = 1; j<mypow2(L); j++){
    //Find binary string value
    binary = number2binary(j, L);
    int size = 0;
    //Check the how many 1's the string contains
    for(int k=0; k<L; k++){
      if(binary[k] == '1'){
	size ++;
      }
    }
    //Add j to the correct place
    states[cumulative_states[size-1]+count[size]] = j;
    //Increase count by 1 in the correct place
    count[size] ++;
  }
}




/*
  A function for creating numerous lists with information about the possible states that can reach a state i.
  Input variables (All of the vectors needs to be empty as input):
  - vector<int> n_from: The number of possible states that can reach state i.
  - vector<int> cumulative_from: A vector where the element at index i represents the index where we find the first element that can reach state i.
  - vector<string> from: A vector containing all the indices (in CRS format) for states in correct order.
  - vector<int> cor_row: A vector containing the corresponding row to the index in from.
  - arma::vec row_ptr: The transition matrix row_ptr.
  - arma::vec col_idx: The transition matrix col_idx.
  - int L: The total length of the binary string.
  Output:
  - Updated versions of n_from, cumulative_from, from, and cor_row.
*/
void possible_transitions_from(vector<int>& n_from, vector<int>& cumulative_from, vector<int>& from, vector<int>& cor_row, int L, arma::vec row_ptr, arma::vec col_idx){
  int index = 0;
  int end_node_int, value;
  double time = 0.;
  string end_node;
  // loop through all vertices using integers for convenience
  for(int i = 0; i < mypow2(L); i++){
    // get binary name of this vertex
    string vertex = number2binary(i, L);
    n_from.push_back(0);
    // loop through bits looking for 1's we can change to 0's
    for(int j = 0; j < L; j++){
      if(vertex[j] == '1'){
	// we've found a 1 that we can change to a 0
	// store the resulting string as the "n_from[i]"th partner for i
	end_node = vertex;
	end_node[j] = '0';
	//Find the int value of this binary string
	end_node_int = binary2int(end_node, L);
	//Find the index value in CRS format
	int idx = row_ptr(end_node_int);
	while(col_idx(idx) != i){
	  idx ++;
	}
	from.push_back(idx);
	cor_row.push_back(end_node_int);
	// increment number of possible for i
	n_from[i]++;
      }
    }
    //Add the correct number to the cumulatice vector
    if(i==0){
      cumulative_from.push_back(0);
    }else{
      int c = cumulative_from[i-1] + n_from[i-1];
      cumulative_from.push_back(c);
    }
  }

}










/*
  A function for creating numerous lists with information about the possible states to go to from a given state.
  Input variables (All of the vectors needs to be empty as input):
  - vector<int> n_partners: The number of possible states to go to from the state represented by the index.
  - vector<int> cumulative_partners: A vector where the element at index i represents the index where t=i starts in the vector possible
  - vector<string> partners: A vector containing all the indices (in CRS format) for states in correct order.
  - int L: The total length of the binary string.
  Output:
  - Updated versions of n_partners, cumulative_partners, and partners.
*/
void possible_transitions(vector<int>& n_partners, vector<int>& cumulative_partners, vector<int>& partners, int L){
  int index = 0;
  int end_node_int;
  string end_node;
  //Loop through all the states
  for(int i = 0; i < mypow2(L); i++){
    string vertex = number2binary(i, L);
    n_partners.push_back(0);
    for(int j = 0; j < L; j++){
      if(vertex[j] == '0'){
	//If we find a 0 in the state it is possible to add a 1. Hence it is possible to go to that state with all equal elements except this place.
	end_node = vertex;
	end_node[j] = '1';
	end_node_int = binary2int(end_node, L);
	partners.push_back(end_node_int);
	n_partners[i]++;
      }
    }
    //Add the correct number to the cumulative vector
    if(i==0){
      cumulative_partners.push_back(0);
    }else{
      int c = cumulative_partners[i-1] + n_partners[i-1];
      cumulative_partners.push_back(c);
    }
  }
}





/*
  The function for calculating all the forward/alpha probabilities.
  Input variables:
  - arma::mat alpha: An all zero matrix for storing the alpha calculations
  - arma::vec A_val: The current transition matrix values
  - arma::vec A_row_ptr: The row pointer vector for the transition matrix
  - arma::vec A_col_idx: The column index vector for the transition matrix
  - vector<string> O: The current observation containing every known state and ? for every unknown state
  - vector<int> n_from, c_from, from, and cor_row: vectors from the possible_transitions_from function
  - vector<int> n_states, c_states, and states: vectors from the states_at_time_t function
  Output:
  - An updated alpha matrix filled with all the values for alpha_t(i)
*/
void forward_prob(arma::mat& alpha, arma::vec A_val, arma::vec A_row_ptr, arma::vec A_col_idx, vector<string> O, vector<int> n_from, vector<int> c_from, vector<int> from, vector<int> cor_row, vector<int> n_states, vector<int> c_states, vector<int> states){
  int T = O.size();
  int n = mypow2(T-1);
  double tmp = 0.;
  int start_state;

  //Add the initial state value
  //alpha(0,0) = 1.;

  int t_start = 0;
  for(int i=0; i<T-1; i++){
    string time_i = O[i];
    if(time_i.length() > T-3){
      t_start = i;
      start_state = binary2int(time_i, T-1);
      break;
    }
  }
  alpha(t_start,start_state) = 1.;

  //Loop thorugh all the possible times except the first.
  for(int t=t_start+1; t<T; t++){
    string time = O[t];
    if(time.length() < T-1){ //check if the state is unknown (a '?') at the given time t
      //Loop through all the states to be in at time t
      int start_i = c_states[t-1];
      for(int i=0; i<n_states[t]; i++){
	int i2 = states[start_i+i];
	int start_idx = c_from[i2];
	tmp = 0.;
	//Loop through all states where it is possible to reach state i
	for(int j=0; j<n_from[i2]; j++){
	  int j2 = from[start_idx+j];
	  tmp += alpha(t-1,cor_row[start_idx+j])*A_val(j2);
	}
	alpha(t,i2) = tmp;
      }
    }else{ //If the state is not '?' it is known.
      string s = O[t];
      int state = binary2int(s, T-1);
      int start_idx = c_from[state];
      //Loop through all states where it is possible to reach state "state"
      for(int j=0; j<n_from[state]; j++){
	int j2 = from[start_idx+j];
	alpha(t, state) += alpha(t-1,cor_row[start_idx+j])*A_val(j2);
      }
      /*
	if(alpha(t,state)==0){
	cout << "alpha(t,state): " << alpha(t,state) << endl;
	cout << t << " " << s << endl;
	for(int g=0; g<O.size();g++){
	cout << O[g] << endl;
	}
	int start_idx = c_from[state];
	for(int j=0; j<n_from[state]; j++){
	int j2 = from[start_idx+j];
	}
	}
      */
    }
  }
}



/*
  The function for calculating all the backwards/beta probabilities.
  Input variables:
  - arma::mat beta: An all zero matrix for storing the beta calculations
  - arma::vec A_val: The current transition matrix values
  - arma::vec A_row_ptr: The row pointer vector for the transition matrix
  - arma::vec A_col_idx: The column index vector for the transition matrix
  - vector<string> O: The current observation containing every known state and ? for every unknown state
  - vector<int> n_partners, c_partners, and partners: vectors from the possible_transitions function
  - vector<int> n_from, c_from, from, and cor_row: vectors from the possible_transitions_from function
  - vector<int> n_states, c_states, and states: vectors from the states_at_time_t function
  Output:
  - An updated beta matrix filled with all the values for beta_t(i)
*/
void backward_prob(arma::mat& beta, arma::vec A_val, arma::vec A_row_ptr, arma::vec A_col_idx, vector<string> O, vector<int> n_partners, vector<int> c_partners, vector<int> partners, vector<int> n_from, vector<int> c_from, vector<int> from, vector<int> cor_row, vector<int> n_states, vector<int> c_states, vector<int> states){
  int T = O.size();
  int n = mypow2(T-1);

  double time2 = 0.;

  //Add all values to the last row
  for(int i=0;i<n;i++){
    beta(T-1, i) = 1.;
  }

  //Loop through all times except the last. In decending order.
  for(int t = T-2; t > 0; t--){
    string time = O[t+1];

    if(time.length() < T-1){ //check if the state is unknown (a '?') at the given time t
      //Loop through all the states at time t
      int start_i = c_states[t-1];
      for(int i=0; i<n_states[t]; i++){
	int i2 = states[start_i+i];
	int start_idx = c_partners[i2];
	//Loop through all the possible states to go to from state i.
	for(int j = A_row_ptr[i2]; j<A_row_ptr[i2+1]; j++){
	  beta(t,i2) += A_val(j) * beta(t+1, A_col_idx(j));
	}
      }
    }else{ //If the state is not '?' it is known.
      int state = binary2int(time, T-1);
      int start_idx = c_from[state];
      for(int i=0; i<n_from[state]; i++){
	int i2 = from[start_idx+i];
	beta(t, cor_row[start_idx+i]) = A_val(i2) * beta(t+1, state);
      }
    }
  }
}


/*
  The function for calculating all the ksi probabilities.
  Input variables:
  - arma::mat ksi: The current ksi matrix (as given in the adapted_baum_welch function)
  - arma::mat alpha: The output matrix from forward_prob
  - arma::mat beta: The output matrix from backward_prob
  - arma::vec A_val: The current transition matrix values
  - arma::vec A_row_ptr: The row pointer vector for the transition matrix
  - arma::vec A_col_idx: The column index vector for the transition matrix
  - vector<string> O: The current observation containing every known state and ? for every unknown state
  - int n_O: The number of times this observation sequence appears
  - vector<int> n_partners, c_partners, and partners: vectors from the possible_transitions function
  - vector<int> n_from, c_from, from, and cor_row: vectors from the possible_transitions_from function
  - vector<int> n_states, c_states, and states: vectors from the states_at_time_t function
  Output:
  - An updated version of the ksi matrix
*/
void ksi_prob(arma::vec& ksi, arma::mat alpha, arma::mat beta, arma::vec A_val, arma::vec A_row_ptr, arma::vec A_col_idx, vector<string> O, int n_O, vector<int> n_partners, vector<int> c_partners, vector<int> partners, vector<int> n_from, vector<int> c_from, vector<int> from, vector<int> cor_row, vector<int> n_states, vector<int> c_states, vector<int> states){
  int T = O.size();
  int n = mypow2(T-1);
  double prob_obs = 0.;
  bool begin_calc = false;
  int idx;
  double prod;

  double time = 0.;

  prob_obs = alpha(T-1, n-1);

  double prod1 = n_O/prob_obs;

  /*
    if(isnan(prob_obs)){
    cout << "Prob_obs: " << prob_obs << endl;
    //for(int g=0; g<O.size();g++){
    //  cout << O[g] << endl;
    //}
    }
    if(isnan(prod1)){
    cout << "Prod1: " << prod1 << endl;
    }
  */
  int t_start = 0;
  for(int i=0; i<T-1; i++){
    string time_i = O[i];
    if(time_i.length() > T-3){
      t_start = i;
      break;
    }
  }


  //Loop thorugh all the times from 0 to T-1
  for(int t=t_start; t<T-1; t++){
    string time_i = O[t];
    string time_j = O[t+1];
    if(time_i.length() > T-3 && time_j.length() > T-3){ //Check if both of the observations are known.
      int i = binary2int(time_i, T-1);
      int j = binary2int(time_j, T-1);
      idx = A_row_ptr(i);
      while(A_col_idx(idx) != j){
	idx ++;
      }
      ksi(idx) += n_O * 1;
    }
    else if(time_i.length() > T-3){ //Check if the first observation is known.
      int i = binary2int(time_i, T-1);
      int start_idx = c_partners[i];
      prod = prod1 * alpha(t,i);
      for(int j = A_row_ptr[i]; j<A_row_ptr[i+1]; j++){
	ksi(j) += A_val(j)*beta(t+1,A_col_idx(j)) * prod;
      }
    }
    else if(time_j.length() > T-3){ //Check if the second observation is known
      int j = binary2int(time_j, T-1);
      prod = prod1 * beta(t+1,j);
      int start_idx = c_from[j];
      for(int i=0; i<n_from[j]; i++){
	int i2 = from[start_idx+i];
	ksi(i2) += prod * alpha(t,cor_row[start_idx+i])*A_val(i2);
      }
    }
    else{ //If none of the observations are known do this.
      int start_i = c_states[t-1];
      for(int i=0; i<n_states[t]; i++){
	int i2 = states[start_i+i];
	int start_idx = c_partners[i2];
	prod = prod1 * alpha(t,i2);
	for(int j = A_row_ptr[i2]; j<A_row_ptr[i2+1]; j++){
	  ksi(j) += prod * A_val(j)*beta(t+1,A_col_idx(j));
	}
      }
    }
  }
}






/*
  This is the main function of the algorithm. Here we do the iterations and find the maximum likelihood estimate of the transition matrix A.
  Input variables:
  - arma::vec A_val: The current transition matrix values
  - arma::vec A_row_ptr: The row pointer vector for the transition matrix
  - arma::vec A_col_idx: The column index vector for the transition matrix
  - vector<string> O: All the independent observations given as a flat vector
  - vector<int> n_O: A vector which states how many observation(s) we have of observation i
  - int max_itr: The maximum number of iterations we will allow before stopping
  - double eps: The convergence criteria
  - int L: The number of traits/length of the binary string
  - boool single_A: If true we will only look at a hypercube with one pathway (This is for testing cutting of edges)
  - boool double_A: If true we will only look at a hypercube with two possible pathways (This is for testing cutting of edges)
  Output:
  - The output is the maximum likelihood estimate of the transition matrix A
*/
void adapted_baum_welch(arma::vec& A_val, arma::vec A_row_ptr, arma::vec A_col_idx, vector<string> O, vector<int> n_O, int max_itr,int& itr, double eps, int L, bool single_A = false, bool double_A = false){
  double time_alpha = 0.;
  double time_beta = 0.;
  double time_ksi = 0.;
  double time_update = 0.;

  double time = 0.;

  int n = mypow2(L);
  //int itr = 0;
  int T = L +1;
  int total_obs = O.size()/(L+1);
  double change = eps;


  if(n_O[0] == 0){
    for(int i=0;i<O.size();i++){
      n_O[i] = 1;
    }
  }

  vector<int> n_from;
  vector<int> c_from;
  vector<int> from;
  vector<int> cor_row;

  vector<int> n_partners;
  vector<int> c_partners;
  vector<int> partners;

  vector<int> n_states;
  vector<int> c_states;
  vector<int> states;

  states.resize(mypow2(L), 0);

  states_at_time_t(n_states, c_states, states, L);
  possible_transitions_from(n_from, c_from, from, cor_row, L, A_row_ptr, A_col_idx);
  possible_transitions(n_partners, c_partners, partners, L);


  //Loop through either to the maximum iteration is reach or the transition matrix have converged enough
  while(max_itr > itr && change >= eps){
    //cout << "itr: " << itr << "\n";
    //Store a copy of the previous itr
    arma::vec A_prev = A_val;

    arma::vec ksi_sum(mypow2(L-1)*L, arma::fill::zeros);

    //Loop thorugh all of the observation sequences
    for(int i=0; i< total_obs; i++){
      arma::mat beta(T, n, arma::fill::zeros);
      arma::mat alpha(T, n, arma::fill::zeros);
      vector<string> o = std::vector<string>(O.begin() + i*(L+1), O.begin() + (i+1)*L + i+1);
      int n_o = n_O[i];

      auto t_alpha = std::chrono::high_resolution_clock::now();
      forward_prob(alpha, A_val, A_row_ptr, A_col_idx, o, n_from, c_from, from, cor_row, n_states, c_states, states);
      auto t_alpha_end = std::chrono::high_resolution_clock::now();
      double alpha_ti = std::chrono::duration<double>(t_alpha_end - t_alpha).count();
      time_alpha += alpha_ti;

      auto t_beta = std::chrono::high_resolution_clock::now();
      backward_prob(beta, A_val, A_row_ptr, A_col_idx, o, n_partners, c_partners, partners, n_from, c_from, from, cor_row, n_states, c_states, states);
      auto t_beta_end = std::chrono::high_resolution_clock::now();
      double beta_ti = std::chrono::duration<double>(t_beta_end - t_beta).count();
      time_beta += beta_ti;

      auto t_ksi = std::chrono::high_resolution_clock::now();
      ksi_prob(ksi_sum, alpha, beta, A_val, A_row_ptr, A_col_idx, o, n_o, n_partners, c_partners, partners, n_from, c_from, from, cor_row, n_states, c_states, states);
      auto t_ksi_end = std::chrono::high_resolution_clock::now();
      double ksi_ti = std::chrono::duration<double>(t_ksi_end - t_ksi).count();
      time_ksi += ksi_ti;
    }

    int at_partner = 0;

    auto t_update = std::chrono::high_resolution_clock::now();
    //Loop thorugh and update all the transition probabilities
    for(int k=0; k<mypow2(L); k++){
      int n_end_vertices = n_partners[k];
      double ksi_sum2 = 0.;
      int r = A_row_ptr(k);
      for(int i=0; i<A_row_ptr(k+1)-r; i++){
	ksi_sum2 += ksi_sum(r+i);
      }
      int r2 = A_row_ptr(k+1)-r;
      double st = 0.;
      int u=0;
      for(int l=A_row_ptr(k); l<A_row_ptr(k+1); l++){
	if(ksi_sum2 == 0){
	  A_val(l) = 1./r2;
	}else{
	  A_val(l) = ksi_sum(l) / ksi_sum2;
	}
	at_partner ++;
      }
    }





    auto t_update_end = std::chrono::high_resolution_clock::now();
    double update_ti = std::chrono::duration<double>(t_update_end - t_update).count();
    time_update += update_ti;

    itr ++;
    arma::mat change_mat = abs(A_val-A_prev);
    change = change_mat.max();
  }

  //cout << "This is time_alpha: " << time_alpha << endl;
  //cout << "This is time_beta: " << time_beta << endl;
  //cout << "This is time_ksi: " << time_ksi << endl;
  //cout << "This is time_update: " << time_update << endl;
}





//This function is just printed from this website: https://sort-care.github.io/C++-Random-Number/
double random_zero_to_one(){

  //  uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  // std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
  // rng.seed(ss);
  std::uniform_real_distribution<double> unif(0, 1);
  const int nSimulations = 10;
  double currentRandomNumber = unif(rng);
  return currentRandomNumber;
}


/*
  A function for finding the binary length og an interger.
  Input variables:
  - int n: The interger n you want to find the binary length of
  Output:
  - The binary length of the interger n
*/
int num2binlength(int n){
  string binary = "";
  for(int i=0; n>0; i++){
    binary += n%2;
    n = n/2;
  }
  return binary.length();
}


/*
  A function for making random walkers move on a transition matrix A.
  Input variables:
  - arma::mat rw: A matrix for storing the results of the ransom walkers
  - int n_walkers: The number of ransom walkers
  - arma::vec A_val: The transition matrix values
  - arma::vec A_row_ptr: The row pointer vector for the transition matrix
  - arma::vec A_col_idx: The column index vector for the transition matrix
  - int n_traits: The number of traits in the system (lenght of binary string)
  Output:
  - The matrix rw containing the probabilities of gaining trait i at time t
*/
void random_walkers(arma::mat& rw, int n_walkers, arma::vec A_val, arma::vec A_row_ptr, arma::vec A_col_idx, int n_traits){
  vector<int> n_partners;
  vector<int> c_partners;
  vector<int> partners;

  possible_transitions(n_partners, c_partners, partners, n_traits);
  double time2 = 0.;
  //Loop through the number of random walkers
  for(int i=0; i<n_walkers; i++){
    int t = 0;
    int state = 0;
    //Move from the first state (0) and stop once we have gained all the traits
    while(t < n_traits){
      double random_choice = random_zero_to_one();
      int choice = 0;
      int start_idx = c_partners[state];
      double prob = A_val(row_col_to_idx(state, partners[start_idx], A_row_ptr, A_col_idx));
      //Loop until we reach the random choice
      while(random_choice > prob){
	choice += 1;
	prob += A_val(row_col_to_idx(state, partners[start_idx+choice], A_row_ptr, A_col_idx));
      }
      int change = num2binlength(partners[start_idx+choice]-state)-1;
      rw(change, t) += 1;
      t += 1;
      state = partners[start_idx+choice];
    }

  }
  rw = rw/n_walkers;
}



/*
  A function for making a unform transiton matrix given a hypercube with a given dimension.
  Input variables:
  - arma::vec A_val: The vector for storing the transition matrix values
  - arma::vec A_row_ptr: The row pointer vector for the transition matrix
  - arma::vec A_col_idx: The column index vector for the transition matrix
  - int n_traits: The number of traits in the system
  Output:
  - The updated version of A with uniform transiton from each node
*/
void uniform_transition_matrix(arma::vec& A_val, arma::vec& A_row_ptr, arma::vec& A_col_idx, int n_traits){
  vector<int> n_partners;
  vector<int> c_partners;
  vector<int> partners;
  possible_transitions(n_partners, c_partners, partners, n_traits);
  int k = 0, c=0;
  //Loop through all the transition probabilities (>0) and update the spesific values to be uniform
  for(int i=0; i<mypow2(n_traits); i++){
    int n_end_vertices = n_partners[i];
    int r = A_row_ptr(i);
    A_row_ptr(i+1) = n_end_vertices + r;
    c = 0;
    for(int j=0; j<n_end_vertices; j++){
      int j2 = partners[k];
      A_val(r+c) = 1./n_end_vertices;
      A_col_idx(r+c) = j2;
      k ++, c++;
    }
  }
}






/*
  A function for importing cross-sectional data.
  Input variables:
  - string file_name: The full name of the .txt file
  - vector<string> data: An empty vector for storing the data
  Output:
  All the cross-sectional data in the "data" vector
*/
void import_data(string file_name, vector<string>& data, int *len){
  std::ifstream in(file_name);
  if(!in.is_open())
    {
      cout << "Can't find input file: " << file_name << "!\n";
      myexit(0);
    }
  std::string str;
  while(std::getline(in, str)){
    if(str.size()>0){
      (*len) = str.size();
      data.push_back(str);
    }
  }
  in.close();
}

/*
  A function for counting the number of 1's in a binary string.
  Input variables:
  - string binary: The binary string which we want to count the number of 1's in
  Output:
  int count: The number of 1's in the binary string
*/
int count_nr_1(string binary){
  int count = 0;
  for(int i=0; i<binary.size(); i++){
    if(binary[i] == '1'){
      count += 1;
    }
  }
  return count;
}


/*
  A function for importing longitudinal data. (This importing needs to be updated. How to get data_count to be as it should?)
  Input variables:
  - string file_name: The full name of the .txt file
  - vector<string> data: An empty vector for storing the data
  - vector<int> data_count: An empty vector for storing the count of each observation
  - int L: The number of traits
  Output:
  All the longitudinal data in the "data" vector and the count in data_count
*/
void import_data_longitudinal(string file_name, vector<string>& data, vector<int>& data_count, int *L){

  int c;
  string str;
  string start = "0";
  string end = "1";
  string word;

  std::ifstream in(file_name);
  if(!in.is_open())
    {
      cout << "Can't find input file: " << file_name << "!\n";
      myexit(0);
    }

  std::getline(in, str);
  stringstream spls(str);
  spls >> word;
  (*L) = word.size();
  in.seekg(0);

  for(int l=0; l<(*L)-1; l++){
    start += "0";
    end += "1";
  }
  while(std::getline(in, str)){
    int k = 0;
    if(str.size()>0){
      stringstream ss(str);
      string token;
      while(ss >> token){
	if(k==0){
	  c = count_nr_1(token);
	  if(c > 0){
	    data.push_back(start);
	  }
	  if(c != 1){
	    for(int j1=0; j1<c-1; j1++){
	      data.push_back("?");
	    }
	  }
	  data.push_back(token);
	  k++;
	}
	else if(k == 1){
	  int c2 = count_nr_1(token);
	  for(int j2=0; j2<c2-c-1; j2++){
	    data.push_back("?");
	  }
	  data.push_back(token);
	  if(c2 < token.size()){
	    for(int j3=c2; j3<(*L)-1; j3++){
	      data.push_back("?");
	    }
	    data.push_back(end);
	  }
	  k = 0;
	}
      }
      data_count.push_back(1);
    }
  }
  in.close();
}




/*
  A function for bootstraping cross-sectional data.
  Input variables:
  - vector<string> data: An empty vector for storing the data
  - vector<string> new_data: An empty vector for storing the bootstrap results
  Output:
  The updated new_data containing the bootstrap results
*/
void bootstrap(vector<string> data, vector<string>& new_data){
  int size = data.size();
  for(int i=0; i<size; i++){
    int random_idx = random_zero_to_one() * size;
    new_data[i] = data[random_idx];
  }
}


/*
  A function for bootstraping longitudinal data.
  Input variables:
  - vector<string> data: An empty vector for storing the data
  - vector<string> new_data: An empty vector for storing the bootstrap results
  - int L: The number of traits
  Output:
  The updated new_data containing the bootstrap results
*/
void bootstrap2(vector<string> data, vector<string>& new_data, int L){
  int size = data.size()/(L+1);
  for(int i=0; i<size-1; i++){
    int random_idx = random_zero_to_one() * size;
    for(int j=0; j<L+1; j++){
      new_data[i*(L+1)+j] = data[(random_idx)*(L+1)+j];
    }
  }
}


/*
  A function for bootstraping longitudinal data.
  Input variables:
  - vector<string> data: An empty vector for storing the data
  - vector<string> data_unique: An empty vector for storing the unique observations
  - vector<int> data_count: An empty vector for storing the count of each observation
  Output:
  The updated data_unique and data_count
*/
void unique_and_count(vector<string> data, vector<string>& data_unique, vector<int>& data_count){
  data_count.clear();
  data_unique.clear();
  sort(data.begin(), data.end());
  data_unique.push_back(data[0]);
  int count = 1;
  int n_unique = 1;
  for(int i=1; i<data.size(); i++){
    if(data[i] == data_unique[n_unique-1]){
      count ++;
      //if(i == data.size()-1){
      //  data_count.push_back(count);
      //}

    }else{
      data_unique.push_back(data[i]);
      data_count.push_back(count);
      count = 1;
      n_unique ++;
    }
  }
  data_count.push_back(count);

}



/*
  A function for adding questionmarks to the data.
  Input variables:
  - vector<string> data_unique: A vector storing the unique observations
  - vector<string> data_bw: An empty vector for storing the data with questionmarks
  - int L: The number of traits
  Output:
  The updated data_bw
*/
void add_questionmarks(vector<string> data_unique, int L, vector<string>& data_bw){
  data_bw.clear();
  string start = "0";
  string end = "1";
  for(int l=0; l<L-1; l++){
    start += "0";
    end += "1";
  }

  for(int i=0; i<data_unique.size(); i++){
    int count_1 = count_nr_1(data_unique[i]);
    if(count_1 > 0){
      data_bw.push_back(start);
    }
    for(int j1=0; j1<count_1-1; j1++){
      data_bw.push_back("?");
    }
    data_bw.push_back(data_unique[i]);
    for(int j2=0; j2<L-count_1-1; j2++){
      data_bw.push_back("?");
    }
    if(count_1 < L){
      data_bw.push_back(end);
    }
  }
}


/*
  A function for finding the mean matrix given a cude
  Input variables:
  - arma::cube C: The cube storing all the matrices
  - arma::mat mean_C: An all zero matrix
  - int n_slices: The number of matrices in the cube
  Output:
  Update of the mean_C with the mean results
*/
void mean_cube_slices(arma::cube C, arma::mat& mean_C, int n_slices){
  for(int i=0; i<n_slices; i++){
    mean_C += C.slice(i);
  }
  mean_C = mean_C/n_slices;
}


/*
  A function for finding the standard deviation matrix given a cude
  Input variables:
  - arma::cube C: The cube storing all the matrices
  - arma::mat sd_C: An all zero matrix
  - int n_slices: The number of matrices in the cube
  - arma::mat mean: The mean matrix of all the matrices in C
  Output:
  Update of the sd_C with the standard deviation results
*/
void sd_cube_slices(arma::cube C, arma::mat& sd_C, int n_slices, arma::mat mean){
  for(int i=0; i<n_slices; i++){
    sd_C += pow(C.slice(i) - mean, 2);
  }
  sd_C = sd_C/n_slices;
}


void graph_visualization(string file_name, int n_walkers, arma::vec A_val, arma::vec A_row_ptr, arma::vec A_col_idx, int n_traits){

  vector<int> n_partners;
  vector<int> c_partners;
  vector<int> partners;

  possible_transitions(n_partners, c_partners, partners, n_traits);
  std::ofstream myfile;
  myfile.open(file_name);
  //Loop through the number of random walkers
  for(int i=0; i<n_walkers; i++){
    int t = 0;
    int state = 0;
    //Move from the first state (0) and stop once we have gained all the traits
    while(t < n_traits){
      double random_choice = random_zero_to_one();
      int choice = 0;
      int start_idx = c_partners[state];
      double prob = A_val(row_col_to_idx(state, partners[start_idx], A_row_ptr, A_col_idx));
      //Loop until we reach the random choice
      while(random_choice > prob){
	choice += 1;
	prob += A_val(row_col_to_idx(state, partners[start_idx+choice], A_row_ptr, A_col_idx));
      }
      int change = num2binlength(partners[start_idx+choice]-state)-1;
      t += 1;
      int choice2 = partners[start_idx+choice];
      myfile << number2binary(state, n_traits) << " " << number2binary(choice2, n_traits) << endl;
      state = partners[start_idx+choice];
    }

  }
}


/*
  A function for doing the bootstrap
  Input variables:
  - string file_name: The full name of the .txt file
  - int n_boot: The total number of bootstraps you want to perform
  - int L: The numberof traits
  - arma::cube mean: ???
  - arma::mat sd: ???
  - string name: The name you want the output file to have
  Output:
  Two files containing the mean of the random walkers and the standard deviation
*/
#ifndef _USE_CODE_FOR_R
void run_inference(vector<string>& data, int L, int n_boot, string name, double& time, int rw_boot) {
  #else
List run_inference(vector<string>& data, int L, int n_boot, string name, double& time, int rw_boot) {
  #endif
    
  vector<string> data_bw;
    
  double time_itr = 0.;

  arma::cube mean(L,L,n_boot+1,arma::fill::zeros);
  arma::cube sd(L,L,n_boot+1,arma::fill::zeros);


  vector<string> data_unique;
  vector<int> data_count;

  unique_and_count(data, data_unique, data_count);
  add_questionmarks(data_unique, L, data_bw);


  arma::vec A_val(pow(2,L-1)*L, arma::fill::zeros);
  arma::vec A_row_ptr(pow(2,L)+1, arma::fill::zeros);
  arma::vec A_col_idx(pow(2,L-1)*L, arma::fill::zeros);
  uniform_transition_matrix(A_val, A_row_ptr, A_col_idx, L);

  int itr = 0;
  auto t1 = std::chrono::high_resolution_clock::now();
  adapted_baum_welch(A_val, A_row_ptr, A_col_idx, data_bw, data_count, 1000, itr, pow(10, -3), L, false, false);
  auto t2 = std::chrono::high_resolution_clock::now();
  double duration_seconds = std::chrono::duration<double>(t2 - t1).count(); //Measure time
  time_itr += duration_seconds/itr;
  time += duration_seconds;

  #ifndef _USE_CODE_FOR_R
  graph_visualization("graph_viz_" + name + ".txt", 10000, A_val, A_row_ptr, A_col_idx, L);
  #else
  #endif
  
  arma::mat rw(L,L, arma::fill::zeros);
  random_walkers(rw, 10000, A_val, A_row_ptr, A_col_idx, L);

  cout << "Mean ordering matrix from sampling:\n";
  cout << rw << endl;
  mean.slice(0) = rw;
  sd.slice(0) = rw;

  //Write the maximum likelihood values to file

  vector<int> n_partners;
  vector<int> c_partners;
  vector<int> partners;
  possible_transitions(n_partners, c_partners, partners, L);

  #ifndef _USE_CODE_FOR_R
  std::ofstream myfile3;
  myfile3.open("transitions_" + name + ".txt");
  myfile3 << "From " << "To " << "Probability" << endl;
  int k = 0, c=0;
  for(int i=0; i<mypow2(L); i++){
    int n_end_vertices = n_partners[i];
    int r = A_row_ptr(i);
    c = 0;
    for(int j=0; j<n_end_vertices; j++){
      int j2 = partners[k];
      myfile3 << i << " "<< j2 <<" "<< A_val(r+c) << endl;
      k ++, c++;
    }
  }
  #else
  NumericVector from_v, to_v, prob_v;
  int k = 0, c=0;
  for(int i=0; i<mypow2(L); i++){
    int n_end_vertices = n_partners[i];
    int r = A_row_ptr(i);
    c = 0;
    for(int j=0; j<n_end_vertices; j++){
      int j2 = partners[k];
      from_v.push_back(i);
      to_v.push_back(j2);
      prob_v.push_back(A_val(r+c));
      k ++, c++;
    }
  }
  #endif

  cout << "Running bootstrap resamples\n";
  int outboot;
  for(int i=0; i<n_boot; i++){
    outboot = 0;
    if(n_boot < 10) outboot = 1;
    else if(i % (int)(n_boot/10) == 0) outboot = 1;
    if(outboot == 1)
    {
      cout << "Bootstrap number: " << i+1 << "\n";
    }
    vector<string> new_data = data;
    bootstrap(data, new_data);
    unique_and_count(new_data, data_unique, data_count);
    add_questionmarks(data_unique, L, new_data);


    arma::vec A_val(pow(2,L-1)*L, arma::fill::zeros);
    arma::vec A_row_ptr(pow(2,L)+1, arma::fill::zeros);
    arma::vec A_col_idx(pow(2,L-1)*L, arma::fill::zeros);
    uniform_transition_matrix(A_val, A_row_ptr, A_col_idx, L);

    itr = 0;
    auto t3 = std::chrono::high_resolution_clock::now();
    adapted_baum_welch(A_val, A_row_ptr, A_col_idx, new_data, data_count, 1000, itr, pow(10, -3), L, false, false);
    auto t4 = std::chrono::high_resolution_clock::now();
    double duration_seconds2 = std::chrono::duration<double>(t4 - t3).count(); //Measure time
    time_itr += duration_seconds2/itr;
    time += duration_seconds2;
    if(rw_boot == 1){
      arma::mat rw(L,L, arma::fill::zeros);
      random_walkers(rw, 10000, A_val, A_row_ptr, A_col_idx, L);
      mean.slice(i+1) = rw;
      sd.slice(i+1) = rw;
    }
  }

  //cout << "Average rutime per iteration with " << n_boot << " bootstrap(s): " << time_itr/(n_boot+1) << endl;

  arma::mat mean_slices(L,L, arma::fill::zeros);
  mean_cube_slices(mean, mean_slices, n_boot+1);

  arma::mat sd_slices(L,L, arma::fill::zeros);
  sd_cube_slices(sd, sd_slices, n_boot+1, mean_slices);


    #ifndef _USE_CODE_FOR_R
  std::ofstream myfile;
  myfile.open("mean_" + name + ".txt");

  for(int m=0; m<L; m++){
    for(int n=0; n<L; n++){
      myfile << scientific << mean_slices(m,n) << " ";
    }
    myfile << "\n";
  }
  #else
  NumericMatrix mean_m(L, L);
  for(int m=0; m<L; m++){
    for(int n=0; n<L; n++){
      mean_m(m, n) = mean_slices(m,n);
    }
  }
  #endif

    #ifndef _USE_CODE_FOR_R
  std::ofstream myfile2;
  myfile2.open("sd_" + name + ".txt");

  for(int m=0; m<L; m++){
    for(int n=0; n<L; n++){
      myfile2 << scientific << sd_slices(m,n) << " ";
    }
    myfile2 << "\n";
  }
  #else
  NumericMatrix sd_m(L, L);
  for(int m=0; m<L; m++){
    for(int n=0; n<L; n++){
      sd_m(m, n) = sd_slices(m,n);
    }
  }
 List Lflux = List::create(Named("From") = from_v,
			    Named("To") = to_v,
			   Named("Probability") = prob_v);
  DataFrame Lfluxdf(Lflux);
  
  List Lout = List::create(Named("transitions") = Lfluxdf, Named("means") = mean_m, Named("sds") = sd_m);

  return Lout;
  #endif

}




#ifndef _USE_CODE_FOR_R
void run_inference_longitudinal(vector<string>& data, vector<int>& data_count, int L, int n_boot, string name, double& time, int rw_boot){
#else
List run_inference_longitudinal(vector<string>& data, vector<int>& data_count, int L, int n_boot, string name, double& time, int rw_boot){
#endif
  
  int itr = 0;
			     
  arma::cube mean(L,L,n_boot+1,arma::fill::zeros);
  arma::cube sd(L,L,n_boot+1,arma::fill::zeros);

  arma::vec A_val(pow(2,L-1)*L, arma::fill::zeros);
  arma::vec A_row_ptr(pow(2,L)+1, arma::fill::zeros);
  arma::vec A_col_idx(pow(2,L-1)*L, arma::fill::zeros);
  uniform_transition_matrix(A_val, A_row_ptr, A_col_idx, L);
  auto t1 = std::chrono::high_resolution_clock::now();
  adapted_baum_welch(A_val, A_row_ptr, A_col_idx, data, data_count, 1000, itr, pow(10, -3), L, false, false);
  auto t2 = std::chrono::high_resolution_clock::now();
  double duration_seconds = std::chrono::duration<double>(t2 - t1).count(); //Measure time
  time += duration_seconds;
  arma::mat rw(L,L, arma::fill::zeros);
  random_walkers(rw, 10000, A_val, A_row_ptr, A_col_idx, L);

  #ifndef _USE_CODE_FOR_R
  graph_visualization("graph_viz_" + name + ".txt", 10000, A_val, A_row_ptr, A_col_idx, L);
  #else
  #endif
  
  cout << "Mean ordering matrix from sampling:\n";
  cout << rw << endl;
  mean.slice(0) = rw;;
  sd.slice(0) = rw;

  //Write the maximum likelihood values to file
  vector<int> n_partners;
  vector<int> c_partners;
  vector<int> partners;
  possible_transitions(n_partners, c_partners, partners, L);

  #ifndef _USE_CODE_FOR_R
  std::ofstream myfile3;
  myfile3.open("transitions_" + name + ".txt");
  myfile3 << "From " << "To " << "Probability" << endl;
  int k = 0, c=0;
  for(int i=0; i<mypow2(L); i++){
    int n_end_vertices = n_partners[i];
    int r = A_row_ptr(i);
    c = 0;
    for(int j=0; j<n_end_vertices; j++){
      int j2 = partners[k];
      myfile3 << i << " "<< j2 <<" "<< A_val(r+c) << endl;
      k ++, c++;
    }
  }
  #else
  NumericVector from_v, to_v, prob_v;
  int k = 0, c=0;
  for(int i=0; i<mypow2(L); i++){
    int n_end_vertices = n_partners[i];
    int r = A_row_ptr(i);
    c = 0;
    for(int j=0; j<n_end_vertices; j++){
      int j2 = partners[k];
      from_v.push_back(i);
      to_v.push_back(j2);
      prob_v.push_back(A_val(r+c));
      k ++, c++;
    }
  }
  #endif

  int outboot;
  for(int i=0; i<n_boot; i++){
    outboot = 0;
    if(n_boot < 10) outboot = 1;
    else if(i % (int)(n_boot/10) == 0) outboot = 1;
    if(outboot == 1)
    {
      cout << "Bootstrap number: " << i+1 << "\n";
    }
    vector<string> new_data = data;
    bootstrap2(data, new_data, L);

    arma::vec A_val(pow(2,L-1)*L, arma::fill::zeros);
    arma::vec A_row_ptr(pow(2,L)+1, arma::fill::zeros);
    arma::vec A_col_idx(pow(2,L-1)*L, arma::fill::zeros);
    uniform_transition_matrix(A_val, A_row_ptr, A_col_idx, L);
    auto t3 = std::chrono::high_resolution_clock::now();
    adapted_baum_welch(A_val, A_row_ptr, A_col_idx, new_data, data_count, 1000,itr, pow(10, -3), L, false, false);
    auto t4 = std::chrono::high_resolution_clock::now();
    double duration_seconds2 = std::chrono::duration<double>(t4 - t3).count(); //Measure time
    time += duration_seconds2;
    if(rw_boot == 1){
      arma::mat rw(L,L, arma::fill::zeros);
      random_walkers(rw, 10000, A_val, A_row_ptr, A_col_idx, L);
      mean.slice(i+1) = rw;
      sd.slice(i+1) = rw;
    }

  }

  arma::mat mean_slices(L,L, arma::fill::zeros);
  mean_cube_slices(mean, mean_slices, n_boot+1);

  arma::mat sd_slices(L,L, arma::fill::zeros);
  sd_cube_slices(sd, sd_slices, n_boot+1, mean_slices);


  #ifndef _USE_CODE_FOR_R
  std::ofstream myfile;
  myfile.open("mean_" + name + ".txt");

  for(int m=0; m<L; m++){
    for(int n=0; n<L; n++){
      myfile << scientific << mean_slices(m,n) << " ";
    }
    myfile << "\n";
  }
 #else
  NumericMatrix mean_m(L, L);
  for(int m=0; m<L; m++){
    for(int n=0; n<L; n++){
      mean_m(m, n) = mean_slices(m,n);
    }
  }
  #endif

   #ifndef _USE_CODE_FOR_R
  std::ofstream myfile2;
  myfile2.open("sd_" + name + ".txt");

  for(int m=0; m<L; m++){
    for(int n=0; n<L; n++){
      myfile2 << scientific << sd_slices(m,n) << " ";
    }
    myfile2 << "\n";
  }
   #else
  NumericMatrix sd_m(L, L);
  for(int m=0; m<L; m++){
    for(int n=0; n<L; n++){
      sd_m(m, n) = sd_slices(m,n);
    }
  }
 List Lflux = List::create(Named("From") = from_v,
			    Named("To") = to_v,
			   Named("Probability") = prob_v);
  DataFrame Lfluxdf(Lflux);
  
  List Lout = List::create(Named("transitions") = Lfluxdf, Named("means") = mean_m, Named("sds") = sd_m);

  return Lout;
  #endif
}



/*
  A function for finding the sum of the elements in a vector of doubles.
  Input variables:
  - vector<double> vec: The vector to take the sum of
  Output:
  The sum of the vector
*/
double sum_vector(vector<double> vec){
  double s = 0.;
  for(int i=0; i<vec.size(); i++){
    s += vec[i];
  }
  return s;
}

void help(void)
{
  printf("Options [defaults]:\n\n--obs file.txt\t\tobservations file [NA]\n--label labelstring\tString to label file output [[filename]-out]\n--seed N\t\tRandom number seed [1]\n--nboot N\t\tNumber of bootstrap resamples [100]\n--fullsample\t\tSimulate walkers for each resample [0]\n--longitudinal\t\tLongitudinal data [OFF]\n--help\t\t\tDisplay this message\n\n");
}

/*
  The command line arguments needs to be as follow:
  data_file_name.txt L number_bootstraps name_output_file cross_sectional rw_bootstrap

  If you just want the maximum likelihood results set number_bootstrap = 0.
  If you have cross_sectional data set cross_sectional = 1, if longitudinal data set cross_sectional = 0.
  If you want to add summary data for each bootstrap set rw_bootstrap = 1, else set rw_bootstrap = 0.
*/
int main(int argc, char** argv){

  int L, n_boot, rw_boot;
  char obsfile[1000];
  char labelstr[1000];
  int cross_sectional;
  int filelabel;
  int fullsample;
  int seed;
  double time;
  int i;
  int oldform = 0;
  vector<string> data;
  vector<int> data_count;
  
  n_boot = 100;
  fullsample = 0;
  strcpy(obsfile, "");
  strcpy(labelstr, "");
  seed = 1;
  filelabel = 0;
  cross_sectional = 1;

  printf("== HyperHMM ==\n\nPlease cite Moen, M.T. and Johnston, I.G., 2023. HyperHMM: efficient inference of evolutionary and progressive dynamics on hypercubic transition graphs. Bioinformatics, 39(1), p.btac803.\n\n");

  if(argc == 7)
    {
      if(!(atoi(argv[5]) != 0 && atoi(argv[5]) != 1) || atoi(argv[2]) <= 0 || atoi(argv[3]) < 0 || (atoi(argv[6]) != 0 && atoi(argv[6]) != 1)) {
	oldform = 1;
      }
    }

  if(oldform == 1)
    {
      //    cout << "Usage:\n\t./hyperhmm.ce [datafile] [number of features] [number of bootstrap resamples] [output file label] [cross-sectional data (0 or 1)] [simulate random walkers for each sample (0 or 1)]\n";
      // 1 = datafile, 2 = L, 3 = nboot, 4 = output label, 5 = cross-sectional, 6 = walkers for each resample
      printf("It looks like you've used the old command-line argument specification, which is preserved for backwards compatibility. The newer format may allow more flexibility:\n");
      help();
      n_boot = atoi(argv[3]);
      rw_boot = atoi(argv[6]);
      cross_sectional = atoi(argv[5]);
      sprintf(obsfile, "%s", argv[1]);
      sprintf(labelstr, "%s", argv[4]);  
    }
  else
    {
      // deal with command-line arguments
	
      for(i = 1; i < argc; i+=2)
	{
	  if(strcmp(argv[i], "--obs\0") == 0) strcpy(obsfile, argv[i+1]);
	  else if(strcmp(argv[i], "--label\0") == 0) { filelabel = 1; strcpy(labelstr, argv[i+1]); }
	  else if(strcmp(argv[i], "--seed\0") == 0) seed = atoi(argv[i+1]);
	  else if(strcmp(argv[i], "--nboot\0") == 0) n_boot = atoi(argv[i+1]);
	  else if(strcmp(argv[i], "--fullsample\0") == 0) { rw_boot = 1; i--; }
	  else if(strcmp(argv[i], "--longitudinal\0") == 0) { cross_sectional = 0; i--; }
	  else if(strcmp(argv[i], "--help\0") == 0) { help(); return 0; }
	  else { printf("Didn't understand argument %s\n", argv[i]); i--; }
	}
      if(filelabel == 0)
	{
	  //sprintf(labelstr, "%s-out", obsfile);
	  snprintf(labelstr, sizeof(labelstr), "%.*s-out", (int)(sizeof(labelstr) - 5), obsfile);
	}
    }

  rng.seed(seed);
  
  if(strcmp(obsfile, "") == 0)
    {
      printf("I need at least an input file!\n\n");
      help();
      myexit(0);
    }
  
  // printf("%f %f %f %f\n", random_zero_to_one(), random_zero_to_one(), random_zero_to_one(), random_zero_to_one());

  printf("Attempting to run HyperHMM with the following parameters:\n");
  printf("  obs: %s\n", obsfile);
  printf("  cross-sectional: %i\n", cross_sectional);
  printf("  label: %s\n", labelstr);
  printf("  seed: %i\n", seed);
  printf("  nboot: %i\n", n_boot);
  printf("  fullsample: %i\n\n", rw_boot);
        
  if(cross_sectional == 1){
    import_data(obsfile, data, &L);
    cout << "From " << obsfile << " I read " << data.size() << " entries and am assuming " << L << " features\n\n";
    run_inference(data, L, n_boot, labelstr, time, rw_boot);
  }
  else if(cross_sectional == 0){
    import_data_longitudinal(obsfile, data, data_count, &L);
    cout << "From " << obsfile << " I read " << data.size() << " entries and am assuming " << L << " features\n\n";
    run_inference_longitudinal(data, data_count, L, n_boot, labelstr, time, rw_boot);
  }


  //cout << "The time it took to run the algorithm with " << n_boot << " bootstrap(s): " << time/(n_boot+1) << endl;
  //cout << "Average rutime per iteration with " << n_boot << " bootstrap(s): " << time/(n_boot+1) << endl;
  return 0;
}

