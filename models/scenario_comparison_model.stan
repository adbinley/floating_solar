//input data

data {
  
  int<lower=0> N; //number of observations
  int<lower=1> n_WQ; //binary categories for water quality
  int<lower=1> n_SOC; //binary categories for social benefits
  
  
  array[N] real VI_value; //vector of observations of VI
  array[N] int WQ_value; //vector of observations of water quality (binary)
  array[N] int SOC_value; //vector of observations of social benefit (binary)
  
}

// parameters
parameters {
  
  //real intercept;
  row_vector[n_WQ] WQ_intercept;
  row_vector[n_SOC] SOC_intercept;
  //real slope;
  real<lower = 0> sigma;
  
}

// The model to be estimated
model {
  
  //intercept ~ std_normal();
  WQ_intercept ~ std_normal();
  SOC_intercept ~ std_normal();
  //slope ~ std_normal();
  sigma ~ exponential(1);
  
  for (i in 1:N){
    
    VI_value[i] ~ normal(WQ_intercept[WQ_value[i]] + SOC_intercept[SOC_value[i]], sigma);
  
  }

}

