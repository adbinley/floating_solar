//input data

data {
  
  int<lower=0> N; //number of observations
  int<lower=1> n_order; //number of orders
  
  array[N] real log_axial_length; //vector of observations of axial length
  array[N] real log_bodymass; //vector of observations of bodymass (log transformed)
  array[N] int order; //vector of observations of order
  
}

// parameters
parameters {
  
  //real intercept;
  real slope;
  row_vector[n_order] order_intercept;
  real<lower = 0> sigma;
  
}

// The model to be estimated
model {
  
  //intercept ~ std_normal();
  order_intercept ~ std_normal();
  slope ~ std_normal();
  sigma ~ exponential(1);
  
  for (i in 1:N){
    
    log_axial_length[i] ~ normal(order_intercept[order[i]] + slope*log_bodymass[i], sigma);
  
  }

}

