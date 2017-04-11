// 01/05/2016
// State-Space model of inverts at Lees Ferry


data {
  int Nmd;                  // Number of drift months from start   
  int Nmb;                  // Number of benthic months from start
  
  int Nmonth;               // Number of months (12)
  
  int Nsamps;               // Number of drift samples 
  
  vector[Nmd] d_trip;       // Index for month from start (drift)
  int month[Nmd];        // Index for month in year (1 - 12)
  
  int DC[Nsamps, Nmd];   // Drift conc. 
  
  matrix[Nsamps, Nmd] log_Q; // log Q
} 

parameters {
  vector[Nmd] beta_o;
  // beta_b
  real beta_Q;
  vector[Nmonth] eta_month;
  real<lower = 0> sig_month;
}

transformed parameters {
  matrix[Nsamps, Nmd] lamda;
  
  for(j in 1:Nsamps){
     for(k in 1:Nmd){
       // lamda[j,k] = exp(beta_o + beta_b * log_tN[d_trip[k]] + beta_Q * log_Q[j,k] + eta);
       lamda[j,k] = exp(beta_o[k] + beta_Q * log_Q[j,k]);// + eta_month[month[k]]);
     }
   }
}

model {
  beta_o ~ normal(0, 10);
  sig_month ~ normal(0, 10);
  
  for(i in 1:Nmonth){
    eta_month[i] ~ normal(0, sig_month);
  }
  
   for(j in 1:Nsamps){
     for(k in 1:Nmd){
       DC[j,k] ~ poisson(lamda[j,k]);
     }
   }
}

// generated quantities {
// }
