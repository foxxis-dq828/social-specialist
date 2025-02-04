# Load necessary library
library(dplyr)

# Define the Softmax Function
sf <- function(ll_values, tau) {
  # Prevent overflow by subtracting the max before exponentiation
  ll_values <- ll_values - max(ll_values)
  exp_neg_tau_ll <- exp(ll_values * tau)
  probabilities <- exp_neg_tau_ll / sum(exp_neg_tau_ll)
  return(probabilities)
}

#define alpha and beta
integral_function<-function(c,n1,i,tau){
  w_list <- c(0.5, 1-n1, 1-n1, 1-n1, n1, n1, 0.5, n1)
  w<-w_list[i]
  alpha= 1+w*c
  beta= 1+(1-w)*c
  ll1<-pbeta(1/6,alpha,beta)
  ll2<-pbeta(2/6,alpha,beta)-pbeta(1/6,alpha,beta)
  ll3<-pbeta(3/6,alpha,beta)-pbeta(2/6,alpha,beta)
  ll4<-pbeta(4/6,alpha,beta)-pbeta(3/6,alpha,beta)
  ll5<-pbeta(5/6,alpha,beta)-pbeta(4/6,alpha,beta)
  ll6<-pbeta(6/6,alpha,beta)-pbeta(5/6,alpha,beta)
  #create a list
  ll_values <- c(ll1,ll2,ll3,ll4,ll5,ll6)
  sf_list<-sf(ll_values,tau)
  return(sf_list)
}



#adjusting the training data
#condition1
choice_d1<-data%>%
  group_by(d1)%>%
  summarise(count = n())%>%
  rename(d=d1)%>%
  mutate(condition=1)
choice_d1
#condition2
choice_d2<-data%>%
  group_by(d2)%>%
  summarise(count = n())%>%
  rename(d=d2)%>%
  mutate(condition=2)
choice_d2
#condition3
choice_d3<-data%>%
  group_by(d3)%>%
  summarise(count = n())%>%
  rename(d=d3)%>%
  mutate(condition=3)
choice_d3
#condition4
choice_d4<-data%>%
  group_by(d4)%>%
  summarise(count = n())%>%
  rename(d=d4)%>%
  mutate(condition=4)
choice_d4
#condition5
choice_d5<-data%>%
  group_by(d5)%>%
  summarise(count = n())%>%
  rename(d=d5)%>%
  mutate(condition=5)
choice_d5
#condition6
choice_d6<-data%>%
  group_by(d6)%>%
  summarise(count = n())%>%
  rename(d=d6)%>%
  mutate(condition=6)
choice_d6
#condition7
choice_d7<-data%>%
  group_by(d7)%>%
  summarise(count = n())%>%
  rename(d=d7)%>%
  mutate(condition=7)
choice_d7
#condition8
choice_d8<-data%>%
  group_by(d8)%>%
  summarise(count = n())%>%
  rename(d=d8)%>%
  mutate(condition=8)
choice_d8
#join this condition dataframe
all_data<-list(choice_d1,choice_d2,choice_d3,choice_d4,choice_d5,choice_d6,choice_d7,choice_d8)
all_choice <- Reduce(function(x, y) full_join(x, y, by = intersect(names(x), names(y))), all_data)
all_choice

# Define the Objective Function for Optimization
objective_function <- function(par) {
  # Extract ce and n from the parameter vector (ensure positivity)
  c <- exp(par[1])
  n1 <- logistic(par[2])
  tau <- exp(par[3])
  
  # Initialize log-likelihood sum
  ll_log_sum <- 0
  
  # Iterate over each row in the combined choice data
  for (m in 1:nrow(all_choice)) {
    likert <- all_choice$d[m]
    condition <- all_choice$condition[m]
    count <- all_choice$count[m]
    sf_lllist <- integral_function(c,n1,condition,tau)
    # Get the probability of the observed Likert value
    ll_here <-  sf_lllist[likert]
    # Accumulate the negative log-likelihood
    ll_log_sum <- ll_log_sum - count * log(ll_here)
  }
  return(ll_log_sum)
}


# Perform Optimization Using L-BFGS-B Method
initial_params<-c(0,0,0)
result <- optim(
  par = initial_params,
  fn = objective_function,
)

c <- exp(result$par[1])
n1 <- logistic(result$par[2])
tau <- exp(result$par[3])

# View Optimization Results
print(result)

# recalculate the parameter
q <- function(n){
  (1+2*n)/(2+2*n) - n1}
res <- uniroot(q, interval = c(1, 10))
n <-res$root