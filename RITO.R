###### This code summarizes the functions and variables used in the RITO algorithm, taking the SSALT of
###### single stress voltage as an example, and provides comments on the main steps. You can adjust the 
###### corresponding parameters and environment according to your own needs.

##### The parameter initialization for ALT of Voltage
#### We have 3 groups
alpha=4
beta0=20.4
beta1=-2.74
x1=log(80)
x2=log(100)
x3=log(120)
theta1=exp(beta0+x1*beta1)
theta2=exp(beta0+x2*beta1)
theta3=exp(beta0+x3*beta1)

##### The decomposed trace function for Weibull CSALT
TR=function(t,theta,alpha,Ni,x){
  ### t: inspection times
  ### theta: scale paramter
  ### Ni: number of unit
  n=length(t)
  F=0
  p=c()
  ###interval
  P=c()
  for (i in 1:n) {
    p=c(p,pweibull(t[i],alpha,theta)-F)
    F=pweibull(t[i],alpha,theta)
    P=c(P,F)
  }
  p=c(p,1-F)
  Ai=alpha*(1-P)*log(1-P)
  Ndp=Ni/p
  coeff=c((1+x^2+1/(alpha^2)*log(t[1]/theta)^2)*Ai[1]^2)
  for (i in 2:n) {
    coe=(1+x^2)*(Ai[i]-Ai[i-1])^2+1/(alpha^2)*(log(t[i]/theta)*Ai[i]-log(t[i-1]/theta)*Ai[i-1])^2
    coeff=c(coeff,coe)
  }
  coeff=c(coeff,(1+x^2+1/(alpha^2)*log(t[i]/theta)^2)*Ai[i]^2)
  trace=sum(coeff*Ndp)
  return(-trace)
}

##### The decomposed trace function (conditional on the past inspections) for Weibull CSALT
TR_conditional=function(t_past,t_fu,theta,alpha,Ni,x){
  ### t_past: the past inspection time
  ### t_fu: the future inspection time
  n=length(t_past)+length(t_fu)
  t=c(t_past,t_fu)
  F=0
  p=c()
  ###interval
  P=c()
  for (i in 1:n) {
    p=c(p,pweibull(t[i],alpha,theta)-F)
    F=pweibull(t[i],alpha,theta)
    P=c(P,F)
  }
  p=c(p,1-F)
  Ai=alpha*(1-P)*log(1-P)
  Ndp=Ni/p
  coeff=c((1+x^2+1/(alpha^2)*log(t[1]/theta)^2)*Ai[1]^2)
  for (i in 2:n) {
    coe=(1+x^2)*(Ai[i]-Ai[i-1])^2+1/(alpha^2)*(log(t[i]/theta)*Ai[i]-log(t[i-1]/theta)*Ai[i-1])^2
    coeff=c(coeff,coe)
  }
  coeff=c(coeff,(1+x^2+1/(alpha^2)*log(t[i]/theta)^2)*Ai[i]^2)
  trace=sum(coeff*Ndp)
  return(-trace)
}

##### We can optimize the future inspections conditional on the past inspections with the paramters. like
optim(c(1500,2000), TR_conditional, x=x3,theta=theta3,t_past=c(500,1000,1200),Ni=10,method = "BFGS")
### With the trace decomposition, it it easy to use the unconstrained optimization solver optim()

##### The Weibull likelihood function for interval censoring
Weibull_likelihood_interval=function(parameter,num_of_failure,check_time){
  ### number_of_failure: the matirx with (i,j) element represents the number of unit
  ### in the jth interval of the ith group
  ### check_time: the matrix with (i,j) element represents the jth inspection time of ith group
  x=c(x1,x2,x3)
  alpha=parameter[1]
  beta0=parameter[2]
  beta1=parameter[3]
  theta1=exp(beta0+x1*beta1)
  theta2=exp(beta0+x2*beta1)
  theta3=exp(beta0+x3*beta1)
  theta=c(theta1,theta2,theta3)
  m=dim(num_of_failure)[1]
  n=dim(num_of_failure)[2]
  
  p=matrix(0,m,n)
  for (i in 1:m) {
    F=0
    for (j in 1:(n-1)) {
      p[i,j]=pweibull(check_time[i,j],alpha,theta[i])-F
      F=pweibull(check_time[i,j],alpha,theta[i])
    }
    p[i,n]=1-F
  }
  logL=0
  for (i in 1:m) {
    logL=logL+sum(num_of_failure[i,]*log(p[i,]))
  }
  return(-logL)
}


##### An example for Weibull MLE
check_time=matrix(c(1000,1500,5000,1000,1500,3000,1000,1500,2000),3,3,byrow=T)
num_of_failure=matrix(c(0,1,16,3,0,3,15,2,3,8,8,1),3,4,byrow=T)
check_time=matrix(c(1000,1500,1000,1500,1000,1500),3,2,byrow=T)
num_of_failure=matrix(c(0,1,19,0,2,18,3,8,9),3,3,byrow=T)
optim(c(3.5,20,-2),Weibull_likelihood_interval,num_of_failure=num_of_failure,check_time=check_time,method = "BFGS")


##### Function to find the next inspection time
nonzerominposition=function(matrix){
  non_zero_positions <- which(matrix != 0, arr.ind = TRUE)
  non_zero_elements <- matrix[non_zero_positions]
  min_non_zero <- min(non_zero_elements)
  position=which(matrix ==min_non_zero, arr.ind = TRUE)
  return(position)
}


theta=c(theta1,theta2,theta3)
x=c(x1,x2,x3)

##### Function to generate the exact failure time
vol_Failure_times=function(theta,alpha,num_in_group){
  failure_time=matrix(0,3,num_in_group)
  for (i in 1:3) {
    failure_time[i,]=sort(rweibull(num_in_group,alpha,theta[i]))
  }
  return(failure_time)
}

##### Function to convert the exact failure time matrix to matrix of number of failures in each interval
num_in_inspection=function(failure_time,inspection_time){
  n=dim(inspection_time)[2]
  num_matrix=matrix(0,3,n+1)
  for (i in 1:3) {
    num_matrix[i,1]=sum(failure_time[i,]<=inspection_time[i,1])
    for (j in 2:n) {
      num_matrix[i,j]=sum((failure_time[i,]>inspection_time[i,j-1])&(failure_time[i,]<=inspection_time[i,j]))
    }
    num_matrix[i,j+1]=sum(failure_time[i,]>inspection_time[i,j])
  }
  return(num_matrix)
}

##### The main function for RITO simulation
vol_trace_simulation=function(theta,x,alpha,num_in_group,num_of_inspection){
  failure_time=vol_Failure_times(theta,alpha,num_in_group)
  ### Failure time generation and initialization
  inspection_initial=matrix(c(1000,1500,1000,1500,1000,1500),3,2,byrow=T)
  num_matrix_initial=num_in_inspection(failure_time,inspection_initial)
  ### Comparision with fixed-time design
  fix_inspection=matrix(rep(500*(1:num_of_inspection),3),3,num_of_inspection,byrow = T)
  fix_num_matrix=num_in_inspection(failure_time,fix_inspection)
  para_fix=optim(c(3.5,20,-2),Weibull_likelihood_interval,num_of_failure=fix_num_matrix,check_time=fix_inspection,method = "BFGS")$par
  para_est=optim(c(3.5,20,-2),Weibull_likelihood_interval,num_of_failure=num_matrix_initial,check_time=inspection_initial,method = "BFGS")$par
  theta_est=exp(para_est[2]+para_est[3]*x)
  inspection_time=rep(0,num_of_inspection)
  for (i in 1:3) {
    inspection=c(1000,1500,optim(seq(4000/i,6000/i,2000/((num_of_inspection-3)*i)), TR_conditional, x=x[1],theta=theta_est[i],alpha=para_est[1],t_past=c(1000,1500),method = "BFGS")$par)
    inspection_time=rbind(inspection_time,inspection)
  }
  inspection_time=inspection_time[-1,]
  already_inspection_time=matrix(0,3,num_of_inspection)
  already_inspection_time[,c(1,2)]=matrix(c(1000,1500,1000,1500,1000,1500),3,2,byrow=T)
  ### Matrix to record the already recorded inspection times
  while (already_inspection_time[1,num_of_inspection]==0) {
    ### Condition is for convenient
    position=nonzerominposition(inspection_time-already_inspection_time)
    ### Next inspection
    if(inspection_time[position]>max(already_inspection_time)){
      already_inspection_time[position]=inspection_time[position]
      ### For the spectial case that after a inspection, the next suggested time is even more earlier
    }else{
      already_inspection_time[position]=max(already_inspection_time)
    }
    num_matrix=num_in_inspection(failure_time,inspection_time)
    ### Record the failures
    para_est=optim(c(3.5,20,-2),renew_Weibull_likelihood_interval,x=x,num_matrix=num_matrix,already_inspection_time=already_inspection_time,method = "BFGS")$par
    theta_est=exp(para_est[2]+para_est[3]*x)
    ### Update the parameters
    inspection_time=rep(0,num_of_inspection)
    ### Update the future inspections
    for (i in 1:3) {
      n=sum(already_inspection_time[i,]!=0)
      if(n!=num_of_inspection){
        initial_inspection=already_inspection_time[i,n]+(1:(num_of_inspection-n))*2500/((num_of_inspection-3)*i)
        inspection=c(already_inspection_time[i,1:n],optim(initial_inspection, TR_conditional, x=x[1],theta=theta_est[i],alpha=para_est[1],t_past=already_inspection_time[i,1:n],method = "BFGS")$par)
      }
      else{
        inspection=already_inspection_time[i,1:n]
      }
      inspection_time=rbind(inspection_time,inspection)
    }
    inspection_time=inspection_time[-1,]
  }
  tr=-TR(inspection_time[1,],theta1,alpha,x1)-TR(inspection_time[2,],theta2,alpha,x2)-TR(inspection_time[3,],theta3,alpha,x3)
  ### Calculate the final trace
  return(c(para_est,para_fix,tr))
}

##### Calculation for the information function
Group_information_matrix=function(t,theta,alpha,x,Ni){
  n=length(t)
  F=0
  p=c()
  ###interval
  P=c()
  for (i in 1:n) {
    p=c(p,pweibull(t[i],alpha,theta)-F)
    F=pweibull(t[i],alpha,theta)
    P=c(P,F)
  }
  p=c(p,1-F)
  Ai=alpha*(1-P)*log(1-P)
  Li=log(t/theta)
  Ndp=Ni/p
  information=matrix(0,3,3)
  coeff_b0=c(Ai[1])
  coeff_a=c(Li[1]*Ai[1])
  for (j in 2:n) {
    coeff_b0=c(coeff_b0,Ai[j]-Ai[j-1])
    coeff_a=c(coeff_a,Li[j]*Ai[j]-Li[j-1]*Ai[j-1])
  }
  coeff_b0=c(coeff_b0,-Ai[j])
  coeff_b1=x*coeff_b0
  coeff_a=(-1/alpha)*c(coeff_a,-Li[j]*Ai[j])
  coeff=rbind(coeff_b0,coeff_b1,coeff_a)
  for (i in 1:3) {
    for (j in 1:3) {
      information[i,j]=sum(coeff[i,]*coeff[j,]*Ndp)
    }
  }
  return(information)
}

##### Determinant
DET=function(t,theta,alpha,x){
  inspection=matrix(t,3,length(t)/3,byrow=T)
  information=matrix(0,3,3)
  for (i in 1:3) {
    information=information+Group_information_matrix(inspection[i,],theta[i],alpha,x[i],Ni)
  }
  return(det(information))
}


