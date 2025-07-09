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
Ni=20

##### The decomposed trace function for Weibull CSALT
TR=function(t,theta,alpha,x){
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
TR_conditional=function(t_past,t_fu,theta,alpha,x){
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
optim(c(1500,2000), TR_conditional, x=x3,theta=theta3,t_past=c(500,1000,1200),method = "BFGS")
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
  para_fix=optim(c(3.5,20,-2),Weibull_likelihood_interval,num_of_failure=fix_num_matrix,check_time=fix_inspection,method = "BFGS",hessian = T)$par
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

##### We can try a simulation like below
set.seed(999)
num_in_group=20
vol_trace_simulation(theta=theta,x=x,alpha=4,num_in_group=20,num_of_inspection=6)
##### We have the following outcome
### [1]     4.068165    21.033252    -2.871754     4.933873    19.374498    -2.523628 19287.048608

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

##### Coffin-Manson-Arrhenius Model (Multi-stress)
alpha=4
deltaT=c(100,100,120,120)
Tmax=c(373,393,373,393)
K=8.623e-05
b=1.9
Ea=0.17
A=10^4
theta_CMA=A*deltaT^(-1.9)*exp(Ea/(K*Tmax))
Ni=10
num_in_group=10

TR_CMA=function(t,theta,alpha,deltaT,Tmax){
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
  coeff=c((1+log(deltaT)^2+(1/(K*Tmax))^2+1/(alpha^2)*log(t[1]/theta)^2)*Ai[1]^2)
  for (i in 2:n) {
    coe=(1+log(deltaT)^2)*(Ai[i]-Ai[i-1])^2+1/(alpha^2)*(log(t[i]/theta)*Ai[i]-log(t[i-1]/theta)*Ai[i-1])^2
    coe1=(1+(1/(K*Tmax))^2)*(Ai[i]-Ai[i-1])^2+1/(alpha^2)*(log(t[i]/theta)*Ai[i]-log(t[i-1]/theta)*Ai[i-1])^2
    coeff=c(coeff,coe+coe1)
  }
  coeff=c(coeff,(1+log(deltaT)^2+(1/(K*Tmax))^2+1/(alpha^2)*log(t[i]/theta)^2)*Ai[i]^2)
  trace=sum(coeff*Ndp)
  return(-trace)
}

TR_CMA_conditional=function(t_past,t_fu,theta,alpha,deltaT,Tmax){
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
  coeff=c((1+log(deltaT)^2+(1/(K*Tmax))^2+1/(alpha^2)*log(t[1]/theta)^2)*Ai[1]^2)
  for (i in 2:n) {
    coe=(1+log(deltaT)^2)*(Ai[i]-Ai[i-1])^2+1/(alpha^2)*(log(t[i]/theta)*Ai[i]-log(t[i-1]/theta)*Ai[i-1])^2
    coe1=(1+(1/(K*Tmax))^2)*(Ai[i]-Ai[i-1])^2+1/(alpha^2)*(log(t[i]/theta)*Ai[i]-log(t[i-1]/theta)*Ai[i-1])^2
    coeff=c(coeff,coe+coe1)
  }
  coeff=c(coeff,(1+log(deltaT)^2+(1/(K*Tmax))^2+1/(alpha^2)*log(t[i]/theta)^2)*Ai[i]^2)
  trace=sum(coeff*Ndp)
  return(-trace)
}

###Varify the effectiveness of TR_CMA
optim(c(300,400,500), TR_CMA_conditional,t_past=c(100,200),theta=theta_CMA[1],alpha=4,deltaT=deltaT[1],Tmax=Tmax[1],method = "BFGS")
###$par [1] 306.5916 369.3287 429.1506
log_deltaT=log(deltaT)
inverse_KTmax=1/(K*Tmax)


### Likelihood function
Weibull_likelihood_interval_CMA=function(parameter,num_of_failure,check_time){
  alpha=parameter[1]
  A=exp(parameter[2])
  b=parameter[3]
  Ea=parameter[4]
  theta=A*deltaT^(-b)*exp(Ea/(K*Tmax))
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

check_time=matrix(c(100,150,300,100,150,300,100,150,300,100,150,300),4,3,byrow=T)
num_of_failure=matrix(c(0,1,10,9,1,2,15,2,1,3,15,1,2,7,11,0),4,4,byrow=T)
optim(c(3.5,9,1.9,0.17),Weibull_likelihood_interval_CMA,num_of_failure=num_of_failure,check_time=check_time,method = "BFGS")

### Optimize the unconstrained trace
optim(c(100,150,200,300,400,500), TR_CMA,theta=theta_CMA[1],alpha=4,deltaT=deltaT[1],Tmax=Tmax[1],method = "BFGS")
### [1] 252.9793 307.0770 347.9708 384.7369 421.9655 465.8053
### [1] -150036.9

optim(c(120,150,200,300,400,450), TR_CMA,theta=theta_CMA[2],alpha=4,deltaT=deltaT[2],Tmax=Tmax[2],method = "BFGS")
### [1] 193.3036 234.6562 265.9095 294.0111 322.4677 355.9813
### [1] -137330.4

optim(c(120,150,200,300,400,450), TR_CMA,theta=theta_CMA[3],alpha=4,deltaT=deltaT[3],Tmax=Tmax[3],method = "BFGS")
### [1] 178.9093 217.1698 246.0910 272.0933 298.4242 329.4292
### [1] -152329.7

optim(c(120,150,200,250,300,310), TR_CMA,theta=theta_CMA[4],alpha=4,deltaT=deltaT[4],Tmax=Tmax[4],method = "BFGS")
### [1] 136.7086 165.9571 188.0573 207.9310 228.0569 251.7520
### [1] -137593.6
###### Sum=577290.6

### Fix-time design
 
 TR_CMA(c(100,150,200,250,300,350),theta_CMA[1],4,deltaT[1],Tmax[1])
 ### [1] -122960.9
 TR_CMA(c(100,150,200,250,300,350),theta_CMA[2],4,deltaT[2],Tmax[2])
 ### [1] -132262
 TR_CMA(c(100,150,200,250,300,350),theta_CMA[3],4,deltaT[3],Tmax[3])
 ### [1] -145664.4
 TR_CMA(c(100,150,200,250,300,350),theta_CMA[4],4,deltaT[4],Tmax[4])
 ### [1] -124633.9
###### Sum=525521.2
###### Relative efficiency=525521.2/577290.6=0.9103235

CMA_Failure_times=function(theta,alpha,num_in_group){
  failure_time=matrix(0,4,num_in_group)
  for (i in 1:4) {
    failure_time[i,]=sort(rweibull(num_in_group,alpha,theta[i]))
  }
  return(failure_time)
}

num_in_inspection_CMA=function(failure_time,inspection_time){
  n=dim(inspection_time)[2]
  num_matrix=matrix(0,4,n+1)
  for (i in 1:4) {
    num_matrix[i,1]=sum(failure_time[i,]<=inspection_time[i,1])
    for (j in 2:n) {
      num_matrix[i,j]=sum((failure_time[i,]>inspection_time[i,j-1])&(failure_time[i,]<=inspection_time[i,j]))
    }
    num_matrix[i,j+1]=sum(failure_time[i,]>inspection_time[i,j])
  }
  return(num_matrix)
}

CMA_trace_simulation=function(theta,alpha,deltaT,Tmax,num_in_group,num_of_inspection){
  failure_time=CMA_Failure_times(theta,alpha,num_in_group)
  inspection_initial=matrix(c(100,150,100,150,100,150,100,150),4,2,byrow=T)
  num_matrix_initial=num_in_inspection_CMA(failure_time,inspection_initial)
  fix_inspection=matrix(50+rep(50*(1:num_of_inspection),4),4,num_of_inspection,byrow = T)
  fix_num_matrix=num_in_inspection_CMA(failure_time,fix_inspection)
  para_fix=optim(c(3.5,9,1.9,0.17),Weibull_likelihood_interval_CMA,num_of_failure=fix_num_matrix,check_time=fix_inspection,method = "BFGS",hessian = T)
  lowerlimit1=para_fix$par-1.96*sqrt(diag(solve(para_fix$hessian)))
  upperlimit1=para_fix$par+1.96*sqrt(diag(solve(para_fix$hessian)))
  para_fix=para_fix$par
  para_est=optim(c(3.5,9,1.9,0.17),Weibull_likelihood_interval_CMA,num_of_failure=num_matrix_initial,check_time=inspection_initial,method = "BFGS")$par
  theta_est=exp(para_est[2])*deltaT^(-para_est[3])*exp(para_est[4]/(K*Tmax))
  inspection_time=rep(0,num_of_inspection)
  for (i in 1:4) {
    inspection=c(100,150,optim(c(200,250,300,350), TR_CMA_conditional,t_past=c(100,150),theta=theta_est[i],alpha=4,deltaT=deltaT[i],Tmax=Tmax[i],method = "BFGS")$par)
    inspection_time=rbind(inspection_time,inspection)
  }
  inspection_time=inspection_time[-1,]
  
  already_inspection_time=matrix(0,4,num_of_inspection)
  already_inspection_time[,c(1,2)]=matrix(c(100,150,100,150,100,150,100,150),4,2,byrow=T)
  while (already_inspection_time[1,num_of_inspection]==0) {
    position=nonzerominposition(inspection_time-already_inspection_time)
    ### Next inspection
    if(inspection_time[position]>max(already_inspection_time)){
      already_inspection_time[position]=inspection_time[position]
      ### For the spectial case that after a inspection, the next suggested time is even more earlier
    }else{
      already_inspection_time[position]=max(already_inspection_time)
    }
    num_matrix=num_in_inspection_CMA(failure_time,inspection_time)
    para_est=optim(c(3.5,9,1.9,0.17),renew_Weibull_likelihood_interval_CMA,num_matrix=num_matrix,already_inspection_time=already_inspection_time,method = "BFGS")$par
    theta_est=exp(para_est[2])*deltaT^(-para_est[3])*exp(para_est[4]/(K*Tmax))
    
    inspection_time=rep(0,num_of_inspection)
    for (i in 1:4) {
      n=sum(already_inspection_time[i,]!=0)
      if(n!=num_of_inspection){
        initial_inspection=already_inspection_time[i,n]+(1:(num_of_inspection-n))*50
        inspection=c(already_inspection_time[i,1:n],optim(initial_inspection, TR_CMA_conditional,t_past=already_inspection_time[i,1:n],theta=theta_CMA[i],alpha=para_est[1],deltaT=deltaT[i],Tmax=Tmax[i],method = "BFGS")$par)
      }
      else{
        inspection=already_inspection_time[i,1:n]
      }
      inspection_time=rbind(inspection_time,inspection)
    }
    inspection_time=inspection_time[-1,]
  }
  final_para=optim(c(4.5,8.8,1.77,0.17),renew_Weibull_likelihood_interval_CMA,num_matrix=num_matrix,already_inspection_time=already_inspection_time,method = "BFGS",hessian = T)
  lowerlimit=final_para$par-1.96*sqrt(diag(solve(final_para$hessian)))
  upperlimit=final_para$par+1.96*sqrt(diag(solve(final_para$hessian)))
  tr=-TR_CMA(inspection_time[1,],theta_CMA[1],4,deltaT[1],Tmax[1])-TR_CMA(inspection_time[2,],theta_CMA[2],4,deltaT[2],Tmax[2])-TR_CMA(inspection_time[3,],theta_CMA[3],4,deltaT[3],Tmax[3])-TR_CMA(inspection_time[4,],theta_CMA[4],4,deltaT[4],Tmax[4])
  return(list(inspection_time,num_matrix,fix_num_matrix,c(para_est,para_fix,tr),round(c(lowerlimit,upperlimit),4),round(c(lowerlimit1,upperlimit1),4)))
}

renew_Weibull_likelihood_interval_CMA=function(parameter,num_matrix,already_inspection_time){
  alpha=parameter[1]
  A=exp(parameter[2])
  b=parameter[3]
  Ea=parameter[4]
  theta=A*deltaT^(-b)*exp(Ea/(K*Tmax))
  m=dim(num_matrix)[1]
  logL=0
  for (i in 1:m) {
    p=c()
    F=0
    n=sum(already_inspection_time[i,]!=0)
    for (j in 1:n) {
      p=c(p,pweibull(already_inspection_time[i,j],alpha,theta[i])-F)
      F=pweibull(already_inspection_time[i,j],alpha,theta[i])
    }
    p=c(p,1-F)
    num_of_censoring=num_in_group-sum(num_matrix[i,1:n])
    logL=logL+sum(c(num_matrix[i,1:n],num_of_censoring)*log(p))
  }
  return(-logL)
}

##### Illustrative Example
#2499(20)
set.seed(2499)
num_in_group=20; Ni=20
CMA_trace_simulation(theta_CMA,4,deltaT,Tmax,num_in_group=20,6)

#inspection  100  150 281.0735 339.9512 392.4031 439.0734
#inspection  100  150 223.2738 267.9430 306.0326 344.8885
#inspection  100  150 214.5208 252.5675 284.5697 319.9307
#inspection  100  150 181.6986 214.5208 235.1073 260.6031
### num matrix of RITO
#[,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    0    0    4    3   11    1    1
#[2,]    1    1    5    8    4    0    1
#[3,]    2    3    9    3    1    2    0
#[4,]    3    6    5    5    1    0    0
### num matrix of fixed-time
#[,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    0    0    0    3    3    6    8
#[2,]    1    1    5    5    7    0    1
#[3,]    2    3    8    4    2    1    0
#[4,]    3    6    8    3    0    0    0

### RITO CI             [3.3889,5.0457] [7.0394,14.3649] [1.8924,3.0902] [0.1453,0.2822]
### Fixed-time Estimate [3.0679,4.6246] [6.8428,15.0928] [1.8460,3.2474] [0.1341,0.2925]
###Relative Efficiency=1.131133e+06/2/577290=0.9796922
