setwd("/Users/lizhuo/Downloads")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(MASS)
library(mvtnorm)
library(MCMCpack)
library(mnormt)
library(stringr)
##############################################################################################
### Bayesian analysis of H1b




# find unique SOC names -- job category
job <- sort(table(df16$SOC_NAME))
# top jobs
job_title <- c("SOFTWARE DEVELOPERS, APPLICATIONS", "COMPUTER SYSTEMS ANALYSTS",
               "COMPUTER OCCUPATIONS, ALL OTHER", "SOFTWARE DEVELOPERS, SYSTEMS SOFTWARE",
               "COMPUTER SYSTEMS ANALYST", "MANAGEMENT ANALYSTS", "ACCOUNTANTS AND AUDITORS",
               "NETWORK AND COMPUTER SYSTEMS ADMINISTRATORS", "MECHANICAL ENGINEERS", 
               "FINANCIAL ANALYSTS", "DATABASE ADMINISTRATORS", "MARKET RESEARCH ANALYSTS AND MARKETING SPECIALISTS",
               "OPERATIONS RESEARCH ANALYST", "ELECTRONICS ENGINEERS, EXCEPT COMPUTER", 
               "COMPUTER AND INFORMATION SYSTEMS MANAGERS", "ELECTRICAL ENGINEERS", "PHYSICIANS AND SURGEONS, ALL OTHER",
               "MEDICAL SCIENTISTS, EXCEPT EPIDEMIOLOGISTS", "WEB DEVELOPERS", "STATISTICIANS")
# recategory job title
df16$Job <- df16$SOC_NAME

df16$Job[is.na(df16$Job)] <- df16$SOC_NAME[is.na(df16$Job)]

# find number of certified and uncertified
table(df16$CASE_STATUS)
# certified 569646
# certified-withdrawn 47092
# denied 9175
# withdrawn 21890 
# total 647803

table(df15$CASE_STATUS)
# certified 547278
# certified-withdrawn 41071
# denied 10923
# withdrawn 19455
# total 618727

table(df14$CASE_STATUS)
# 455144
# 36350
# 11896
# 16034
# 519427

table(df13$CASE_STATUS)
# 382951
# 35432
# 12126
# 11590
# 442114

table(df12$CASE_STATUS)
# 352668
# 31118
# 21096
# 10725
# 415607

table(df11$CASE_STATUS)
# 307936
# 11596
# 29130
# 10105
# 358767



##############################################################################################
# create summary dataset for yearly certificate data
df_cert <- matrix(data=c(569646, 47092, 9175, 21890, 647803, 
                         547278, 41071, 10923, 19455, 618727,
                         455144, 36350, 11896, 16034, 519427,
                         382951, 35432, 12126, 11590, 442114,
                         352668, 31118, 21096, 10725, 415607,
                         307936, 11596, 29130, 10105, 358767), byrow=T, nrow=6)
colnames(df_cert) <- c("certified", "certwit", "denied", "withdrawn", "total")
df_cert <- as.data.frame(df_cert)
df_cert$P_c <- (df_cert$certified)/df_cert$total
df_cert$P_d <- df_cert$denied/df_cert$total
df_cert$P_w <- (df_cert$certwit+df_cert$withdrawn)/df_cert$total




##############################################################################################
# For FY 2016 by state
df16$State[df16$State=="NA"] <- NA
df16 <- df16 %>% drop_na()
df16 <- na.omit(df16)
df16 <- df16[!is.na(df16$State), ]

table <- table(df16$State, df16$CASE_STATUS)
table <- as.data.frame(table)
t_cert <- subset(table, table$Var2=="CERTIFIED")
t_cewi <- subset(table, table$Var2=="CERTIFIED-WITHDRAWN")
t_deni <- subset(table, table$Var2=="DENIED")
t_widr <- subset(table, table$Var2=="WITHDRAWN")
sum(t_cert[,1]==t_cewi[,1])
sum(t_cewi[,1]==t_deni[,1])
sum(t_deni[,1]==t_widr[,1])
FY16 <- cbind(t_cert, t_cewi[,3], t_deni[,3], t_widr[,3])
FY16$Var2 <- NULL

write_csv(FY16, "FY16.csv")

temp <- read_csv("FY16.csv")

# name the columns
colnames(FY16) <- c("State", "Cert", "Cewi", "Deni", "Widr")
FY16 <- FY16[-28,]
rownames(FY16) <- FY16[,1]
FY16[,1] <- NULL

df2 <- NA
row_sum <- apply(FY16, 1, sum)
df2$CE <- FY16$Cert/row_sum
df2$DE <- FY16$Deni/row_sum
df2$WD <- (FY16$Cewi+FY16$Widr)/row_sum
df2$Pro <- row_sum/sum(FY16)
df2 <- as.data.frame(df2)
df2$NA. <- NULL

df <- df2
df[,4]<- sum(FY16)*df[,4]
df[,1]<- df[,4]*df2[,1]
df[,2]<- df[,4]*df2[,2]
df[,3]<- df[,4]*df2[,3]
df <- df[,1:3]
df2 <- as.data.frame(df2)

n=100000

a1 <- as.data.frame(matrix(NA, nrow=n, ncol=52))
a2 <- as.data.frame(matrix(NA, nrow=n, ncol=52))
a3 <- as.data.frame(matrix(NA, nrow=n, ncol=52))
for (i in 1:52){
  a1[,i] <- rbeta(n, as.numeric(df[i,][1])+1, as.numeric(df[i,][2])+as.numeric(df[i,][3])+2)
  a2[,i] <- rbeta(n, as.numeric(df[i,][2]+1), as.numeric(df[i,][1])+as.numeric(df[i,][3])+2)
  a3[,i] <- rbeta(n, as.numeric(df[i,][3]+1), as.numeric(df[i,][2])+as.numeric(df[i,][1])+2)
}
par(mfrow=c(1,1))
hist(apply(a1-a2, 1, sum), breaks=14, main="simple model")
##############################################################################################
### Hierarchical model
y1 <- FY16$Cert                # number of applications that are certified
y2 <- FY16$Deni                # number of applications denied
y3 <- FY16$Cewi+FY16$Widr      # number of applications withdrawn

# posterior
f <- function(y1, y2, y3, b1, b2, mu1, mu2, t1, t2, rho){
  theta1<-theta2<-theta3<-out <- rep(0, 16)
  for (j in 1:52){
    theta1[j] <- (exp(b1[j])*exp(b2[j]))/((1+exp(b1[j]))*(1+exp(b2[j])))
    theta2[j] <- (exp(b2[j])/((1+exp(b1[j]))*(1+exp(b2[j]))))
    theta3[j] <- (1/(1+exp(b2[j])))
    s <- matrix(c(t1^2, rho*t1*t2, rho*t1*t2, t2^2),ncol=2,nrow=2,byrow=T)
    out[j] <- ddirichlet(c(theta1[j],theta2[j],theta3[j]), c(y1[j],y2[j],y3[j]))*(dmvnorm(c(b1[j],b2[j]),mean=c(mu1, mu2),s))
  }
  return(prod(out))
}

## Simulation
N=200000
# starting point
theta1=apply(a1, 2, median)
theta2=apply(a2, 2, median)
theta3=apply(a3, 2, median)
alpha1 <- theta1/(theta1+theta2)
alpha2 <- 1-theta3
beta1 <- as.data.frame(log(1/(1-alpha1)))
beta2 <- as.data.frame(log(1/(1-alpha2)))
tao <- c(sqrt(var(beta1)),sqrt(var(beta2)))
rho <- cor(beta1, beta2)
mu <- cbind(apply(beta1, 2, mean), apply(beta2, 2, mean))
beta1 <- rbind(t(beta1), matrix(NA, nrow=N, ncol=52))
beta2 <- rbind(t(beta2), matrix(NA, nrow=N, ncol=52))
tao <- rbind(tao, matrix(NA, nrow=N, ncol=2))
rho <- c(0, numeric(N))
mu <- rbind(mu, matrix(NA, nrow=N, ncol=2))
accp <- numeric(N)



# MCMC
for ( i in 1:(N)){
  # proposal
  p_beta <- matrix(nrow=52,ncol=2)
  s <- matrix(c(0.0001,0,0,0.0001), nrow=2, ncol=2, byrow=T)
  for (j in 1:52){
    mean <- c(beta1[i,j], beta2[i,j])
    p_beta[j,] <- mvrnorm(1, mu=mean,s)
  }
  p_mu <- rmvnorm(1, mu[i,], s)
  p_tao <- c(tao[1,1]+runif(1, -0.0025, 0.0025), tao[1,2]+runif(1, -0.0025, 0.0025))
  p_rho <- rho[i]+runif(1, -0.0025, 0.0025)
  
  if (p_tao[1]<=0 || p_tao[2]<=0 || abs(p_rho)>=1){
    break
  }
  else{
    r <- f(y1, y2, y3, t(p_beta[,1]), t(p_beta[,2]), p_mu[1], p_mu[2], p_tao[1], p_tao[2], p_rho)/f(y1, y2, y3, beta1[i,], beta2[i,], mu[i,1], mu[i,2], tao[i,1], tao[i,2], rho[i])
    accp[i] <- min(1,r)
    accept <- rbinom(1,1, min(1,r))
    if (accept){
      beta1[i+1,]=t(p_beta[,1])
      beta2[i+1,]=t(p_beta[,2])
      mu[i+1,]=p_mu
      tao[i+1,]=p_tao
      rho[i+1]=rho[i]
    }
    else{
      beta1[i+1,]=beta1[i,]
      beta2[i+1,]=beta2[i,]
      mu[i+1,]=mu[i,]
      tao[i+1,]=tao[i,]
      rho[i+1]=rho[i]
    }
  }
  
  i=+1
  
}

# average acceptance rate of Metropolis Hasting
mean(accp)
# find alpha
alpha1 <- exp(beta1)/(1+exp(beta1))
alpha2 <- exp(beta2)/(1+exp(beta2))
#for (i in 1:52){
#  print (quantile(alpha1[,i],probs=c(.025,.25,.5,.75,.975), na.rm=T))
#}
par(mfrow=c(1,1))
plot(beta1[,1], type="l")
plot(alpha1[,1], type="l")

par(mfrow=c(4, 4))
for (i in 1:52){
  plot(beta1[,i], type="l")
}

par(mfrow=c(4, 4))
for (i in 1:52){
  plot(beta2[,i], type="l")
}


par(mfrow=c(4, 4))
for (i in 1:52){
  plot(alpha1[,i], type="l")
}

par(mfrow=c(4, 4))
for (i in 1:52){
  plot(alpha2[,i], type="l")
}

colnames(alpha1) <- names
alpha1 <- as.data.frame(alpha1)
par(mfrow=c(1,1))
plot(alpha1$VIRGINIA, type="l", main="Trace plot of MCMC Simulation", ylab="VIRGINIA", xlab="ITERATION")
plot(alpha1$CALIFORNIA, type="l", main="California")
plot(alpha1$`DISTRICT OF COLUMBIA`, type="l", main="DC")
plot(alpha1$`NEW YORK`, type="l", main="New York")

#### Save the 10k iteration results.....1 hour to run
# convergence not great, like half? has pretty good convergence...but thank god Virginia converges///
ite10kbeta1 <- as.data.frame(beta1)
ite10kbeta2 <- as.data.frame(beta2)
ite10kalpha1 <- as.data.frame(alpha1)
ite10kalpha2 <- as.data.frame(alpha2)
ite10kmu <- as.data.frame(mu)
ite10ktao <- as.data.frame(tao)
ite10krho <- as.data.frame(rho)

write_csv(ite10kbeta1, "10kb1.csv")
write_csv(ite10kbeta2, "10kb2.csv")
write_csv(ite10kalpha1, "10ka1.csv")
write_csv(ite10kalpha2, "10ka2.csv")
write_csv(ite10kmu, "10kmu.csv")
write_csv(ite10ktao, "10ktao.csv")
write_csv(ite10krho, "10kr.csv")


#### Save the 20k iteration results.....1 hour to run
# convergence not great, like half? has pretty good convergence...but thank god Virginia converges///
ite20kbeta1 <- as.data.frame(beta1)
ite20kbeta2 <- as.data.frame(beta2)
ite20kalpha1 <- as.data.frame(alpha1)
ite20kalpha2 <- as.data.frame(alpha2)
ite20kmu <- as.data.frame(mu)
ite20ktao <- as.data.frame(tao)
ite20krho <- as.data.frame(rho)

write_csv(ite20kbeta1, "10kb1.csv")
write_csv(ite20kbeta2, "10kb2.csv")
write_csv(ite20kalpha1, "10ka1.csv")
write_csv(ite20kalpha2, "10ka2.csv")
write_csv(ite20kmu, "10kmu.csv")
write_csv(ite20ktao, "10ktao.csv")
write_csv(ite20krho, "10kr.csv")


#################################################################
### graphics and other steps MODEL CHECKING!!!!!!! DONT FORGET!
theta1=alpha1*alpha2
theta3=1-alpha2
theta2=1-theta1-theta3
names <- rownames(FY16)
colnames(theta1) <- names
colnames(theta2) <- names
colnames(theta3) <- names

## take out first half as burn-in
# with N=20k
t1 <- theta1[((N+1)/2):(N+1),]
t2 <- theta2[((N+1)/2):(N+1),]
t3 <- theta3[((N+1)/2):(N+1),]

for (i in 1:52){
  theta1[,i] <- theta1[,i]*as.numeric(df2[i, 4])
  theta2[,i] <- theta2[,i]*as.numeric(df2[i,4])
}

par(mfrow=c(1,1))
hist(apply(theta1[50001:N+1,]-theta2[50001:N+1,], 1, sum), xlab="Simulation of Hierachical Model", main="Histgram of Hierarhical Model", breaks=17)


theta1=alpha1*alpha2
theta3=1-alpha2
theta2=1-theta1-theta3
names <- rownames(FY16)e
colnames(theta1) <- names
theta1 <- as.data.frame(theta1)

hist(theta1$VIRGINIA, main="Virginia")
abline(v=df2$CE[48], col="orange", lwd=5, lty=4)

hist(beta1[,48], main="Virginia")

## simulate y
ind <- sample(100000:N, 5000)
t1 <- theta1[ind,]
t2 <- theta2[ind,]
t3 <- theta3[ind,]


yrep <- matrix(NA, nrow=2000, ncol=3)
for (i in 1:2000){
  yrep[i,] <- rmultinom(1, size=18901, prob=c(t1[i,48], t2[i,48], t3[i,48]))
}
hist(yrep[,1], main="Histogram of Replicated Data for Virginia", xlab="Virginia", breaks=20)
abline(v=16670, col="brown4", lwd=3, lty=4)

for (i in 1:2000){
  yrep[i,] <- rmultinom(1, size=sum(FY16[5,]), prob=c(t1$CALIFORNIA[i], t2[i,5], t3[i,5]))
}
hist(yrep[,1], main="Histogram of Replicated Data for California", xlab="Cal", breaks=20)
abline(v=FY16[5,], col="brown4", lwd=3, lty=4)
sum(yrep[,1]>=FY16[5,1])/2000

for (i in 1:2000){
  yrep[i,] <- rmultinom(1, size=sum(FY16[9,]), prob=c(t1[i,9], t2[i,9], t3[i,9]))
}
hist(yrep[,1], main="Histogram of Replicated Data for DC", xlab="DC", breaks=20)
abline(v=FY16[9,], col="brown4", lwd=3, lty=4)
sum(yrep[,1]>=FY16[9,1])/2000

for (i in 1:2000){
  yrep[i,] <- rmultinom(1, size=sum(FY16[33,]), prob=c(t1[i,33], t2[i,33], t3[i,33]))
}
hist(yrep[,1], main="Histogram of Replicated Data for New York", xlab="New York", breaks=20)
abline(v=FY16[33,1], col="brown4", lwd=3, lty=4)
sum(yrep[,1]>=FY16[33,1])/2000


quantile(theta1$VIRGINIA[100000:200000],probs=c(.025,.25,.5,.75,.975), na.rm=T)
quantile(theta1$CALIFORNIA[100000:200000],probs=c(.025,.25,.5,.75,.975), na.rm=T)
quantile(theta1$`DISTRICT OF COLUMBIA`[100000:200000],probs=c(.025,.25,.5,.75,.975), na.rm=T)
quantile(theta1$`NEW YORK`[100000:200000],probs=c(.025,.25,.5,.75,.975), na.rm=T)
#2.5%       25%       50%       75%     97.5% 
#0.8774073 0.8804264 0.8819733 0.8834897 0.8864227 
t1 <- as.data.frame(t1)


p1 <- density(t1$VIRGINIA)
p2 <- density(t1$`DISTRICT OF COLUMBIA`)
p3 <- density(t1$CALIFORNIA)
p4 <- density(t1$`NEW YORK`)
p5 <- density(t1$MASSACHUSETTS)
plot(p1, col=rgb(0.647059, 0.164706, 0.164706,0.8), xlim=c(0.84,0.9), ylim=c(0,400), main="Simulated Certified Rate of Different States")
polygon(p1, col=rgb(0.647059, 0.164706, 0.164706,0.8), border="black")
lines(p5, col=c( 0.960784, 0.870588, 0.701961, 1/4))
polygon(p5, col=c( 0.960784, 0.870588, 0.701961, 1/4), border="black")
lines(p2, col=rgb(0.545098,0,0, 1/4))
polygon(p2, col=rgb(0.545098,0,0, 1/4), border="black")
lines(p3, col=rgb(0.545098, 0.270588, 0.0745098,1/4))
polygon(p3, col=rgb(0.545098, 0.270588, 0.0745098,1/4), border="black")
lines(p4, col=rgb(0.823529, 0.705882, 0.54902, 1/4))
polygon(p4, col=rgb(0.823529, 0.705882, 0.54902, 1/4), border="black")

legend("topright", legend=c("VIRGINIA", "MASSACHUSETTS","DC", "CALIFORNIA", "NY"), 
       fill=c(rgb(0.647059, 0.164706, 0.164706,0.8), rgb( 0.960784, 0.870588, 0.701961, 1/4), rgb(0.545098,0,0, 1/4), rgb(0.545098, 0.270588, 0.0745098,1/4), rgb(0.823529, 0.705882, 0.54902, 1/4)))






##############################################################################################
### then analysing the job.....

# select the top 24 job titles with the most applicants, recategory the other jobs as "other"
job <- sort(table(df16$SOC_NAME))
job <- as.data.frame(job)
job <- job[job$Freq!=0,]

## top SOC job category
top <- as.character(tail(job$Var1, 24))


df16$SOC_NAME <- as.character(df16$SOC_NAME)
df16$Job <- ifelse(df16$SOC_NAME %in% top, df16$SOC_NAME, "Other")

# create contingency table
table_job <- table(df16$Job, df16$CASE_STATUS)



##############################################################################################
### actually just focus on virginia and DC jobs
VA <- subset(df16, df16$State=="VIRGINIA")
DC <- subset(df16, df16$State=="DISTRICT OF COLUMBIA")

VD <- rbind(VA)

job <- sort(table(VD$SOC_NAME))
job <- as.data.frame(job)
job <- job[job$Freq!=0,]
top <- as.character(tail(job$Var1, 24))
top


VD$SOC_NAME <- as.character(VD$SOC_NAME)
VD$Job <- ifelse(VD$SOC_NAME %in% top, VD$SOC_NAME, "Other")

# create contingency table
table_job <- as.data.frame(table(VD$Job, VD$CASE_STATUS))

t_certv <- subset(table_job, table_job$Var2=="CERTIFIED")
t_cewiv <- subset(table_job, table_job$Var2=="CERTIFIED-WITHDRAWN")
t_deniv <- subset(table_job, table_job$Var2=="DENIED")
t_widrv <- subset(table_job, table_job$Var2=="WITHDRAWN")
sum(t_certv[,1]==t_cewiv[,1])
sum(t_cewiv[,1]==t_deniv[,1])
sum(t_deniv[,1]==t_widrv[,1])
VD16 <- cbind(t_certv, t_cewiv[,3], t_deniv[,3], t_widrv[,3])
VD16$Var2 <- NULL

colnames(VD16) <- c("Job", "Cert", "Cewi", "Deni", "Widr")
rownames(VD16) <- VD16[,1]
VD16[,1] <- NULL


##############################################################################################
### use previous year data as prioir
VA15 <- subset(df15, df15$State=="VIRGINIA")
DC15 <- subset(df15, df15$State=="DISTRICT OF COLUMBIA")

VD15 <- rbind(VA15)
VD15$SOC_NAME <- as.character(VD15$SOC_NAME)
VD15$Job <- ifelse(VD15$SOC_NAME %in% top, VD15$SOC_NAME, "Other")

# create contingency table
table_job15 <- as.data.frame(table(VD15$Job, VD15$CASE_STATUS))

t_certv15 <- subset(table_job15, table_job15$Var2=="CERTIFIED")
t_cewiv15 <- subset(table_job15, table_job15$Var2=="CERTIFIED-WITHDRAWN")
t_deniv15 <- subset(table_job15, table_job15$Var2=="DENIED")
t_widrv15 <- subset(table_job15, table_job15$Var2=="WITHDRAWN")
sum(t_certv15[,1]==t_cewiv15[,1])
sum(t_cewiv15[,1]==t_deniv15[,1])
sum(t_deniv15[,1]==t_widrv15[,1])
VD15 <- NULL
VD15 <- cbind(t_certv15, t_cewiv15[,3], t_deniv15[,3], t_widrv15[,3])
VD15$Var2 <- NULL

colnames(VD15) <- c("Job", "Cert", "Cewi", "Deni", "Widr")
rownames(VD15) <- VD15[,1]
VD15[,1] <- NULL



## year 2016
df216 <- NULL
row_sum <- apply(VD16, 1, sum)
df216$CE <- VD16$Cert/row_sum
df216$DE <- VD16$Deni/row_sum
df216$WD <- (VD16$Cewi+VD16$Widr)/row_sum
df216$Pro <- row_sum/sum(VD16)
df216 <- as.data.frame(df216)

df16 <- df216
df16[,4]<- sum(VD16)*df16[,4]
df16[,1]<- df16[,4]*df216[,1]
df16[,2]<- df16[,4]*df216[,2]
df16[,3]<- df16[,4]*df216[,3]
df16 <- df16[,1:3]
df216 <- as.data.frame(df216)

## year 2015
df215 <- NULL
row_sum <- apply(VD15, 1, sum)
df215$CE <- VD15$Cert/row_sum
df215$DE <- VD15$Deni/row_sum
df215$WD <- (VD15$Cewi+VD16$Widr)/row_sum
df215$Pro <- row_sum/sum(VD15)
df215 <- as.data.frame(df215)

df15 <- df215
df15[,4]<- sum(VD15)*df15[,4]
df15[,1]<- df15[,4]*df215[,1]
df15[,2]<- df15[,4]*df215[,2]
df15[,3]<- df15[,4]*df215[,3]
df15 <- df15[,1:3]
df215 <- as.data.frame(df215)


##############################################################################################
### Hierarchical model for jobs, use Dir as prior
# alpha
alpha <- df16+df15*0.1


ite=5000
t1 <- matrix(NA, nrow=ite, ncol=25)
names <- rownames(df16)
t2 <- matrix(NA, nrow=ite, ncol=25)
t3 <- matrix(NA, nrow=ite, ncol=25)
colnames(t1) <- colnames(t2) <- colnames(t3) <- names

for (i in 1:ite){
  for (j in 1:25){
    t1[i,j] <- rbeta(1, alpha[j,1], sum(alpha[j,])-alpha[j,1])
    t2[i,j] <- rbeta(1, alpha[j,2], sum(alpha[j,])-alpha[j,2])
    t3[i,j] <- rbeta(1, alpha[j,3], sum(alpha[j,])-alpha[j,3])
  }
}


t1 <- as.data.frame(t1)
t2 <- as.data.frame(t2)
t3 <- as.data.frame(t3)

par(oma=c(0,0,0,0))
p1 <- density(t1$STATISTICIANS)
p2 <- density(t1$`DATABASE ADMINISTRATORS`)
p3 <- density(t1$`FINANCIAL ANALYSTS`)
p4 <- density(t1$`MANAGEMENT ANALYSTS`)
p5 <- density(t1$`COMPUTER PROGRAMMERS`)
plot(p1, col=rgb(0.647059, 0.164706, 0.164706,0.8), xlim=c(0.60,0.98), ylim=c(0,90), main="Simulated Certified Rate of Different Jobs")
polygon(p1, col=rgb(0.647059, 0.164706, 0.164706,0.8), border="black")
lines(p5, col=rgb( 0.960784, 0.870588, 0.701961, 1/4))
polygon(p5, col=rgb( 0.960784, 0.870588, 0.701961, 1/4), border="black")
lines(p2, col=rgb(0.545098,0,0, 1/4))
polygon(p2, col=rgb(0.545098,0,0, 1/4), border="black")
lines(p3, col=rgb(0.545098, 0.270588, 0.0745098,1/4))
polygon(p3, col=rgb(0.545098, 0.270588, 0.0745098,1/4), border="black")
lines(p4, col=rgb(0.823529, 0.705882, 0.54902, 1/4))
polygon(p4, col=rgb(0.823529, 0.705882, 0.54902, 1/4), border="black")

legend("topleft", legend=c("Statsiticians", "Programmer","Database adm", "Financial Analyst", "Management Analyst"), 
       fill=c(rgb(0.647059, 0.164706, 0.164706,0.8), rgb( 0.960784, 0.870588, 0.701961, 1/4), rgb(0.545098,0,0, 1/4), rgb(0.545098, 0.270588, 0.0745098,1/4), rgb(0.823529, 0.705882, 0.54902, 1/4)))




dfh16 <- subset(dfh, dfh$YEAR=="2016")
dfh15 <- subset(dfh, dfh$YEAR=="2015")


VD$Job <- as.factor(VD$Job)
VD$PREVAILING_WAGE <- as.numeric(VD$PREVAILING_WAGE)


par(oma=c(0,15,0,0))
boxplot(PREVAILING_WAGE~Job, data=VD, ylim=c(0, 250000), horizontal=T, las=1, cex.axis = 0.7,
        col=ifelse(levels(VD$Job)=="STATISTICIANS" , rgb(0.647059, 0.164706, 0.164706,0.8), rgb(0.545098, 0.270588, 0.0745098,1/4)))


quantile(t1$STATISTICIANS,probs=c(.025,.25,.5,.75,.975), na.rm=T)
quantile(t1$`FINANCIAL ANALYSTS`,probs=c(.025,.25,.5,.75,.975), na.rm=T)
quantile(t1$`MANAGEMENT ANALYSTS`,probs=c(.025,.25,.5,.75,.975), na.rm=T)
quantile(t1$`COMPUTER PROGRAMMERS`,probs=c(.025,.25,.5,.75,.975), na.rm=T)



par(oma=c(0,0,0,0))
hist(t1$STATISTICIANS, main="Statistician")
abline(v=df216$CE[24], col="orange", lwd=5, lty=4)


hist(t1$`COMPUTER PROGRAMMERS`, main="Computer Programmer")
abline(v=df216$CE[5], col="orange", lwd=5, lty=4)



#################################################################
alpha <- df16+df15


ite=100
t1 <- rbind(df216[,1], matrix(NA, nrow=ite, ncol=25))
names <- rownames(df16)
t2 <- rbind(df216[,2], matrix(NA, nrow=ite, ncol=25))
t3 <- rbind(df216[,3], matrix(NA, nrow=ite, ncol=25))
colnames(t1) <- colnames(t2) <- colnames(t3) <- names
rho <- rbind(rep(0.8, 25), matrix(NA, nrow=ite, ncol=25))

for (i in 2:(ite+1)){
  for (j in 1:25){
    t1[i,j] <- rbeta(1, df16[j,1]+df15[j,1]*rho[(i-1),j], sum(df16[j,]+df15[j,]*rho[(i-1),j])-df16[j,1]-df15[j,1]*rho[(i-1),j])
    t2[i,j] <- rbeta(1, df16[j,2]+df15[j,2]*rho[(i-1),j], sum(df16[j,]+df15[j,]*rho[(i-1),j])-df16[j,2]-df15[j,2]*rho[(i-1),j])
    t3[i,j] <- rbeta(1, df16[j,3]+df15[j,3]*rho[(i-1),j], sum(df16[j,]+df15[j,]*rho[(i-1),j])-df16[j,3]-df15[j,3]*rho[(i-1),j])
    sum <- sum(t1[i,j]+t2[i,j]+t3[i,j])
    t1[i,j] <- t1[i,j]/sum
    t2[i,j] <- t2[i,j]/sum
    t3[i,j] <- t3[i,j]/sum
  }

  for (j in 1:25){
    r[j] <- rho[(i-1),j] +runif(1, -0.01, 0.01)
    p1 <- dmultinom(c(df16[j,1]+df15[j,1]*r[j], df16[j,2]+df15[j,2]*r[j], df16[j,3]+df15[j,3]*r[j]), prob=c(t1[i,j], t2[i,j], t3[i,j]))
    po <- dmultinom(c(df16[j,1]+df15[j,1]*rho[(i-1),j], df16[j,2]+df15[j,2]*rho[(i-1), j], df16[j,3]+df15[j,3]*rho[(i-1),j]), prob=c(t1[i,j], t2[i,j], t3[i,j]))
    accp[i] <- min(1,p1/po)
    accept <- rbinom(1,1, min(1,p1/po))
    if (accept){
      rho[i,j]=r[j]
    }
    else{
      rho[i,j]=rho[(i-1),j]
    }
  }
  i=i+1
}


plot(t1[,1], type="l")
t1 <- as.data.frame(t1)
t2 <- as.data.frame(t2)
t3 <- as.data.frame(t3)
t <- na.omit(t1)
t2 <- na.omit(t2)
t3 <- na.omit(t3)

par(mfrow=c(2,2))
plot(t1$STATISTICIANS, type="l")
plot(t1$`MARKET RESEARCH ANALYSTS AND MARKETING SPECIALISTS`, type="l")
plot(t1$`FINANCIAL ANALYSTS`, type="l")
plot(t1$`DATABASE ADMINISTRATORS`, type="l")

par(mfrow=c(1,1))
p1 <- density(t1$STATISTICIANS)
p2 <- density(t1$`DATABASE ADMINISTRATORS`)
p3 <- density(t1$`FINANCIAL ANALYSTS`)
p4 <- density(t1$`MANAGEMENT ANALYSTS`)
p5 <- density(t1$`COMPUTER PROGRAMMERS`)
plot(p1, col=rgb(0.647059, 0.164706, 0.164706,0.8), xlim=c(0.60,0.98), ylim=c(0,80), main="Simulated Certified Rate of Different Jobs")
polygon(p1, col=rgb(0.647059, 0.164706, 0.164706,0.8), border="black")
lines(p5, col=rgb( 0.960784, 0.870588, 0.701961, 1/4))
polygon(p5, col=rgb( 0.960784, 0.870588, 0.701961, 1/4), border="black")
lines(p2, col=rgb(0.545098,0,0, 1/4))
polygon(p2, col=rgb(0.545098,0,0, 1/4), border="black")
lines(p3, col=rgb(0.545098, 0.270588, 0.0745098,1/4))
polygon(p3, col=rgb(0.545098, 0.270588, 0.0745098,1/4), border="black")
lines(p4, col=rgb(0.823529, 0.705882, 0.54902, 1/4))
polygon(p4, col=rgb(0.823529, 0.705882, 0.54902, 1/4), border="black")

legend("topleft", legend=c("Statsiticians", "Programmer","Database adm", "Financial Analyst", "Management Analyst"), 
       fill=c(rgb(0.647059, 0.164706, 0.164706,0.8), rgb( 0.960784, 0.870588, 0.701961, 1/4), rgb(0.545098,0,0, 1/4), rgb(0.545098, 0.270588, 0.0745098,1/4), rgb(0.823529, 0.705882, 0.54902, 1/4)))



quantile(t1$STATISTICIANS,probs=c(.025,.25,.5,.75,.975), na.rm=T)
quantile(t1$`FINANCIAL ANALYSTS`,probs=c(.025,.25,.5,.75,.975), na.rm=T)
quantile(t1$`MANAGEMENT ANALYSTS`,probs=c(.025,.25,.5,.75,.975), na.rm=T)
quantile(t1$`COMPUTER PROGRAMMERS`,probs=c(.025,.25,.5,.75,.975), na.rm=T)


for (i in 1:2000){
  yrep[i,] <- rmultinom(1, size=sum(VD16[24,]), prob=c(t1[i,24], t2[i,24], t3[i,24]))
}
hist(yrep[,1], main="Histogram of Replicated Data for Statistician", xlab="Statistician", breaks=20)
abline(v=VD16[24,1], col="brown4", lwd=3, lty=4)
sum(yrep[,1]>=VD16[24,1])/2000

for (i in 1:2000){
  yrep[i,] <- rmultinom(1, size=sum(VD16[5,]), prob=c(t1[i,5], t2[i,5], t3[i,5]))
}
hist(yrep[,1], main="Histogram of Replicated Data for Computer Programmer", xlab="Computer Programmer", breaks=20)
abline(v=VD16[5,1], col="brown4", lwd=3, lty=4)
sum(yrep[,1]>=VD16[5,1])/2000

for (i in 1:2000){
  yrep[i,] <- rmultinom(1, size=sum(VD16[13,]), prob=c(t1[i,13], t2[i,13], t3[i,13]))
}
hist(yrep[,1], main="Histogram of Replicated Data for Management Analyst", xlab="Management Analyst", breaks=20)
abline(v=VD16[13,1], col="brown4", lwd=3, lty=4)
sum(yrep[,1]>=VD16[13,1])/2000


