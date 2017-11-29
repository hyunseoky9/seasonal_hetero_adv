# This model seeks for heterozygosity advantage in shifting seasons
# through simulating population frequency each time step. 
# depending on the expression level difference due to mutation allele, the 
# survival probability changes along the normal curve.
# In dominancy4, a new, third mutation is introduced to the population
# at time stip 100 with a random change in expression level.
#set.seed(1001)
rm(list=ls())
x = c(499/500,1/500) #frequency
freq = c(x)
t=1000 # generation amount
s = 1 # season either 1 or 2.
switch = 5 #generation amount after a season change.
u1 = 0.333; u2 = 0.666; sig = 0.1665 #means of normal curves for different seasons.
a1 = runif(1,0,0.5) #expression level of a wildtype
d = ((-2*a1+u1+u2)/3 + (-2*a1+u1+u2))/2
#a1 = 0.33; d = 0.17
exp_d = c(0,d)
a = c(a1, a1+d, a1+2*d)
w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig)
w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
w = w1
new_mut = 100 #time step when new mutation comes in

for(i in 1:new_mut){ #loooping time step
  if(i%%switch == 1 && i>1){ # change season after specified generation time
    if(s==1){
      s=2
      w = w2
    } else {
      s=1
      w = w1
    }
  }
  wbar = w[1]*x[1]^2 + w[2]*2*x[1]*x[2] + w[3]*x[2]^2
  x[2] = (w[3]*x[2]^2 + w[2]*x[1]*x[2])/wbar
  x[1] = (w[2]*x[1]*x[2]+w[1]*x[1]^2)/wbar
  freq = cbind(freq,x)
}

#NEW MUTATION

d2 = runif(1,0,0.5) # new mutation expression dist.
exp_d = c(exp_d,d2)
a = c()
alleles_list = c() #keeps record of which index is what genotype
#allele combination orders stored
for(i in 1:length(exp_d)){ #refill all the expression levels
  for(j in i:length(exp_d)){
    a = c(a, a1+exp_d[i]+exp_d[j])
    alleles_list = cbind(alleles_list,c(i,j))
  }
}
w1 = dnorm(a,u1,sig)/dnorm(u1,u1,sig) #remake all the w1&w2
w2 = dnorm(a,u2,sig)/dnorm(u2,u2,sig)
w = w1

#put in frequency for third allele. New mutation arise from
#wild-type
x=c(x[1]-1/500,x[2:length(x)],1/500)
#put all the past freq for third allele (which is 0)
freq = rbind(freq,rep(0,dim(freq)[2]))
freq[,ncol(freq)] = x

#time steps after 3rd mutation
for(i in (new_mut+1):t){ #loooping time step
  if(i%%switch == 1){ # change season after specified generation time
    if(s==1){
      s=2
      w = w2
    } else {
      s=1
      w = w1
    }
  }
  x_square = c() 
  for(j in 1:ncol(alleles_list)){ #get all the factors when x is squared
    if(alleles_list[1,j]==alleles_list[2,j]){
      x_square = c(x_square, x[alleles_list[1,j]]*x[alleles_list[2,j]])
    } else{
      x_square = c(x_square, 2*x[alleles_list[1,j]]*x[alleles_list[2,j]])
    }
  }
  wbar = sum(x_square*w)
  for(j in 1:length(x)){ #update each frequency in x
    factor_sum = 0
    for(k in 1:ncol(alleles_list)){ #for every genotype
      if(all(alleles_list[,k]==j)){
        factor_sum = factor_sum + x_square[k]*w[k]
      } else if(any(alleles_list[,k]==j)){
        factor_sum = factor_sum + 0.5*x_square[k]*w[k]
      }
    }
    x[j] = factor_sum/wbar
  }
  freq = cbind(freq,x)
}

#result output
par(mfrow=c(3,1))
par(mar=c(2,2,2,2))
title = sprintf('bl=1,r=2,g=3')
time_line = 100:1000
plot(freq[1,time_line], ylim=c(0,1), type='l',main = title) #plot of frequency
lines(freq[2,time_line],type='l',col=2)
lines(freq[3,time_line],type='l',col=3)

allele_label = c()
for( i in 1:ncol(alleles_list)){
  label = sprintf('A%s%s',alleles_list[1,i],alleles_list[2,i])
  allele_label = c(allele_label, label)
}
ww = w1*w2
wmean = (w1+w2)/2
plot(wmean)
plot(ww, axes=FALSE, ylim=c(0,max(ww)), col=1) #plot of w1*w2
axis(2)
axis(1, at=seq_along(ww),labels=allele_label)
title('harmonic mean of two seasons (w1*w2)')
#points(w1,col=1)
#points(w2,col=2)
#legend(5,0.6,legend=c('season1 w','season2 w'),pch=1,col=c("black",'red'), cex=0.8)

sum(ww[1],ww[2]/2,ww[3]/2)
sum(ww[2]/2,ww[4],ww[5]/2)
sum(ww[3]/2,ww[5]/2,ww[6])

sum(wmean[1],wmean[2]/2,wmean[3]/2)
sum(wmean[2]/2,wmean[4],wmean[5]/2)
sum(wmean[3]/2,wmean[5]/2,wmean[6])

sum(ww[c(1,2,3)])
sum(ww[c(2,4,5)])
sum(ww[c(3,5,6)])

sum(wmean[c(1,2,3)])
sum(wmean[c(2,4,5)])
sum(wmean[c(3,5,6)])

print(ww[2])
