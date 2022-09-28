#Q1

using Plots
#Initial concentrations
q_river = 100000 #m^3/d
q_1 = 10000 #m^3/d
q_2 = 15000 #m^3/d

do_river = 7500 #mg/m^3
do_1 = 5000 #mg/m^3
do_2 = 5000 #mg/m^3

cbod_river = 5000 #mg/m^3
cbod_1 = 50000 #mg/m^3
cbod_2  = 45000 #mg/m^3

nbod_river = 5000 #mg/m^3
nbod_1 = 35000 #mg/m^3
nbod_2 = 35000 #mg/m^3

cs = 10000 #mg/m^3

#Finding initial concentrations given two inflows with DO, CBOD, and NBOD concentrations

function int_conditions(inflow1, inflow2, do1, do2, cbod1, cbod2, nbod1, nbod2)
   c0 = (inflow1 * do1 + inflow2 * do2)/(inflow1+inflow2)
   b0 = (inflow1 * cbod1 + inflow2 * cbod2)/(inflow1+inflow2)
   n0 = (inflow1*nbod1 + inflow2*nbod2)/(inflow1+inflow2)

   return [c0, b0, n0]
   
end

#Finding dissolved oxygen concentration as a function of distance and final CBOD and NBOD concentrations given initial conditions and decay rates. 

function dissolved_ox(u, c, cs, c0, b0, n0, ka, kc, kn, x1, x2)   
   
   
   for i = 0:x2-x1
    a1 = exp(-ka*i/u)
    a2 = (kc/(ka-kc))*(exp(-kc*i/u)-exp(-ka*i/u))
    a3 = (kn/(ka-kn))*(exp(-kn*i/u)-exp(-ka*i/u))
   
    c[x1+i+1] = (cs*(1-a1))+(c0*a1)-(b0*a2)-(n0*a3)

   end
   
   b = b0*exp(-kc*x2/u)
   n = n0*exp(-kn*x2/u)

   return [c[x2],b,n]

end



function total_do(cbod_river, nbod_river, cbod_1, nbod_1, cbod_2, nbod_2)
   c = zeros(51)

   conc_1 = int_conditions(100000,10000,7500,5000,cbod_river,cbod_1,nbod_river,nbod_1)
   d = dissolved_ox(6, c, cs, conc_1[1], conc_1[2], conc_1[3], 0.55, 0.35, 0.25, 0, 15)

   conc_2 = int_conditions(110000,15000,d[1],do_2,d[2],cbod_2,d[3],nbod_2)
   dissolved_ox(6,c,cs,conc_2[1],conc_2[2],conc_2[3],0.55,0.35,0.25,15,50)

   return [c, minimum(c)]

end


plot([0:50],total_do(cbod_river,nbod_river, cbod_1,nbod_1,cbod_2,nbod_2)[1], label="Dissolved Oxygen (mg/m^3)")
vline!([15], label="Discharge")
hline!([6000], label = "Regulation")

#Q3
min_2_treated = 0;
t = 0;
while min_2_treated < 4000
      global min_2_treated = total_do(cbod_river,nbod_river,cbod_1,nbod_1,(1-t)*cbod_2,(1-t)*nbod_2)[2]
      global t = t + 0.001
end

println(t)

#Q4
min_both_treated = 0;
t = 0;
while min_both_treated < 4000
      global min_both_treated = total_do(cbod_river,nbod_river,(1-t)*cbod_1,(1-t)*nbod_1,(1-t)*cbod_2,(1-t)*nbod_2)[2]
      global t = t + 0.001
end

println(t)


#Q5 

#Q6
using Distributions
using PDMats
function sample_correlated_uniform(n, a, b, corr_coef=0.7)
   mvnorm = MvNormal([0, 0], PDMat([1 corr_coef; corr_coef 1])) # set up a multivariate normal with each marginal variance of 1 and the right correlation
   norm_samples = rand(mvnorm, n)' # sample from the multivariate normal, the marginal distributions are a standard normal
   unif_samples = cdf.(Normal(0, 1), norm_samples) # convert samples to a uniform distribution using the pdf of a standard Normal
   samples = (unif_samples .* [a[2] - a[1] b[2] - b[1]]) .+ [a[1] b[1]]
   return samples
 end
 
 
n = 1000;
k = zeros(n^2);
total1 = 0;
success_count1 = 0;
for i = 1:n
   for j = 1:n
      r1 = rand(Uniform(4000,7000))
      r2 = rand(Uniform(3000,8000))
      min1 = total_do(r1,r2,(1-t)*cbod_1,(1-t)*nbod_1,(1-t)*cbod_2,(1-t)*nbod_2)[2]
      if min1 > 4000
         global success_count1 = success_count1 +1
      end
      global total1 = total1 +1
      k[total1] = success_count1/total1
   end

end

println(k[total1])
plot(k)

#Q7
n = 1000;
p = zeros(n^2);
total2 = 0;
success_count2 = 0;
g = sample_correlated_uniform(n, [4000,7000],[3000,8000]);
for i = 1:n
   for j = 1:n
      min2 = total_do(g[i,1],g[j,2],(1-t)*cbod_1,(1-t)*nbod_1,(1-t)*cbod_2,(1-t)*nbod_2)[2]
      if min2 > 4000
         global success_count2 = success_count2 +1
      end
      global total2 = total2 +1
      p[total2] = success_count2/total2
   end

end

println(p[total2])
plot(p)

diff = zeros(n^2)
for i = 1:length(p)
   diff[i] = k[i]-p[i]
end

println(diff[total])
plot(diff)