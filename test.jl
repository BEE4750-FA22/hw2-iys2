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

#Find initial concentration at Water Source 1

function int_conditions(inflow1, inflow2, do1, do2, cbod1, cbod2, nbod1, nbod2)
   c0 = (inflow1 * do1 + inflow2 * do2)/(inflow1+inflow2)
   b0 = (inflow1 * cbod1 + inflow2 * cbod2)/(inflow1+inflow2)
   n0 = (inflow1*nbod1 + inflow2*nbod2)/(inflow1+inflow2)

   return [c0, b0, n0]
   
end

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

c = zeros(51)

function total_do(cbod_river, nbod_river, cbod_1, nbod_1, cbod_2, nbod_2)
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
min = 0
count = 0
while min < 4000
      min = total_do(cbod_river,nbod_river,cbod_1,nbod_1,(1-count)*cbod_2,(1-count)*nbod_2)[2]
      count = count + 0.001
end

#Q4
min = 0
count = 0
while min < 4000
      min = total_do(cbod_river,nbod_river,(1-count)*cbod_1,(1-count)*nbod_1,(1-count)*cbod_2,(1-count)*nbod_2)[2]
      count = count + 0.001
end


#Q5
both 

#Q6
k = zeros(100*100)
n = 100
total = 0
success_count = 0
g = sample_correlated_uniform(n, [4000,7000],[3000,8000],0)
for i = 1:n
   for j = 1:n
      min = total_do(g[i,1],g[j,2],(1-count)*cbod_1,(1-count)*nbod_1,(1-count)*cbod_2,(1-count)*nbod_2)[2]
      if min > 4000
         success_count = success_count +1
      end
      total = total +1
      k[total] = success_count/total
      println(success_count/total)
   end

end

#Q7
p = zeros(1000*1000)
n = 1000
total = 0
success_count = 0
g = sample_correlated_uniform(n, [4000,7000],[3000,8000])
for i = 1:n
   for j = 1:n
      min = total_do(g[i,1],g[j,2],(1-count)*cbod_1,(1-count)*nbod_1,(1-count)*cbod_2,(1-count)*nbod_2)[2]
      if min > 4000
         success_count = success_count +1
      end
      total = total +1
      p[total] = success_count/total
      println(min, "     ", success_count/total)
   end

end