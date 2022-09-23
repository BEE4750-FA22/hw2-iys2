#Q1

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

#Find initial concentration at Water Source 1

c0 = (q_river * do_river + q_1 * do_1)/(q_river+q_1)
b0 = (q_river * cbod_river + q_1 * cbod_1)/(q_river+q_1)
n0 = (q_river*nbod_river + q_1*nbod_1)/(q_river+q_1)
cs = 10000 #mg/m^3

function dissolved_ox(u, c, cs, c0, b0, n0, ka, kc, kn, x1, x2)   
   
   
   for i = x1:x2
    a1 = exp(-ka*i/u)
    a2 = (kc/(ka-kc))*(exp(-kc*i/u)-exp(-ka*i/u))
    a3 = (kn/(ka-kn))*(exp(-kn*i/u)-exp(-kn*i/u))

    c[i+1] = cs*(1-a1)+c0*a1-b0*a2-n0*a3

   end


end

c = zeros(51)

d = dissolved_ox(6, c, cs, c0, b0, n0, 0.55, 0.35, 0.25, 0, 15)


plot(c)