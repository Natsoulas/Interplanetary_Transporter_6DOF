%pwpf matlab (from simulink design)
%model based off of Trond Dagfinn Krovel's Paper: "Optimal Selection of PWPF Parameters"
function [u,DC,f_o] = PWPF_Run(C,K_p,K_m,T_m,U_on,U_off,t_s,timesteps,tol)
%initialize
f_zero = 0;
U_initial = 0;
b = 0;
f = [];
u = [];
r = C*ones(1,timesteps);
uconcat = U_initial;
u = [u uconcat];
counter = 1;
for t = linspace(0,t_s,timesteps)
f_time = f_zero + (K_m*(K_p*r(counter)-u(counter))-f_zero)*(1 - exp(-(t-b)/T_m));
f =[f f_time];
 if f(counter) < U_off
     uconcat = 0;
 elseif f(counter) > U_on - tol && f(counter) < U_on + tol && u(counter) == 0
     uconcat = 1;
     f_zero = f(counter);
     b = t;
 elseif f(counter) > U_off - tol && f(counter) < U_off + tol && u(counter) == 1
     uconcat = 0;
     f_zero = f(counter);
     b = t;
 else
     uconcat = u(end);
 end
 u = [u uconcat];
 counter = counter + 1;
 end
u = u(1:end-1);
t = linspace(0,t_s,timesteps);
h = U_on - U_off;
Ton = -T_m*log(1- h/(U_on - K_m*(C-u(end))));
Toff = -T_m*log(1-h/(K_m*C-U_off));
DC = Ton/(Ton + Toff);
f_o = 1/(Ton + Toff);
%figure
%plot(t,r,t,f,t,u)
end