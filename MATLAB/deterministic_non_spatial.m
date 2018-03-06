%Solving the differential equations and then plotting the solution. 

%Input variables
Vmax = 30;
Km = 30;
k1 = 0.00;
k2 = 0.0;

%Initial conditions
D1_0 = 300; %Constant
DAout_0 = 3000; DAD1_0 = 10; 

%ODE45 options
tstart = 0;  tend = 100;
abstol = 1e-8; reltol = 1e-8;
options = odeset('RelTol',reltol,'AbsTol',abstol);

%Applying the MATLAB scheme ode45 to solve the DEs
[t,Y] = ode45(@synapse,[tstart tend],[DAout_0 DAD1_0],options,D1_0,Vmax,Km,k1,k2);

%Plotting the results
hold on
plot(t,Y(:,1))
plot(t,Y(:,2))
