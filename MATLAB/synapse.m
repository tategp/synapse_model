function dYdt = synapse(t,Y,D1_0,Vmax,Km,k1,k2)

%Set of equations for deterministic_non_spatial.m

dYdt = zeros(2,1);

%Y(1)=[DA]_out, Y(2)=[DAD1]
dYdt(1) = -Vmax*Y(1)/(Km+Y(1))-k1*D1_0*Y(1)+k2*Y(2);
dYdt(2) = k1*D1_0*Y(1)-k2*Y(2);
