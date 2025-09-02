clc
clearvars
close all 

%Do we change both current density and concentration from next cycle?

%positive electrode(anode)
R = 2.5E-6;
dr = R/20;
r = 0:dr:R;
dt = 0.1;
t = 0:dt:3600;
D = 1E-14; %from paper
Nr = numel(r); Nt = numel(t);
Cmax = 30555;
Ci = zeros(11,1); % 10 number of cycles
Ci(1,1) = 0.95*Cmax; %from paper
curr_den = 1.5;
F = 96485;
a = 0.5; %alpha
C_Li = 1000;
k = 1E-7; %MAN MADE VALUE TO AGREE
T = 298;
Rgc = 8.314;


J_sei = zeros(11,1);
J = zeros(11,1);
J(1,1) = curr_den;
Q_curr = zeros(11,1);

for cycle = 1:10
    
    C = zeros(Nr,Nt);
    C(:,1) = Ci(cycle,1);
    curr_den = J(cycle,1);
    over_pot = zeros(Nt,1);
    U_ocp = zeros(Nt,1);
    phi = zeros(Nt,1);

    for j=2:Nt %time loop
        for i=2:Nr-1 %space loop
            C(i,j)=C(i,j-1)+D*dt*(((C(i+1,j-1)-2*C(i,j-1)+C(i-1,j-1))/(dr)^2)+((C(i+1,j-1)-C(i-1,j-1))/(r(i)*dr)));
        end
        C(1,j)=C(2,j);
        C(Nr,j)=C(Nr-1,j)-curr_den*dr/(D*F);
        
        J_0 = k*((Cmax-C(Nr,j))*C(Nr,j)*C_Li)^a;
        over_pot(j,1) = (Rgc*T/(a*F))*asinh(curr_den/(2*J_0));  
        SOC = C(Nr,j)/Ci(cycle,1);
        U_ocp(j,1) = Uocp_interp(SOC);
        phi(j,1) = over_pot(j,1) + U_ocp(j,1);
        if C(Nr,j)/Cmax < 0.50
            break
        end
        
    end
    
    %calculated using above values unlike research paper
    J_sei(cycle,1) = -J_0*(exp(-a*F*over_pot(Nt,1)/(Rgc*T)));
    J(cycle+1,1) = J(cycle,1) + J_sei(cycle,1);
    %problem with units here 
    %also no integration
    Ci(cycle+1,1) = J(cycle+1,1)*3600; %time for one cycle

end

cycles = 1:1:11;
%PLOTTING
plot(cycles, Ci/Ci(1,1), "-");
title("Capacity Fade")


%{
figure(1)
plot(t,over_pot, '-')
title("Over-potential")

figure(2)
plot(t,U_ocp,'-')
title("Open circuit potential")

figure(3)
plot(t,phi,'-')
title("potential at anode")

%{
figure(1)
plot(t, C(1,:), '-', t, C(50,:), '-', t, C(100,:),'-')
title('Concentration')
xlabel('R')
ylabel('C')
%}

% mesh(C)
%}
