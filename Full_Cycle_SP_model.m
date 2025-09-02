
%I AM SO FUCKING CONFUSED
clc
clearvars
close all 

%Do we change both current density and concentration from next cycle?

F = 96485;
C_Li = 1000;
T = 298;
Rgc = 8.314;

%time loop
dt = 0.05;
t = 0:dt:3600;
Nt = numel(t);
voltage = zeros(Nt,1);

%CATHODE
%changed R values according to paper *****************
R_p = 8.5E-6;
dr_p = R_p/20;
r_p = 0:dr_p:R_p;
D_p = 1E-14; %from paper
Nr_p = numel(r_p);
C_p = zeros(Nr_p,Nt);
Cmax_p = 51410;
curr_den_p = 1;
a_p = 0.5; %alpha
k_p = 6.67E-11;
over_pot_p = zeros(Nt,1);
U_ocp_p = zeros(Nt,1);
phi_p = zeros(Nt,1); 

%ANODE
R_n = 12.5E-6;
dr_n = R_n/20;
r_n = 0:dr_n:R_n;
D_n = 3.9E-14; %from paper
Nr_n = numel(r_n);
C_n = zeros(Nr_n,Nt);
Cmax_n = 31833;
curr_den_n = curr_den_p;
a_n = 0.5; %alpha
k_n = 1.764E-11;
over_pot_n = zeros(Nt,1);
U_ocp_n = zeros(Nt,1);
phi_n = zeros(Nt,1);
J_sei = zeros(11,1);
J_n = zeros(11,1);
J_n(1,1) = curr_den_n;
Q_curr = zeros(11,1);


SOC = zeros(Nt,1);
%ANODE CATHODE CHARGE DISCHARGE CYCLE - ASK IF THIS IS CORRECT
Ci_p = 0.5*Cmax_p;
Ci_n = 0.95*Cmax_n; 

Q(1,1) = Ci_n;

counter = 0;
for cycle = 1:10
    
    %change in time - could make change directly in plot
    %plotting
    %voltage of anode, cathode, cell
    C_p = zeros(Nr_p, Nt);
    C_n = zeros (Nr_n, Nt);
    C_p(:,1) = Ci_p; %CHANGE %Ci_p
    C_n(:,1) = Ci_n; %CHANGE

    if counter == 0
        counter = 1;
        if curr_den_p<0
            curr_den_p = -curr_den_p;
            curr_den_n = -curr_den_n;
        end
    
    elseif counter == 1
        counter = 0;
        if curr_den_p>0
            curr_den_p = -curr_den_p;
            curr_den_n = -curr_den_n;
        end
    end

    %for the boundary elements
    J_0_p = k_p*((Cmax_p-C_p(Nr_p,1))*C_p(Nr_p,1)*C_Li)^a_p;  
    J_0_n = k_n*((Cmax_n-C_n(Nr_n,1))*C_n(Nr_n,1)*C_Li)^a_n;
    var = J_0_p;
    var2 = J_0_n;

    over_pot_p(1,1) = (Rgc*T/(a_p*F))*asinh(curr_den_p/(2*J_0_p));
    over_pot_n(1,1) = (Rgc*T/(a_n*F))*asinh(curr_den_n/(2*J_0_n));

    SOC_p = C_p(Nr_p,1)/Cmax_p;
    SOC_n = C_n(Nr_n,1)/Cmax_n;
    SOC(1,1) = SOC_p;

    U_ocp_p(1,1) = Uocp_interp(SOC_p); 
    U_ocp_n(1,1) = Uocp_interp2(SOC_n);

    phi_p(1,1) = over_pot_p(1,1) + U_ocp_p(1,1);
    phi_n(1,1) = over_pot_n(1,1) + U_ocp_n(1,1); %this function will be different for the cathode
    
    voltage(1,1) = phi_n(1,1) - phi_p(1,1);
    
    for j=2:Nt %time loop
    
        for i=2:Nr_p-1 %space loop for anode
            C_p(i,j)=C_p(i,j-1)+D_p*dt*(((C_p(i+1,j-1)-2*C_p(i,j-1)+C_p(i-1,j-1))/(dr_p)^2)+((C_p(i+1,j-1)-C_p(i-1,j-1))/(r_p(i)*dr_p)));
        end
    
        for i=2:Nr_n-1 %space loop for cathode
            C_n(i,j)=C_n(i,j-1)+D_n*dt*(((C_n(i+1,j-1)-2*C_n(i,j-1)+C_n(i-1,j-1))/(dr_n)^2)+((C_n(i+1,j-1)-C_n(i-1,j-1))/(r_n(i)*dr_n)));
        end
    
        C_p(1,j)=C_p(2,j);
        C_p(Nr_p,j)=C_p(Nr_p-1,j)+curr_den_p*dr_p/(D_p*F);  
        
        C_n(1,j)=C_n(2,j);
        C_n(Nr_n,j)=C_n(Nr_n-1,j)-curr_den_n*dr_n/(D_n*F);  
        
        J_0_p = k_p*((Cmax_p-C_p(Nr_p,j))*C_p(Nr_p,j)*C_Li)^a_p;  
        J_0_n = k_n*((Cmax_n-C_n(Nr_n,j))*C_n(Nr_n,j)*C_Li)^a_n;  

        over_pot_p(j,1) = (Rgc*T/(a_p*F))*asinh(curr_den_p/(2*J_0_p));
        over_pot_n(j,1) = (Rgc*T/(a_n*F))*asinh(curr_den_n/(2*J_0_n));
    
        SOC_p = C_p(Nr_p,j)/Cmax_p;
        SOC(j,1) = SOC_p;
        SOC_n = C_n(Nr_n,j)/Cmax_n;
    
        U_ocp_p(j,1) = Uocp_interp(SOC_p); %change functions
        U_ocp_n(j,1) = Uocp_interp2(SOC_n);
    
        phi_p(j,1) = over_pot_p(j,1) + U_ocp_p(j,1);
        phi_n(j,1) = over_pot_n(j,1) + U_ocp_n(j,1);
        
        voltage(j,1) = phi_n(j,1) - phi_p(j,1);


        %might not be correct
        %{
        if (voltage(j,1)<3 || SOC_n<0.05) && counter == 1
            %C_p(Nr_p,j) = 0;
            break
        elseif  (voltage(j,1)>4.2 || SOC_p<0.5) && counter == 0
            %C_n(Nr_n,j) = 0;
            break
        end
        %}
        
        
        
    end
    
    figure(1)
    timevec = 0+(cycle-1)*3600:dt:3600+(cycle-1)*3600; 
    axis([0 36000 3  5])
    plot(timevec,voltage,'-');
    xlabel("Time")
    ylabel("Cell Voltage")
    pause(1)
    hold on
    %{
    if mod(cycle,2) == 0
        figure(2)
        plot(SOC*curr_den_p*3600/dt, voltage);
        title("Cell Voltage vs Capacity")
        xlabel("Capacity (Ah)")
        ylabel("Cell Voltage (V)")
        hold on
    end
    %}
    %{
    figure(2);
    axis([0 36000 -5 5])
    plot(timevec, phi_p,'-');
    pause(0.5)
    hold on
    %}

    %values to start next cycle
    Ci_p = C_p(Nr_p, Nt);
    Ci_n = C_n(Nr_n, Nt);

end
%}



