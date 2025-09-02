clc
clearvars
close all 

%ocv for li cobalt oxide vs x
%Universal values
cycles = 10; % Total number of cycles
cyc_vec = 1:2:cycles+1; % Vector representing cycles
F = 96485; % Faraday's constant
C_Li = 1000; % Lithium concentration
T = 298; % Temperature (Kelvin)
Rgc = 8.314; % Gas constant

% Time parameters
dt = 0.05; % Time step
t = 0:dt:3600; % Time vector
Nt = numel(t); % Number of time steps
voltage = zeros(Nt,1); % Initialize voltage vector

%CATHODE - constants

R_p = 8.5E-6; % Cathode radius
dr_p = R_p/20; % Cathode step size
r_p = 0:dr_p:R_p; % Cathode radius vector
D_p = 1E-14; % Diffusion coefficient for cathode
Nr_p = numel(r_p); % Number of cathode steps
Cmax_p = 51410; % Maximum lithium concentration in cathode
curr_den_p = 1.5; % Initial current density for cathode
a_p = 0.5; % Alpha parameter for cathode
k_p = 6.67E-11; % Rate constant for cathode


%ANODE - constants

R_n = 12.5E-6; % Anode radius
dr_n = R_n/20; % Anode step size
r_n = 0:dr_n:R_n; % Anode radius vector
D_n = 3.9E-14; % Diffusion coefficient for anode
Nr_n = numel(r_n); % Number of anode steps
Cmax_n = 31833; % Maximum lithium concentration in anode
curr_den_n = curr_den_p; % Initial current density for anode
a_n = 0.5; % Alpha parameter for anode
k_n = 1.764E-11; % Rate constant for anode


% Parameters for parasitic reaction
Us_ocp = 0.4; % Open circuit potential
Js = 0.75E-2; % Exchange current density
a_s = 0.5; % Alpha parameter for parasitic reaction
rho_sei = 1690; % SEI density
cond_sei = 5E-6; % SEI conductivity
M_sei = 0.162; % SEI molar mass

% Parameters for lithium plating
a_pl = 0.5; % Alpha parameter for plating
U_pl = 0; % Plating potential
c_e = 1000; % Lithium concentration
k_pl = 2.23E-7; % Rate constant for plating
rho_pl = 534; % Plating density
cond_pl = 1.1E7; % Plating conductivity
M_pl = 0.007; % Plating molar mass
J_pl = F*k_pl*(c_e^(0.5)); % Plating current density

%Parameters for lithium stripping
m_pl = zeros(Nt,1); % Plating thickness vector

%CATHODE - variables
C_p = zeros(Nr_p,Nt); % Cathode concentration matrix
over_pot_p = zeros(Nt,1); % Overpotential for cathode
U_ocp_p = zeros(Nt,1); % Open circuit potential for cathode
phi_p = zeros(Nt,1); % Electric potential for cathode

%ANODE - variables
C_n = zeros(Nr_n,Nt); % Anode concentration matrix
over_pot_n = zeros(Nt,1); % Overpotential for anode
U_ocp_n = zeros(Nt,1); % Open circuit potential for anode
phi_n = zeros(Nt,1); % Electric potential for anode
SOC = zeros(Nt,1); % State of charge

%Parasitic reaction - variables
i_tot = zeros(cycles,1); % Total current vector
i_int = zeros(Nt,1); % Intercalation current vector
i_sei = zeros(Nt,1); % SEI current vector
i_pl = zeros(Nt,1); % Plating current vector
i_str = zeros(Nt,1); % Stripping current vector
Q_lost = zeros(cycles,1); % Lost charge vector

%Resistance model - variables
t_film = zeros(cycles,1); % Film thickness vector
t_sei = zeros(cycles,1); % SEI thickness vector
t_li = zeros(cycles,1); % Lithium thickness vector

%Initialisation
C_p(:,1) = 0.5*Cmax_p; % Initial cathode concentration
C_n(:,1) = 0.95*Cmax_n; % Initial anode concentration
i_tot(1,1) = curr_den_p; % Initial total current
t_film(1,1) = 10E-9; % Initial film thickness
t_sei(1,1) = 10E-9; % Initial SEI thickness

% Start of simulation
counter = 0; % Counter for cycle type
for cycle = 1:cycles
    % Current density remains constant for one cycle

    % Determine type of cycle
    if counter == 0 % Discharge cycle
        counter = 1; 
        if curr_den_p<0
            curr_den_p = i_tot(cycle,1);
            curr_den_n = i_tot(cycle,1);
        end
    
    elseif counter == 1 % Charge cycle
        counter = 0;
        if curr_den_p>0
            curr_den_p = -i_tot(cycle,1);
            curr_den_n = -i_tot(cycle,1);
        end
    end

    
    % Propagation loop
    for j=1:Nt % Time loop
        
        if j~=1 % Not at boundary
            % Cathode concentration update
            for i=2:Nr_p-1 % space loop for anode
                C_p(i,j)=C_p(i,j-1)+D_p*dt*(((C_p(i+1,j-1)-2*C_p(i,j-1)+C_p(i-1,j-1))/(dr_p)^2)+((C_p(i+1,j-1)-C_p(i-1,j-1))/(r_p(i)*dr_p)));
            end
            % Anode concentration update
            for i=2:Nr_n-1 % space loop for cathode
                C_n(i,j)=C_n(i,j-1)+D_n*dt*(((C_n(i+1,j-1)-2*C_n(i,j-1)+C_n(i-1,j-1))/(dr_n)^2)+((C_n(i+1,j-1)-C_n(i-1,j-1))/(r_n(i)*dr_n)));
            end
    
            % Boundary conditions    
            C_p(1,j)=C_p(2,j);
            C_p(Nr_p,j)=C_p(Nr_p-1,j)+curr_den_p*dr_p/(D_p*F); 
        
            C_n(1,j)=C_n(2,j);
            C_n(Nr_n,j)=C_n(Nr_n-1,j)-curr_den_n*dr_n/(D_n*F);
        end
        
        % Some variables for current calculation
        J_0_p = k_p*((Cmax_p-C_p(Nr_p,j))*C_p(Nr_p,j)*C_Li)^a_p;  
        J_0_n = k_n*((Cmax_n-C_n(Nr_n,j))*C_n(Nr_n,j)*C_Li)^a_n;
        SOC_p = C_p(Nr_p,j)/Cmax_p;
        SOC_n = C_n(Nr_n,j)/Cmax_n;
        SOC(j,1) = SOC_n;
        U_ocp_p(j,1) = Uocp_interp2(SOC_p); %:)
        U_ocp_n(j,1) = Uocp_interp(SOC_n);
 
        cond_film = cond_sei;
        R_film = t_film(cycle,1)/cond_film;

        if curr_den_p <= 0 % Charge cycle
            
            if j~=1 && phi_n(j-1,1)< 0
                % Plating
                
                result = Newton_Raphson_pl(curr_den_n, J_0_n, a_n, Js, Us_ocp, U_ocp_n(j,1),a_s, J_pl, U_pl, R_film);
                phi_n(j,1) = result(1,1);
                i_int(j,1) = result(1,2);
                i_sei(j,1) = result(1,3)*10; 
                i_pl(j,1) = result(1,4); 
                
            
            elseif j~=1 && phi_n(j-1,1) > 0 && max(m_pl) > 0
                %Stripping
                
                result = Newton_Raphson_str(0, J_0_n, a_n, Js, Us_ocp, U_ocp_n(j,1),a_s, J_pl, R_film, m_pl(j-1,1), max(m_pl)); %CHANGE
                phi_n(j,1) = result(1,1);
                i_int(j,1) = result(1,2);
                i_sei(j,1) = result(1,3)*10;
                i_str(j,1) = result(1,4); 
           
            else
            
            
                result = Newton_Raphson2(curr_den_n, J_0_n, a_n, Js, Us_ocp, U_ocp_n(j,1),a_s);
                phi_n(j,1) = result(1,1);
                i_int(j,1) = result(1,2);
                i_sei(j,1) = result(1,3)*10; 
            
            end
            %phi_n(j,1) = (Rgc*T/(a_n*F))*asinh(i_tot(cycle,1)/(2*J_0_n)) + U_ocp_n(j,1);
            over_pot_p(j,1) = (Rgc*T/(a_p*F))*asinh(i_tot(cycle,1)/(2*J_0_p));
            phi_p(j,1) = over_pot_p(j,1) + U_ocp_p(j,1);
    
            % Calculate cell voltage
            voltage(j,1) = phi_p(j,1) - phi_n(j,1);
            %Loss of Li ions
            Q_lost(cycle,1) = Q_lost(cycle,1) + i_sei(j,1)*dt + i_pl(j,1)*dt + i_str(j,1)*dt; %negative
            
            % Thickness of SEI layer + plating/stripping
            t_sei(cycle,1) = t_sei(cycle,1) - (i_sei(j,1)*M_sei*dt)/(F*rho_sei);
            t_li(cycle,1) = t_li(cycle,1) - (i_pl(j,1)*M_pl*dt)/(F*rho_pl);
            t_film(cycle,1) = t_sei(cycle,1) + t_li(cycle,1); 
            m_pl(j,1) = 4*pi*R_n*R_n*t_li(cycle,1)*rho_pl; 
            
        elseif curr_den_p > 0 % Discharge cycle

            % Calculate overpotentials
            over_pot_p(j,1) = (Rgc*T/(a_p*F))*asinh(curr_den_p/(2*J_0_p));
            over_pot_n(j,1) = (Rgc*T/(a_n*F))*asinh(curr_den_n/(2*J_0_n));
            
            % Calculate electric potentials
            phi_p(j,1) = over_pot_p(j,1) + U_ocp_p(j,1);
            phi_n(j,1) = over_pot_n(j,1) + U_ocp_n(j,1);
            
            % Calculate cell voltage
            voltage(j,1) = phi_p(j,1) - phi_n(j,1);
    
        end
        

        % Check for voltage cutoff
        
        if (voltage(j,1)<3 || SOC_n<0.05) && counter == 1 %Discharge cutoff
            break
        elseif  (voltage(j,1)>4.2 || SOC_p<0.5) && counter == 0 %Charge cutoff
            break
        end
        
    end
    
    % Reinitialization
    i_tot(cycle+1,1) = i_tot(cycle,1) + Q_lost(cycle,1)/3600;
    t_film(cycle+1,1) = t_film(cycle,1);
    t_sei(cycle+1,1) = t_sei(cycle,1);
    t_li(cycle+1,1) = t_li(cycle,1);
    m_pl = zeros(Nt,1);
    C_p(:,1) = C_p(Nr_p, j); 
    C_n(:,1) = C_n(Nr_n, j); 

    % Plot variation of voltage with time
    %{
    figure(1)
    timevec = 0+(cycle-1)*3600:dt:3600+(cycle-1)*3600; 
    plot(timevec,i_sei,'-');
    xlabel("Time (s)")
    ylabel("Cell Voltage (V)")
    pause(1)
    hold on
    %}

    % Plot variation of voltage with SOC
    %{
    if mod(cycle,2) == 0
        figure(2)
        timevec2 = transpose(0:dt:3600);
        answer = SOC./timevec2;
        plot(answer*abs(curr_den_p), voltage);
        title("Cell Voltage vs SOC")
        xlabel("SOC")
        ylabel("Cell Voltage (V)")
        pause(1)
        hold on
    end
    %}
    
end

%{
figure(3)
plot(t_sei(2:2:end));
figure(4)
plot(i_tot(2:2:end));
%}

