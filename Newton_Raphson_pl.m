function [result] = Newton_Raphson_pl(curr_den,J_0,a,Js,Us_ocp, U_ocp, a_s, J_pl, U_pl, R_film)
    %some constants
    F = 96485;
    T = 298;
    R = 8.314;
    
    %Solving for V
    %i_int + i_sei - i_tot = 0
    %i_int = 2*J_0*sinh(a*F*(V - U_ocp)/(R*T))
    %i_sei = -Js*exp(-a_s*F*(V - Us_ocp)/(R*T))
    %df/dt =  2*J_0*(a*F/(R*T))*cosh(a*F*(V - U_ocp)/(R*T)) + Js*(a_s*F/(R*T))*exp(-a_s*F*(V - Us_ocp)/(R*T))

    %{
    %TESTING
    curr_den = 1;
    J_0 = 3.870106843021043E-06;
    a = 0.5;
    Js = 0.75E-07;
    Us_ocp = 0.4;
    U_ocp = 3.63857020650;
    a_s = 0.5;
    %}

    V = U_ocp+0.5;
    
    while true
        val =  2*J_0*sinh(a*F*(V - U_ocp - R_film*curr_den)/(R*T)) - Js*exp(-a_s*F*(V - Us_ocp - R_film*curr_den)/(R*T)) + min(0, 2*J_pl*sinh(a*F*(V - U_pl- R_film*curr_den)/(R*T))) - curr_den;
        derivative = 2*J_0*(a*F/(R*T))*cosh(a*F*(V - U_ocp- R_film*curr_den)/(R*T)) + Js*(a_s*F/(R*T))*exp(-a_s*F*(V - Us_ocp - R_film*curr_den)/(R*T));
        if 2*J_pl*sinh(a*F*(V - U_pl- R_film*curr_den)/(R*T))<0
            derivative = derivative + 2*J_pl*(a*F/(R*T))*cosh(a*F*(V - U_pl- R_film*curr_den)/(R*T));
        end
        phi_new = V - (val/derivative);
        if norm(phi_new - V)<1e-10
            break
        end 
        V = phi_new;
    end
    
    V = phi_new;
    i_int = 2*J_0*sinh(a*F*(V - U_ocp - R_film*curr_den)/(R*T));
    i_sei = -Js*exp(-a_s*F*(V - Us_ocp - R_film*curr_den)/(R*T));
    i_pl = min(0, 2*J_pl*sinh(a*F*(V - U_pl- R_film*curr_den)/(R*T))); %might not ever give negative value
    result = [V, i_int, i_sei, i_pl];

%end