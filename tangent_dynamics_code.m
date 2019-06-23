%Suspension of polarizable particles under rotating external field
%This script solves EOM and performs tangent dynamics

%________________________________________________________________________________________________________________%

clear all ;
close all ;
clc ;

%_______________________________________________________________________________________________________________%
% Inputs

N = 4 ; % number of particles
w = 1.253 ; % Field frequency in Hz
IC = 0 ; % initial phase lag

%________________________________________________________________________________________________________________%
% Constants

r0 = 1 ; % particle radius
a = 2*r0 ; % particle diameter

mu = 1 ; % drag coefficient

F0 = 1 ; % magnetic force scale

A = 1 ; % contact force constant
alpha = 36 ; % contact potential exponent

C1 = F0 / mu * a^4 ; % constant # 1
C2 = A * F0 / mu * a ^ alpha ; % constant # 2

KB = 1 ; % Boltzmann constant
T = 2e-4 ; % Temperature

w_c = 2 * F0 / mu / a ; % critical frequency of a pair of particles

TimeStep = 1E-2 ; %time step

Ma = w / w_c ; % Mason number
Maxtime = 500 * 2 * pi / Ma ; %simulation time
t_w = 10 * 2 * pi / Ma ; % size of the window used for time averaging of magneto-static energy
t_w_p = 0.2 * 2* pi / Ma ; % size of the window used for dumping positions

%________________________________________________________________________________________________________________%
% Initialization

P = zeros(2*N+1,1) ; % The (2N+1)th component refers to field angle

for i = 1 : N
    P(2*i-1,1) = (i-1)*a ;
    P(2*i) = 0 ;
end

P(2*N+1,1) = IC ;

E = zeros(2*N+1,1) ; % error
E(1) = 1 ;

%_______________________________________________________________________________________________________________%
% Solving the ODEs of motion

error = 0 ; % logarithem of initial separation

P_LLE = [] ; % array of separation

for counter = 1 : Maxtime * ( 1 / TimeStep ) + 1
    
    t = (counter - 1) * TimeStep ;
    
    %*****************************************************************************************************************%
    % Calculating dP / dt where dP / dt = f(P) and average MS over tw
    
    P_dot = zeros(2*N+1,1) ;
    
    P_dot(2*N+1,1) = w ;
    
    energy_ins = 0 ;
    
    for i = 1 : N
        for j = 1: N
            if j ~= i
                r = sqrt ( (P(2*j-1)-P(2*i-1))^2 + (P(2*j)-P(2*i))^2 ) ;
                phi = P(2*N+1) - atan2( P(2*j)-P(2*i) , P(2*j-1)-P(2*i-1) ) ;
                del_x = P(2*j-1) - P(2*i-1) ;
                del_y = P(2*j) - P(2*i) ;
                F_r = 3 * ( cos(phi) )^2 -1 ;
                F_t = sin(2*phi) ;
                P_dot(2*i-1) = P_dot(2*i-1) + (C1 / r^5) * ( del_x * F_r + del_y * F_t ) - C2 * del_x / r^(alpha+1) ;
                P_dot(2*i) = P_dot(2*i) + (C1 / r^5) * ( del_y * F_r - del_x * F_t ) - C2 * del_y / r^(alpha+1) ;
                energy_ins = energy_ins + (2 * r0^3 / r^3) * ( 1 - 3 * (cos(phi))^2 ) ;
            end
        end
        P_dot(2*i-1) = P_dot(2*i-1) + randn * (2*mu*KB*T / (TimeStep))^(1/2) / mu ;
        P_dot(2*i) = P_dot(2*i) + randn * (2*mu*KB*T / (TimeStep))^(1/2) / mu ;
    end
    %*****************************************************************************************************************%
    % Calculating dE / dt where dE / dt = J * E , J is Jacobian of EOM and E is separation of two neighbor trajectories
    
    J = zeros(2*N+1,2*N+1) ;
    
    for i = 1 : N
        for j = 1 : N
            if j == i % diagonal and (2N+1) blocks of Jacobian
                for k = 1 : N
                    if k ~= j
                        
                        r = sqrt ( (P(2*k-1)-P(2*i-1))^2 + (P(2*k)-P(2*i))^2 ) ;
                        psi = P(2*N+1) ;
                        theta = atan2( P(2*k)-P(2*i) , P(2*k-1)-P(2*i-1) ) ;
                        phi = psi - theta ;
                        del_x = P(2*k-1) - P(2*i-1) ;
                        del_y = P(2*k) - P(2*i) ;
                        
                        J(2*i-1,2*i-1) = J(2*i-1,2*i-1) + ( F0/2/mu/r^7 ) * ( -2*a^alpha*A * (alpha*del_x^2-del_y^2) * r^(4-alpha) + ...
                            a^4 * (-2*del_x-del_y)*(-2*del_x+del_y) + a^4 * ( cos(2*theta) * ((12*del_x^2-7*del_y^2)*cos(2*psi)+16*del_x*del_y*sin(2*psi)) + ...
                            ( -16*del_x*del_y*cos(2*psi) + (12*del_x^2-7*del_y^2)*sin(2*psi) ) * sin(2*theta) ) ) ;
                        
                        J(2*i-1,2*i) = J(2*i-1,2*i) + ( F0/2/mu/r^7 ) * ( del_x*del_y * ( 5*a^4-2*a^alpha*A*(1+alpha)*r^(4-alpha) ) + ...
                            a^4 * ( 19*del_x*del_y*cos(2*phi) - 8*(-del_x-del_y)*(-del_x+del_y)*sin(2*phi) ) ) ;
                        
                        J(2*i,2*i) = J(2*i,2*i) + ( F0/mu/r^7 ) * ( a^alpha*A * (del_x^2-alpha*del_y^2) * r^(4-alpha) - ...
                            a^4/2 * (-del_x-2*del_y)*(-del_x+2*del_y) + a^4/2 * ( cos(2*theta) * ((12*del_y^2-7*del_x^2)*cos(2*psi)-16*del_x*del_y*sin(2*psi)) + ...
                            ( 16*del_x*del_y*cos(2*psi) + (12*del_y^2-7*del_x^2)*sin(2*psi) ) * sin(2*theta) ) ) ;
                        
                        J(2*i,2*i-1) = J(2*i,2*i-1) + ( F0/2/mu/r^7 ) * ( del_x*del_y * ( 5*a^4-2*a^alpha*A*(1+alpha)*r^(4-alpha) ) + ...
                            a^4 * ( 19*del_x*del_y*cos(2*phi) - 8*(-del_x-del_y)*(-del_x+del_y)*sin(2*phi) ) ) ;
                        
                        J(2*i-1,2*N+1) = J(2*i-1,2*N+1) + ( F0*a^4/mu/r^5 ) * ( 2*del_y*cos(2*phi) - 3*del_x*sin(2*phi) ) ;
                        
                        J(2*i,2*N+1) = J(2*i,2*N+1) + ( F0*a^4/mu/r^5 ) * ( -2*del_x*cos(2*phi) - 3*del_y*sin(2*phi) ) ;
                        
                    end
                end
            else % off-diagonal blocks of Jacobian
                
                r = sqrt ( (P(2*j-1)-P(2*i-1))^2 + (P(2*j)-P(2*i))^2 ) ;
                psi = P(2*N+1) ;
                theta = atan2( P(2*j)-P(2*i) , P(2*j-1)-P(2*i-1) ) ;
                phi = psi - theta ;
                del_x = P(2*j-1) - P(2*i-1) ;
                del_y = P(2*j) - P(2*i) ;
                
                J(2*i-1,2*j-1) = J(2*i-1,2*j-1) + ( F0/mu/r^7 ) * ( a^alpha*A * (alpha*del_x^2-del_y^2) * r^(4-alpha) + ...
                    a^4 * (-2*del_x-del_y)*(-2*del_x+del_y) + a^4 * ( 3*(-4*del_x^2+del_y^2) * (cos(phi))^2 - 2*del_y * ...
                    ( -del_y*cos(2*phi) + 4*del_x*sin(2*phi) ) ) ) ;
                
                J(2*i-1,2*j) = J(2*i-1,2*j) + ( F0/mu/r^7 ) * ( del_x*del_y * ( 5*a^4+a^alpha*A*(1+alpha)*r^(4-alpha) ) + ...
                    a^4 * ( -del_x*del_y/2 * (15+19*cos(2*phi)) + 4*(-del_x-del_y)*(-del_x+del_y)*sin(2*phi) ) ) ;
                
                J(2*i,2*j-1) = J(2*i,2*j-1) + ( F0/mu/r^7 ) * ( del_x*del_y * ( 5*a^4+a^alpha*A*(1+alpha)*r^(4-alpha) ) + ...
                    a^4 * ( -del_x*del_y/2 * (15+19*cos(2*phi)) + 4*(-del_x-del_y)*(-del_x+del_y)*sin(2*phi) ) ) ;
                
                J(2*i,2*j) = J(2*i,2*j) + ( F0/2/mu/r^7 ) * ( -2*a^alpha*A * (del_x^2-alpha*del_y^2) * r^(4-alpha) + ...
                    a^4 * (-del_x-2*del_y)*(-del_x+2*del_y) + a^4 * ( cos(2*theta) * ( (7*del_x^2-12*del_y^2)*cos(2*psi)+16*del_x*del_y*sin(2*psi) ) + ...
                    ( -16*del_x*del_y*cos(2*psi) + (-12*del_y^2+7*del_x^2)*sin(2*psi) ) * sin(2*theta) ) ) ;
                
            end
        end
    end
    
    % dE / dt computation
    
    E_dot = J * E ;
    %*****************************************************************************************************************%
    % Updating P (positions)
    
    for i = 1 : N
        P(2*i-1) = P(2*i-1) + P_dot(2*i-1) * TimeStep ;
        P(2*i) = P(2*i) + P_dot(2*i) * TimeStep ;
    end
    P(2*N+1) = P(2*N+1) + P_dot(2*N+1) * TimeStep ;
    
    %*****************************************************************************************************************%
    % Updating E (error i.e. divergence of two neighbor trajectories)
    
    for i = 1 : N
        E(2*i-1) = E(2*i-1) + E_dot(2*i-1) * TimeStep ;
        E(2*i) = E(2*i) + E_dot(2*i) * TimeStep ;
    end
    E(2*N+1) = E(2*N+1) + E_dot(2*N+1) * TimeStep ;
    
    %*****************************************************************************************************************%
    %Re-normalizing E
    
    E_mag = norm(E) ;
    E = E / E_mag ;
    
    %*****************************************************************************************************************%
    %Keeping track of the stretch/average energy/symmetry order
    %parameter
    
    P_LLE(counter) = error ;
    error = error + log(E_mag) ;  
end
%________________________________________________________________________________________________________________%
%Plot of time-averaged energy/symmetry order parameter/ log(e) vs. time

time = 0 : Maxtime/2/pi*Ma/(size(P_LLE,2)-1) : Maxtime/2/pi*Ma ;

figure
p = plot(time,P_LLE,'b') ;
xlabel('t')
ylabel('log(\epsilon)')
axis ([0 floor(t/2/pi*Ma)+1 -1 1.1*max(P_LLE)])
str = sprintf(strcat('Ma=',num2str(Ma,'%.4f'),',','IC=',num2str(IC,'%.4f'),' - log(separation)'));
title(str);
saveas(p,[str,'.png'])