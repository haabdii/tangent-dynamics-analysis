% System under study:
% Suspension of polarizable particles under rotating magnetic field
% This script solves EOM and performs tangent dynamics

%________________________________________________________________________________________________________________%

clear all ;
close all ;
clc ;

%_______________________________________________________________________________________________________________%
% Inputs

N = 4 ; % number of particles
w = 1.22 ; % Field frequency
IC = 0 ; % initial phase lag between the field orienation and the chain orientation

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

w_c = 2 * F0 / mu / a ; % critical frequency for a pair of particles
Ma = w / w_c ; % Mason number

TimeStep = 1E-2 ; %time step
Maxtime = 500 * 2 * pi / Ma ; %simulation time
%t_w = 10 * 2 * pi / Ma ; % size of the window used for time averaging of magneto-static energy
%t_w_p = 0.2 * 2* pi / Ma ; % size of the window used for dumping positions

%________________________________________________________________________________________________________________%
% Initialization

P = init_position(N, a, IC) ; % initializing the particle positions
E = init_E(N) ; % initializing the separation vector
%_______________________________________________________________________________________________________________%
% Solving the ODEs of motion

error = 0 ; % logarithem of initial separation

P_LLE = [] ; % array of separation

for counter = 1 : Maxtime * ( 1 / TimeStep ) + 1
    
    t = (counter - 1) * TimeStep ;
    
    %*****************************************************************************************************************%
    % Calculating particle velocities, dP/dt = -Force/drag_coefficient
    
    P_dot = velocity_calculator (N, w, P, C1, C2, alpha) ;
    %*****************************************************************************************************************%
    % Calculating dE / dt where dE / dt = J * E , J is Jacobian of EOM and E is separation of two neighbor trajectories
    
    J = Jacobian_calculator(N, P, F0, mu, a, alpha, A) ;
    %*****************************************************************************************************************%
    % dE / dt computation
    
    E_dot = J * E ;
    %*****************************************************************************************************************%
    % Updating P (positions)
    
    P = update_position (N, TimeStep, P, P_dot) ;
    %*****************************************************************************************************************%
    % Updating E (separation i.e. divergence of two neighbor trajectories)
    
    E = separation_update(N, E, E_dot, TimeStep) ;
    %*****************************************************************************************************************%
    %Re-normalizing E
    
    E_mag = norm(E) ;
    E = E / E_mag ;
    %*****************************************************************************************************************%
    %Keeping track of the stretch
    
    P_LLE(counter) = error ;
    error = error + log(E_mag) ;
end
%________________________________________________________________________________________________________________%
%Plot of log(e) vs. time

time = 0 : Maxtime/2/pi*Ma/(size(P_LLE,2)-1) : Maxtime/2/pi*Ma ;

figure
p = plot(time,P_LLE,'b') ;
xlabel('t')
ylabel('log(\epsilon)')
axis ([0 floor(t/2/pi*Ma)+1 -1 1.1*max(P_LLE)])
str = sprintf(strcat('Ma=',num2str(Ma,'%.4f'),',','IC=',num2str(IC,'%.4f'),' - log(separation)'));
title(str);
saveas(p,[str,'.png'])
