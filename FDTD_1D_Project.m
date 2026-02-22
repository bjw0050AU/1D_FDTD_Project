% This code implements a 1-D scalar wave equation following the
% finite-difference time domain (FDTD) method of Taflove and Hagness
% Only finding Electric Field!

clear;      %clear workspace
clc;        %clear console
close all;

%To-Do: Improve accuracy and values used in simulation
% Physical constants and user-defined values
mu_0 = pi*4e-7; %permeability of free space in H/m
eps_0 = 8.854e-12; %permittivity of free space in F/m
eps_r = 10; %relative permittivity of the material
c = 1/sqrt(mu_0*eps_0); %speed of light in free space in m/s
f = 20e6; %frequency of the input wave in Hz (20 MHz)
lambda = c/f; %free space wavelength in m

% Set up the space and time dimensions
dx = 0.2; 
%user-defined spatial step in m 
% (make smaller as it will force wave to propagate further)

x = 0:dx:(dx*1000); 
xmax = 1000*dx;
%all values of distance (x) in the simulation space
%max value of distance (x) in the simulation space defined as dx*10000


S = 0.99; %user-defined Courant stability factor (Can't be 1! Can be less.)
dt = S*dx/c; %time step in seconds
fs = 1/dt;
t=0; %start time in seconds
tmax = dt*(length(x)-1); %time it takes for all points to compute spatially
n = tmax/dt; %find number of timesteps


% Put dielectric material into the simulation space
eps = ones(length(x),1).*eps_0; %initialize permittivity everywhere
s1 = 500; %location index of the material boundary
eps(s1:s1+15) = eps(s1:s1+15)*eps_r;%add dielectric material past boundary


% Initialize electric field
E = zeros(length(x),3); %E=0 everywhere and for all previous time

% Set locations for the virtual electric field probe
x1 = s1+4; %index of field probe E1 (material)
t1 = 6*(1/f); 
iteration = 1;
amp_meas = 0; %variable for E0 amp measured
%time of the snapshot in seconds (time for full-wave to complete)


%create a blank figure for the FDTD animation
h=figure;


% Begin the FDTD update loop
while t<tmax %update until the max time value is reached
    %indexing is (space step, time step)
    %implement the 1-D scalar update equation
    E(2:end-1,3) = (1/(mu_0*eps(iteration)))*((dt/dx)^2)*((E(3:end,2)-2*(E(2:end-1,2))...
        +E(1:end-2,2))) + 2*E(2:end-1,2) - E(2:end-1,1);
    %E-field of the incoming wave (turns off after 5 cycles)
    if t<=5/f
        E(1,3) = sin(2*pi*f*t);
        E(end,3) = 0;
    elseif t>5/f
        E(1,3) = 0;
        E(end,3) = 0;
    end
    if (t>5/f) & (amp_meas == 0) 
        %figure;
        E0_magn = max(E(:,3)); %find initial field Amp.
        %plot(abs(fft(E)/length(E))); 
        amp_meas = 1;
        m = size(E(s1:end));
        %plot normalized fft (for input sinusoid)
    end
    %update the plot animation every 5 steps
    if mod(round(t/dt),5)==0
        figure(h)
        plot(x,E(:,3))
        xline(x(s1))
        xline(x(s1+15))
        ylim([-3 3])
        title([sprintf('t=%f',t*1e6) '\mus'])
        drawnow
    end
    %take a snapshot at the user-specified time t=t1
    if (t>t1-dt)&&(t<t1+dt)
        figure
        plot(x,E(:,3))
        xline(x(s1))
        ylim([-3 3])
        title([sprintf('t=%f',t*1e6) '\mus'])
    end
    %store values for the virtual field probes
    if (iteration > 500)
        Et(round(t/dt) + 1 - 498) = E(x1,3); %map transmitted field (reaches medium)
    end
    if (iteration > 550) & (amp_meas ~= 2) & (amp_meas ~= 3)
        Et_magn = max(Et); %find initial field Amp.
        amp_meas = 2;
    elseif (iteration > 950) & (amp_meas ~= 3)
        Er_magn = max(E(105:255,3)); %find initial field Amp.
        amp_meas = 3;
    end
        %zeroing condition to prevent superposition for 
        % reflect & initial being the transmitted value
    %move forward one time step
    t = t+dt;
    E(:,1) = E(:,2);
    E(:,2) = E(:,3);
    iteration = iteration + 1;
end

%add code here to produce figures after simulation ends
%Et = reflected+initial into medium
Tran_Coeff = (E0_magn-Er_magn)/E0_magn;
eta_0 = 120*pi;
eta_r = (Tran_Coeff / (2-Tran_Coeff)); %if negative acts capacitive

Et_filtered = abs(Et(23:end)); %filter zeros in Et matrix (pre-saved)
indices = find(Et_filtered < 0.1); %find zeros

wavelength_mat = (indices(3) - indices(1))*dx; %find material lambda

Up_mat = f*wavelength_mat; %find material propagation velocity






