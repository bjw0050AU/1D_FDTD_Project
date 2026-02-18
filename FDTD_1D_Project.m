% This code implements a 1-D scalar wave equation following the
% finite-difference time domain (FDTD) method of Taflove and Hagness
% Only finding Electric Field!

clear;      %clear workspace
clc;        %clear console

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
t=0; %start time in seconds
tmax = dt*1e5; 
%Simulation stops after time t=tmax (in seconds)--reach Xmax 5 times
%Must at least allow time for wave to completely travel


% Put dielectric material into the simulation space
eps = ones(length(x),1).*eps_0; %initialize permittivity everywhere
s1 = 500; %location index of the material boundary
eps(s1:end) = eps(s1:end)*eps_r;%add dielectric material past boundary


% Initialize electric field
E = zeros(length(x),3); %E=0 everywhere and for all previous time

% Set locations for the virtual electric field probes
x0 = 5; %index of field probe E0 (free space)
x1 = s1; %index of field probe E1 (material)
t1 = (1/f); 
iteration = 1;
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
    %update the plot animation every 5 steps
    if mod(round(t/dt),5)==0
        figure(h)
        plot(x,E(:,3))
        xline(x(s1))
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
    E0(round(t/dt) + 1) = E(x0,3);
    E1(round(t/dt) + 1) = E(x1,3);
    %move forward one time step
    t = t+dt;
    E(:,1) = E(:,2);
    E(:,2) = E(:,3);
    iteration = iteration + 1;
end
%add code here to produce figures after simulation ends
