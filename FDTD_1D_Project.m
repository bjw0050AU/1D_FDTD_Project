% This code implements a 1-D scalar wave equation following the
% finite-difference time domain (FDTD) method of Taflove and Hagness

% Physical constants and user-defined values
mu_0 = pi*4e-7; %permeability of free space in H/m
eps_0 = 8.854e-12; %permittivity of free space in F/m
eps_r = ; %relative permittivity of the material
c = 1/sqrt(mu_0*eps_0); %speed of light in free space in m/s
f = ; %frequency of the input wave in Hz
lambda = c/f; %free space wavelength in m
% Set up the space and time dimensions
dx = ; %user-defined spatial step in m
xmax = ; %max value of distance (x) in the simulation space
x = 0:deltax:xmax; %all values of distance (x) in the simulation space
tmax = ; %simulation stops after time t=tmax (in seconds)
S = ; %user-defined Courant stability factor (Can't be 1, but can be less)
deltat = S*deltax/c; %time step in seconds
t=0; %start time in seconds
% Put dielectric material into the simulation space
eps = ones(length(x),1).*eps_0; %initialize permittivity everywhere
s1 = ; %location index of the material boundary
eps(s1:end) = eps(s1:end)*eps_r;%add dielectric material past boundary


% Initialize electric field
E = zeros(length(x),3); %E=0 everywhere and for all previous time

% Set locations for the virtual electric field probes
x0 = 1; %index of field probe E0
x1 = ; %index of field probe E1
t1 = ; %time of the snapshot in seconds

%create a blank figure for the FDTD animation
h=figure;

% Begin the FDTD update loop
while t<tmax %update until the max time value is reached
    %implement the 1-D scalar update equation
    E(2:end-1,3) = ;
    %E-field of the incoming wave (turns off after 5 cycles)
    if t<=5/f
        E(1,3) = sin(2*pi*f*t);
    elseif t>5/f
        E(1,3) = 0;
    end
    %update the plot animation every 5 steps
    if mod(round(t/deltat),5)==0
        figure(h)
        plot(x,E(:,3))
        xline(x(s1))
        ylim([-3 3])
        title([sprintf('t=%f',t*1e6) '\mus'])
        drawnow
    end
    %take a snapshot at the user-specified time t=t1
    if (t>t1-deltat)&&(t<t1+deltat)
        figure
        plot(x,E(:,3))
        xline(x(s1))
        ylim([-3 3])
        title([sprintf('t=%f',t*1e6) '\mus'])
    end
    %store values for the virtual field probes
    E0(round(t/deltat) + 1) = E(x0,3);
    E1(round(t/deltat) + 1) = E(x1,3);
    %move forward one time step
    t = t+deltat;
    E(:,1) = E(:,2);
    E(:,2) = E(:,3);
end
%add code here to produce figures after simulation ends
