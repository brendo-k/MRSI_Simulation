%GPU tests]
clear all
dev = gpuDevice;
timer = tic;

dI=3.0; %[ppm]
dS=3.4; %[ppm]
JIS=6.5; %[Hz]
B0 = 3;
gamma=42577000;  %[Hz/T]

np=2048; %[points]
dwelltime=500e-6; %[s]
sw=1/dwelltime; %[Hz]
t=([0:dwelltime:dwelltime*(np-1)]); %[s]

I0=([1 0;0 1]);
Ix=(0.5*[0 1;1 0]);
Iy=(1i/2)*([0 -1;1 0]);
Iz=(1/2)*([1 0;0 -1]);

S0=I0;
Sx=Ix;
Sy=Iy;
Sz=Iz;

%Now set up the key two spin basis matrices:
I0Sx=kron(I0,Sx);
I0Sy=kron(I0,Sy);
I0Sz=kron(I0,Sz);

IxS0=kron(Ix,S0);
IyS0=kron(Iy,S0);
IzS0=kron(Iz,S0);

%Now set up the F matrices:
Fx=IxS0+I0Sx;
Fy=IyS0+I0Sy;
Fz=IzS0+I0Sz;

IS=IxS0*I0Sx + IyS0*I0Sy + IzS0*I0Sz;

%Now set up the equilibrium Density matrix:
d0=Fz;

%Set up the hamiltonian operator for the 90 and 180 degree pulses:
H90=Fy*pi/2;

%Set up the free evolution hamiltonian operators for the delay periods:
Hevol_dI=(dI*B0*gamma*2*pi/10e6)*IzS0;
Hevol_dS=(dS*B0*gamma*2*pi/10e6)*I0Sz;
Hevol_JIS=(JIS*2*pi)*IS;
Hevol=Hevol_dI+Hevol_dS+Hevol_JIS;

dtemp=expm(-1i*H90)*d0*expm(1i*H90);                 %90 degree excitation
S = (zeros(100, length(t)));
for j = 1:1
    for n=1:length(t)
        %get FID point:
        S(j,n)=trace((Fx+1i*Fy)*dtemp);
        %evolve by timestep;
        dtemp=expm(-1i*t*Hevol)*dtemp*expm(1i*t*Hevol);
        S(j,n)=trace((Fx+1i*Fy)*dtemp);
    end
end
wait(dev)
toc(timer);
