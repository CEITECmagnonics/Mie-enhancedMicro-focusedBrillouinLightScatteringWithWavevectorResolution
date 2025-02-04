% Created on Mon Oct 26
% This script calculates density of states for BLS signal calculations
% If you do not now how to call Python from Matlab start with 
% https://www.mathworks.com/help/matlab/matlab_external/create-object-from-python-class.html
% @author: Ondrej Wojewoda

function [ff, f00fkx] = SpinWaveGreen(kx, Bext, d, N, Nf, n, mu)
pathToSWT = fileparts(which('SpinWaveToolkit.py'));
if count(py.sys.path,pathToSWT) == 0
    insert(py.sys.path,int32(0),pathToSWT);
end
py.importlib.import_module('SpinWaveToolkit')


kx = kx*sqrt(2);
kxi = py.numpy.linspace(min(kx), max(kx), int32(length(kx)));

% Dispersion of two first thickness modes of 30 nm NiFe thin film
material = py.SpinWaveToolkit.Material(697756.87844, 8.0822e-12, 70e-4, 0, 30.43778*2*pi*1e9, 0); % evaporated NiFe, fitted from exp data in Fig. 4
theta = pi/2;
phii = linspace(0, pi, 100);
weff = 3e-6;
boundaryCond = 1;

f00P = zeros(length(phii), length(kx));
lifetimeP = zeros(length(phii), length(kx));

i=0;
for phi=phii
    i=i+1;
    NiFeChar = py.SpinWaveToolkit.DispersionCharacteristic(Bext, material, d, kxi, theta, phi, weff, boundaryCond);
    f00P(i,:) =  double(NiFeChar.GetDispersion(n))*1e-9/(2*pi);
    lifetimeP(i,:) =  double(NiFeChar.GetLifetime(n))*1e9;
end
f00P = squeeze(f00P);
[RHO, PHI] = meshgrid(kx, phii);
[X, Y] = pol2cart(PHI, RHO);
xi = linspace(-max(kx)/sqrt(2), max(kx)/sqrt(2), N);
yi = linspace(-max(kx)/sqrt(2), max(kx)/sqrt(2), N);
[Xi, Yi] = ndgrid(xi, yi);

Ff00 = scatteredInterpolant(X(:), Y(:), f00P(:), 'linear', 'none');
f00 = Ff00(Xi, Yi);
Flifetime = scatteredInterpolant(X(:), Y(:), lifetimeP(:), 'linear', 'none');
lifetime = Flifetime(Xi, Yi);

ff = linspace(min(f00(:))-2, max(f00(:))+2, Nf);
f00fkx = zeros(length(ff),N,N);
i=0;

for f=ff
    i=i+1;
    % Bose-Einstein distribution
    hbar = 1.0545718e-34;
    kb = 1.38064852e-23;
    T = 300;
    mu=mu;
    B=1;
    BE = 1./(B*exp((hbar*2*pi*(abs(f*1e9)-mu))/(kb*T))-1);
    f00fkx(i,:,:) = sqrt(BE)./(4*pi^2*abs(f00-f).^2+(2./lifetime).^2);        
end
end

