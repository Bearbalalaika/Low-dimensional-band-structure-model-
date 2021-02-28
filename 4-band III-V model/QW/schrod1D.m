%(Hermitian 031229)
%[E, WF] = schrod1D(x, V, mass, levels)
%x = coordinate [nm]
%V = potential [eV]
%mass = mass (m0)
%levels = number of levels
%
%E = energy (eV)
%WF = wavefunctions

function [E, WF] = schrod1D(z, V, mass, levels)

LZ = length(z);

%Adjust the size of mass
if length(mass) < LZ
   mass = ones(1, LZ)*mass;
end

%Constants and unit converision
physconst;
V = V'*eV;
dz = mean(diff(z))*1E-9;
mass = me*mass';

%Produce the Hamiltonian operator
M = finitDiff1D(1./mass, zeros(size(V)), V, -hbar^2/2, 0, 1, dz);

%Find the egienvalues and eigenvectors
OPTIONS.disp = 0;
[WF0, E] = eigs(M, levels, 'SA', OPTIONS);

%Sort the energies and corresponding wavefunctions
[E, SI]= sort(diag(E)/eV);
for i = 1:length(E)
   WF(:,i) =  WF0(:,SI(i));
end