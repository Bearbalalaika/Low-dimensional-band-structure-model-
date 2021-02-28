function N = Occupation(ncb, nvb,[Q.CB.E],[Q.VB.E]) % occupation unmber taken from experimental data     
T=10.6*0.0862; % exciton effective temperature [meV]
u=1.549; % effective fermi energy [eV]
deltaQ=41; % differentce of Gs energy position [meV] experiment-theory
Ebg=1.519; % band gap energy [eV]
N=1/(exp(((Q.CB.E(ncb)+Q.VB.E(nvb)+deltaQ)/1000+Ebg-u)/(T*0.001))+1);