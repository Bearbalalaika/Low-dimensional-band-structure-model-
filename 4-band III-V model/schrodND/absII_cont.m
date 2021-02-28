function [E, I, En0, Mn0] = absII_cont(Q, e)

E = 0:0.1:800;
I = zeros(size(E));

[Mn0, En0] = absQD_II(Q.CB.WF, Q.VB.WF, e, Q.CB.E, Q.VB.E);
%
Eb = 0;
if isfield(Q, 'X')
    disp('Exciton binding included')
    Eb = Q.X.Ex; 
end
    
En0 = En0-Eb;
%}
%En = real(reshape(En0, 1, size(En0, 1)*size(En0, 2)));
En = abs(real(reshape(En0, 1, size(En0, 1)*size(En0, 2))));
Mn = reshape(Mn0, 1, size(Mn0, 1)*size(Mn0, 2));
%indx = round(En*10)+1;
%I(indx) = I(indx) + abs(Mn);

Mn
En
for i = 1:length(En)
    I(round(En(i)*10)+1) = I(round(En(i)*10)+1) + Mn(i);  
    
end

I = gaussBlur(I, 25);


%figure; plot(E, I, 'k'); hold on; stem(En, Mn/30, 'b');