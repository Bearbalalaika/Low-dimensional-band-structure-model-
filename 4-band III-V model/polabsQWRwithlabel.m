%clear all;
%load InfinitelongVQWR_30nm_Diameter_v8.mat;
%Q = exciton_pert(Q);
[mlh, Elh] = absQD_II(Q.CB.WF, Q.VB.WF, [0 0 1], Q.CB.E, Q.VB.E);
[mhh, Ehh] = absQD_II(Q.CB.WF, Q.VB.WF, [1 1 0], Q.CB.E, Q.VB.E);
figure; 

Transition = zeros(100,5);

for i=1:10
    for j=1:2:10
        stem(1519.2+Elh(i,j), mlh(i,j), 'b');
        text(1519.2+Elh(i,j), mlh(i,j)+0.01, ['\color{blue}' num2str([i j])]);
        %text(1.5192+Elh(i,j)/1000, mlh(i,j)+0.15, ['\color{black}' num2str(1.5192+Elh(i,j)/1000) 'eV']);
        hold on;
        stem(1519.2+Ehh(i,j), mhh(i,j), 'r');

    end
end

%k = 1
%for i=1:10
%    for j=1:2:10
%        Transition(k,1) = 1.5192+Elh(i,j)/1000; %transition energy
%        Transition(k,2) = mlh(i,j);             %light holes
%        Transition(k,3) = 1.5192+Ehh(i,j)/1000; %transition energy
%        Transition(k,4) = mhh(i,j);             %heavy holes
%        k=k+1;
%    end
%end  

set(gca,'fontsize',28,'linewidth',3);


%save heksVQWR_28nm-15nm_003Al Transition -ascii -tabs    