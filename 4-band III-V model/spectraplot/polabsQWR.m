%clear all;
%load InfinitelongVQWR_16nm_Diameter_v1.mat;
%Q = exciton_pert(Q);
[mlh, Elh] = absQWR_II(Q.CB.WF, Q.VB.WF, [1 1 1], Q.CB.E, Q.VB.E);
[mhh, Ehh] = absQWR_II(Q.CB.WF, Q.VB.WF, [1 1 0], Q.CB.E, Q.VB.E);
figure; stem(1519.2+Elh, mlh, 'bo');
hold on; stem(1519.2+Ehh, mhh, 'ro');