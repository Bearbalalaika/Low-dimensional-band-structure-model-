%clear all;
%load UNEQUAL_QD1_9nm_QD2_0d5nm_DQD_4dz_Barrier_4nm.mat;
Q = exciton_pert(Q);
out = absII_cont_pol(Q);
axis([1.585 1.62 -1 1]); box on;    
