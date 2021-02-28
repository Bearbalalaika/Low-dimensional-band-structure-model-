%clear all;
%load triple_3p5_3p5_3p5_0p20_0p15_0p20_modified_IV.mat;
Q = exciton_pert(Q);
out = absII_cont_pol(Q);


%axis([1.580 1.630 0 1.1]); 
set(gca,'fontsize',28,'linewidth',3);

%h = get(gca, 'Children');
%set(h, 'LineWidth', 3);



box on;
set(gca, 'FontName', 'Times');

%%
