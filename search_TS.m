clear;
load('simulation_data.mat');
H = H(:,1:990);
H(sum(H,2)==0,:) = [];
H = cycleRemoval(H);

addpath('C:\Users\mpacenti\Documents\Banihashemi-Pacenti.m');

g = 6; %girth

%parameters and tables loading
K=5;
% TOTEX = {[2,3]	[3,3;4,3;4,4]};
TOTEX = {[2,3,4]	[4,4;5,4;5,5]};
amax = 10;
bmax = 5;

[MAT_CONT, I] = newTSenumerator_irreg23(K, g, amax, bmax, H, TOTEX);