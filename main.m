clc;
clear;
% loop = 2000;
% n_sigma = [0.001 0.1];
% N = [20 50 100];
% 
% % Noiseless case
% for i = 1:3
%     NPT_0(N(i), loop);
% end
% 
% % Noisy case s is known
% for i = 1:3
%     NPT_1(N(i), loop, n_sigma(1));
%     NPT_1(N(i), loop, n_sigma(2));
% end
% 
% % Noisy case s is not known
% for i = 1:3
%     NPT_2(N(i), loop, n_sigma(1));
%     NPT_2(N(i), loop, n_sigma(2));
% end

fig = imread('cheetah.bmp');

rec(fig,50,0.001);







