clear all;
close all;

addpath('Functions');
addpath('plotting_function');

Tms = double(imread('facade_oc1_missing.png'));
%Tms = double(imread('facade_oc2_missing.png'));
%Tms = double(imread('peppers_occlusion1.png'));
%Tms = double(imread('peppers_occlusion2.png'));
Qms = (Tms ~= 0);

tau = [32 32 1];

param.incR{1} = [1 2 4 8 16 24 32];
param.incR{2} = [1 2 4 8 16 32 64 96 128 160 192 256];
param.incR{3} = [1 2 4 8 16 24 32];
param.incR{4} = [1 2 4 8 16 32 64 96 128 160 192 256];
param.incR{5} = [1];
param.incR{6} = [1 3];

param.delta = 1e-5; % shreshold of error between observed elements and reconstructed elements

sc = 255;
[X, histo, histR] = MDT_Tucker_incR(Tms/sc,Qms,tau,param);

Xest = sc*X;

figure(1);
imagesc([uint8(Tms) uint8(Xest)]);


