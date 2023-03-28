
clear all
close all

%addpath('new17Optica')

load VELM_Call.mat

mode.verbVELM=1; % 0: non-verbose; 1: just field plots; 2: verbose mode
dox=4;
ivett=0;

colordef black
[velm] = f_CallVELM(mesh,mode,dox,ivett,fil_str);