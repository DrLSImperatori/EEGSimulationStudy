%% Main Code 
clear all;
close all;
clc;
 
restoredefaultpath
addpath('/home/laurai/matlab/fieldtrip-20190219/')
ft_defaults


addpath(genpath('/home/laurai/matlab/BBCBCode/'))
addpath(genpath('/home/laurai/matlab/EEGSimulationStudy/'))
addpath('/home/laurai/matlab/CST/')

tic
generate_datasets_RIPLMFG_shuffle(100, 20); % generate_datasets_RIPLMFG_shuffle(100, 20); generate_datasets_LIPLMFG_AAFT(100, 20); generate_datasets_RIPLMFG_AAFT(100, 20);
toc


