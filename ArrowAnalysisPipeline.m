% ArrowAnalysisPipeline
% Run entire analysis pipeline for the Motivation in WM paper
% Run this in the /Motivation_in_WM_Analysis folder (or where you
% downloaded the github code)
% 
% 
% See the README for the required packages and scripts, and make sure they
% are all on your matlab path

%% make folders for data + figures

system('mkdir Data');
system('mkdir Figs');
addpath ./Model % where the LBA modelling scripts are
addpath ./ % current folder

%% experiment 1

ArrowMotivAnalysis
ArrowMotivModelAnalysis

%% experiment 2

ArrowMotivCueBlockedAnalysis
ArrowMotivCueBlockedModelling

%% experiment 3

ArrowMotivRetrocueAnalysis
ArrowMotivRetrocueModelling

%% experiment 4

ArrowMotivReallocateAnalysis
ArrowMotivReallocateModelling


%% experiment 5a

ArrowMotivMotorAnalysis('A')
ArrowMotivMotorModelling('A')

%% experiment 5b

ArrowMotivMotorAnalysis('B')
ArrowMotivMotorModelling('B')


%% across experiment stuff

ArrowMotivMotorCorrels

%% trajectory analyses
ArrowTrajAll
ArrowFitTrajInterp

%% LBA modelling
cd ./Model/
ArrowRaceModelFittingAll
ArrowRaceAnalysisAll

cd ..
%% draw figures
ArrowsAllFigs

%% Make Demographics table

ArrowDemographics