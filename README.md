# Motivation_in_WM_Analysis

This matlab code will run the analyses and make the figures from the Motivation in Working Memory paper (Grogan et al.).

It requires:
Matlab (preferably R2018b although other versions should work)
matlib (github.com/sgmanohar/matlib)
MemToolbox (www.memtoolbox.org)
barwitherr.m
regression_line_ci.m


The data is available at https://osf.io/tfnzw/ (and should be saved ../../Experiments/Results/ArrowMotivCue, or change the dataFolder variable in the scripts below), as are the experimental scripts that created them.

The entire pipeline can be run by ArrowAnalysisPipeline.m

The repo contains:

ArrowAnalysisPipeline.m % run entire pipeline
ArrowDataLoad.m % load up data files
ArrowDemographics.m % make table of demographics
ArrowFitTrajInterp.m % modelling analysis for response trajectories
ArrowModelCall.m % apply mixture models to data
ArrowModelFitsTrajNormalise.m % fit model to response trajectories at each time point
ArrowModelLoad.m % load up data for mixture modelling
ArrowMotivAnalysis.m % expt 1 analysis
ArrowMotivCueBlockedAnalysis.m % expt 2 analysis
ArrowMotivCueBlockedModelling.m % expt 2 modelling 
ArrowMotivModelAnalysis.m % expt 1 modelling
ArrowMotivMotorAnalysis.m % expt 5a+b analysis
ArrowMotivMotorCorrels.m % correlate motor + WM tasks
ArrowMotivMotorModelling.m % expt 5a+b modelling 
ArrowMotivReallocateAnalysis.m % expt 4 analysis
ArrowMotivReallocateModelling.m  % expt 4 modelling 
ArrowMotivRetrocueAnalysis.m % expt 3 analysis
ArrowMotivRetrocueModelling.m % expt 3 modelling
ArrowsAllFigs.m % make all figures in paper
ArrowTrajAll.m % Analyse all trajectories
ArrowTrajAnalysis.m % analyse and format trajectories
ArrowTrajLoad.m % load up trajectory data 
col.m % make matrix into column
emptyLegend.m % plot a legend on empty data
JoinUpErrorPoints.m % join up points with lines
MixMax.m % get min and max in one function
p2stars.m % make p values into stars
PlotTrajPrecVsTimeInterp.m % plot arrow trajectories for one expt
SuperTitle.m % one main title on a subplot

Model\ArrowRaceAnalysis.m % analyse one race model fit for one expt
Model\ArrowRaceAnalysisAll.m % analyse race model fits for all expts
Model\ArrowRaceModelFitting.m % fit race model to one expt
Model\ArrowRaceModelFittingAll.m % fit race model to all expts
Model\get_neg_log_like.m % calculate negative log likelihood
Model\raceModelFit.m % fit race model to one person/condition
Model\raceModelWrapper.m % actual race model
Model\sim_race.m % simulate race/LBA model
