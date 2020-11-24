%% Cluster-based control
% Please cite: 
% Nair et al., Cluster-based feedback control of turbulent post-stall separated flows, JFM 2019.
% Kaiser et al., Cluster-based reduced-order modelling of a mixing layer, JFM 2014.

clear;clc;close all;

%% Add libraries 

addpath('src/');

%% Load data

load Force_data.mat;                     % CL,CD,CL_dot,t
M          = 0.3;                        % Mach number
Re         = 23000;                      % Reynolds number
width      = 0.2;                        % width in spanwise direction
area       = 1;                          % control surf int (set as necessary)
force_norm = width*M^2;                  % normalizing force
power_norm = width*M^3;                  % normalizing power

%% Additional parameters

K          = 10;                         % number of clusters (for clustering)
gamma      = 0.025;                      % Power tradeoff

%% Normalize variables

x = [(CL(1:end) - min(CL(1:end)))/(max(CL(1:end)) - min(CL(1:end))) ...
    (CL_dot(1:end) - min(CL_dot(1:end)))/(max(CL_dot(1:end)) - min(CL_dot(12:end))) ...
    (CD(1:end) - min(CD(1:end)))/(max(CD(1:end)) - min(CD(1:end)))];
t = tf(1:end);

%% Prepare Data & options

Data2crom.dt = t(2)-t(1);
Data2crom.t  = t;
Data2crom.ts = x;
params_user.nRepetitions         = 10;
params_user.optimalClustering    = 'sparsity';
params_user.ClusterOrdering      ='transitions';
params_user.save                 = 1;
params_user.verbose              = 0;
params_user.plot                 = 0;
params_user.nClusters      	     = K; 

%% Run CROM (Cluster-based reduced order modeling)

CROMobj = CROM(Data2crom,params_user);
CROMobj.run;

%% Extract clustering centroids

ai         = CROMobj.Data.ts;
c1_Centers = CROMobj.c1_Centroids;
c1_Labels  = CROMobj.c1_Labels;

%% Plotting

Cluster_plot(K,CL,CD,CL_dot,c1_Centers,c1_Labels,ai);

%% Cluster-based control

[J,bi] = CBC(K,t,area,force_norm,power_norm,gamma,c1_Centers);
