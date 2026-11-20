clear
close all
clc

% load results
load('modResults.mat','gf');    % gf (folds × metrics × fileLoop)
load('allFits_cells.mat'); % allData (for grid size)

datSize = size(allData);   % [nNoise × nRep]
nNoise  = datSize(1);
nRep    = datSize(2);

% pick fit metric (logLik or adjR²)
col_logLik = 3; 
metric     = col_logLik;

% average across folds
G = squeeze(gf(:,metric,:));   % [folds × fileLoop]
fit_by_file = mean(G,1);       % [1 × (nNoise*nRep)]

% reshape to [nNoise × nRep], then average over reps
LL = reshape(fit_by_file, datSize);   % [nNoise × nRep]
LLmean = mean(LL,2);                  % [nNoise × 1]

% reshape into ori × pos grid
nOri = 10; nPos = 10;                 % adjust if not 10×10
LLgrid = reshape(LLmean,[nOri nPos]);

% heatmap of fit surface
figure;
imagesc(1:nPos,1:nOri,LLgrid);
axis xy; colorbar;
xlabel('Positional noise index');
ylabel('Orientation noise index');
title('Model fit surface (logLik)');