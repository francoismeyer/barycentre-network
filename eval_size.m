%%______________________________________________________________________________
%%
%% Copyright (C) 2025 Francois G. Meyer <FMeyer@Colorado.Edu>
%%
%% All rights reserved.
%%
%% Redistribution and use in source and binary forms, with or without
%% modification, are permitted provided that the following conditions are met:
%% 
%%   a. Redistributions of source code must retain the above copyright notice,
%%      this list of conditions and the following disclaimer.
%% 
%%   b. Redistributions in binary form must reproduce the above copyright
%%      notice, this list of conditions and the following disclaimer in the
%%      documentation and/or other materials provided with the distribution.
%% 
%%   c. Neither the name of the copyright holders nor the names of any
%%      contributors to this software may be used to endorse or promote products
%%      derived from this software without specific prior written permission.
%% 
%% 
%% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS on an
%% "AS IS" basis. THE COPYRIGHT HOLDERS AND CONTRIBUTORS MAKE NO
%% REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  BY WAY OF EXAMPLE, BUT
%% NOT LIMITATION, THE COPYRIGHT HOLDERS AND CONTRIBUTORS MAKE NO AND
%% DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS FOR
%% ANY PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE
%% ANY THIRD PARTY RIGHTS.
%% 
%% THE COPYRIGHT HOLDERS AND CONTRIBUTORS SHALL NOT BE LIABLE TO LICENSEE OR
%% ANY OTHER USERS OF THIS SOFTWARE FOR ANY INCIDENTAL, SPECIAL, OR
%% CONSEQUENTIAL DAMAGES OR LOSS AS A RESULT OF MODIFYING, DISTRIBUTING, OR
%% OTHERWISE USING THIS SOFTWARE, OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED
%% OF THE POSSIBILITY THEREOF.
%%
%%
%%_Module_Name : eval_size.m
%%
%%_Description : this script evaluates the error between the estimated Frechet mean
%% using the Soules basis, and the edge probability used to generate
%% the SBMs
%%              
%% The core function is barycentreSBM.m. The rest is visualisation to get some plots.
%%              
%% 
%%
%%_References :
%%
%%_Remarks : None
%%
%%_Author :                 Francois G. Meyer
%%
%%_Revisions History: 2025 Initial keying
%%
%%______________________________________________________________________________
%% 

addpath (genpath('./../synt graph'));
addpath (genpath('./../tree search'));

labarre='–\|/';		% for display; ignore

%
% range of sample size = [n0,n1]
%

n0 = 100;
n1 = 1300;                           %% range of graph size

nValues = 10;                         %% how many values of the sample size within the range of sample sizes
                                     %% in plain English: number of ticks on the x-axis

T = 1;				     % number of samples
equalSize   = 0;
equalChance = 0;
M = 4;

maxTrials = 64;                      %% number of experiments for a fixed sample size

logq = log(n1/n0)/nValues;           %% we use a log scale for the sample size; this is the increment on a log scale

smplSz  = zeros (1, nValues);
dl1   = zeros (maxTrials, 1);        %% l1 distances for a trial, for a fixed sample size
histl1 = zeros (maxTrials, nValues); %% l1 distances for all trials, for a fixed sample size
meanl1 = zeros (1, nValues);         %% mean l1 distances (computed over maxTrials) for all sample sizes

fprintf ('\n');
lastsize = 0;

for nF=1:nValues		%%  graph size in float
    
    Nf = n0*exp((nF-1)*logq);

    n = round (Nf);		%% graph size
    
    %% we run maxTrials independent experiments for the sample size N
    
    ix = 1;  
    for ntrials=1:maxTrials

      %%  compute the difference between the barycentre and the population mean, P.

      dl1 (ntrials) =  barycentreSBM (n, M, equalSize, equalChance, T);
      
      %% next two lines for display
      fprintf(1,'\b%c', labarre(ix));
      ix = mod(ix,4)+1;
    end
    
    histl1 (:,nF) = dl1;
    meanl1 (:,nF) = mean (dl1);
    
    % save the sample size
    
    smplSz (nF) = Nf;
    
    fprintf(repmat('\b', 1, lastsize));
    lastsize = fprintf('\n%2d\n\n',nF);
end

% = hamming distance (population graph median, sample graph median)

coeff3   = polyfit (log(smplSz),log (meanl1),1)
linereg3 = polyval(coeff3,log(smplSz))';

%
% figure 2: hamming distance (population graph median, sample graph median)
%

save ("evalsize.mat","coeff3", "linereg3", "smplSz", "meanl1", "histl1");

nValues = 10;
f = figure;f.Position = [2000 1000 1000 1000]; hold on;
lesmerde={'100', '129', '167', '216', '279', '361', '466', '602', '778', '1006'};
p1 = boxplot(log(abs(histl1)), smplSz,'Labels', lesmerde, 'colors',[1 0 0]);
set(p1,'LineWidth',4);
p2      = plot ([1:nValues],linereg3,'blue');
set(p2,'LineWidth',4);
set (gca,'FontSize', 36);
set (gca,'FontName','Times New Roman');
set(get(gca, 'XLabel'),'String','n','FontName','Times New Roman','FontSize',48, 'LineWidth', 4);
set(get(gca, 'YLabel'),...
    'String','$\log_{10}\big(n^{-2}\|{\bf P}-\widehat{{\mu\mkern -10.2mu\mu}}_N[{\mathrm I}\mkern -4mu{\mathrm P}]\|_1\big)$',...
    'Interpreter','latex', 'FontName','Times','FontSize',44, 'LineWidth',1.25);
axis tight;
set (get(gca,'XAxis'), 'LineWidth', 4)
set (get(gca,'YAxis'), 'LineWidth', 4)

print -depsc 'dl1_barycentre_size.eps'

keyboard;


