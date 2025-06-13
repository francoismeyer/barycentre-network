%%==============================================================================
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
%%_Module_Name: barycentreSchool
%%
%%_Description: This is to test the barycentre function with the elementary 
%%              school dataset.  We test the reconstruction of the normalized
%%              graph Laplacian based on a Soules basis and the smallest eigenvalues
%%              of the normalized graph Laplacian.
%%
%%              One can set the number of communities to M = 5 grades, or M = 10 classes.
%%
%%  Input:
%%         [min0,min1]: interval over which we get a sample every 5 minutes
%%         aggMin: length of the window (in minute) over which we aggregate the networks.
%%              
%% 
%%
%%_References:
%%
%%_Remarks: None
%%
%%_Author:                 Francois G. Meyer
%%
%%_Revisions History: 2025 Initial keying
%%
%%==============================================================================


function [code] =  barycentreSchool (min0, min1, aggMin)

  addpath ('./generate graph');
  addpath ('./soules');
  addpath ('./util');

  %% to get nice figures and pdf files
  
  DISPLAY = 1;
  
  %% load the aggregated time series of networks

  [Tensor, leseleves, miniP] = schoolAggreg (min0, min1, aggMin);
  
  T = size(Tensor,3);
  
  n = size (Tensor,1);		%% number of students/nodes

  %% M = 5;			%% number of grades
  M = 10;			%% number of classes

  [Qn, mA, Epd, EA] = barycentre (Tensor, M, n, T);

  q = miniP(end);		%% accross communities edge probability
  
  %% remove spurious entries in the adjacency matrix 

  A = makeitbinary (mA, (M-1)*q);

  if (DISPLAY)
    f = drawgraph(EA, leseleves);
    name = sprintf ('graph_EA_%03d_%03d_%03d', aggMin, min0, min1);
    fig2pdf (f, name);
    drawnow('update');
  
    f = drawgraph(A, leseleves);
    name = sprintf ('bary_graph_%03d_%03d_%03d', aggMin, min0, min1);
    fig2pdf (f, name);
    drawnow('update');

    f = figure;imagesc(A); axis square, axis off; colormap ('jet');colorbar;
    name = sprintf ('A_bary_%03d_%03d_%03d', aggMin, min0, min1);
    fig2pdf (f, name);
    drawnow('update');
  end
  

  if (DISPLAY)
    %% 
    %% the next lines are for display purposes only. 
    %%
    A2  = A;
    A2  = 0.5*(A2+A2');

    f=figure;imagesc(ones(size(A2)) - A2); colormap ('gray');
    axis square, axis off;
    name = sprintf ('full_recon_EA_%03d_%03d_%03d', aggMin, min0, min1);
    fig2pdf (f, name);
    drawnow('update');
  end
  
  if (DISPLAY)
     deleves = diff(leseleves);	      %% detect changes in class number
     where   = find(deleves == 1)+1;
 
    S = Qn;
    f = figure;
    hold on;

    f.Theme = "light"; 		%% white background
    newDefaultColors = jet(M);	%% get gazzilion colors
    set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');

    plot (S(:,1:10),'linewidth', 3); %% plot the 10 eigenvectors all at once

    xline (where, 'Color','k', 'linestyle','-.','linewidth', 1);
    
    lalegende = legend('\boldmath{$\psi$}$_1$',...
		       '\boldmath{$\psi$}$_2$',...
		       '\boldmath{$\psi$}$_3$',...
		       '\boldmath{$\psi$}$_4$',...
		       '\boldmath{$\psi$}$_5$',...
		       '\boldmath{$\psi$}$_6$',...
		       '\boldmath{$\psi$}$_7$',...
		       '\boldmath{$\psi$}$_8$',...
		       '\boldmath{$\psi$}$_9$',...
		       '\boldmath{$\psi$}$_{10}$');
    
    set(lalegende,'Interpreter','latex');
    set(lalegende,'FontSize',22);
    set (lalegende,'FontName','Times New Roman');

    set (gca,'FontSize',20);
    set (gca,'FontName','Times New Roman');

    axis tight

    name = sprintf ('les_soules_%03d_%03d_%03d', aggMin, min0, min1);
    fig2pdf (f, name);
    drawnow('update');
  end

  code = true;

end

function [C] = makeitbinary (A,low)
  C = max (0,A -low);
end
