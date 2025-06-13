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
%%_Module_Name: arrange.m
%%
%%_Description: The function takes a sequence of T = size(G,3) SBM
%% graphs, and use multiway spectral clustering to discover the labels,
%% then realign the graphs so that all classes match. The spectral
%% clustering needs the number of classes = M  
%%
%% The function returns the aligned graph in AG. 
%%
%% code is FALSE if the realignment failed.
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
function [AG, code] = arrange (G,M,p)

  n = size (G,1);		%% graph size
  T = size(G,3);		%% number of graphs

  K = zeros (n,n);
  U = zeros (n,n);
  Ka = zeros (n,n);

  AG = zeros (n,n,T);
  
  for t=1:T

    K = G(:,:,t);

    %% multiway clustering of the graph K

    [X] = multiwaycluster (K, M);

    %% compute the permutation to organize the nodes according to
    %% their communities, with the communities ordered by decreasing
    %% size

    U = aligne (X,M);

    %% permute

    Ka = U*K*U';

    %% make it symmetric
    
    Ka = 0.5*(Ka+Ka');
    
    AG(:,:,t) = Ka;

  end

  %%
  %% sanity check to ensure that all graphs are properly aligned. If not, the rest of the story is useless.
  %%

  Alinear = reshape(AG, n*n, T);

  vA = var (Alinear,0,2);

  lavar = 6 * max (p.*(1-p));
  
  if (max (vA) > lavar)
    %% the realignment failed
    fprintf('\n\t Warning: alignment failed\n');
    code = false;
    keyboard;
  else
    %% the variance after realignment is about the same as the variance in the SBM model
    code = true;			
  end
  
end


