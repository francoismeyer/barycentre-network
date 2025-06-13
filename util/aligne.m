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
%%_Module_Name: aligne.m
%%
%%_Description:   The function reorganizes the nodes according to
%%                their classes given by X. The nodes in the largest
%%                class are grouped together, and then all the nodes
%%                in the second largest class, etc. 
%%
%%  The unitary matrix V is the permutation that achieves this
%%  rellabeling.
%%
%%  the vector X contains the class label for each node index i:
%%  X(i) = class of node i The function groups spatially the nodes into
%%  their respective classes, so that the largest class comes first.
%%
%%
%%  Explanation of the variables: mapClass (c) is the new class
%%  index of class c, so that size(maplcass (c)) >= size(maplcass
%%  (c+1)) maplcass (c)) < maplcass (c+1)
%%
%%    sortX(i) is the new class of node i defined by mapclass:
%%
%%     [--- class 1 ----|--- class 2 ---|- class M -]
%%
%%    psi(i) is the index of node i, if we organize nodes by
%%    increasing new class order -- or decreasing class sizes
%%
%%    V is the permutation matrix that assign node i to its new
%%    location in [1,n] according to psi
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
function [V] = aligne (X,M)
  
  nbins = max(X) + 1;

  [counts,~] = histcounts(X, nbins);

  [~, mapClass] = sort (counts,'descend');
  
  n = size(X,1);
  
  for (t = 1:n)
    sortX(t) = mapClass (X(t));
  end

  [~,psi] = sort(sortX,'ascend');
  
  In = eye(n, n);

  V = In(psi,:);
end
