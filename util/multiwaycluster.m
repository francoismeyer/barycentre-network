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
%%_Module_Name : multiwaycluster.m
%%
%%_Description : computes the eigenvectors of the normalized
%%               Laplacian, and perform a multiway spectral clustering based the
%%               bottom M eigenvectors. See clusterassign.m for details.
%%              
%%_References :
%%
%% [1] Eldén, L., 2024. Multiway Spectral Graph Partitioning: Cut
%% Functions, Cheeger Inequalities, and a Simple Algorithm. SIAM
%% Journal on Matrix Analysis and Applications, 45(1), pp.112-133.  
%%
%% [2] Damle, A., Minden, V. and Ying, L., 2018. Simple, direct and
%% efficient multi-way spectral clustering, Inf. Inference J. IMA,
%% 8(1), pp.181-203.
%%
%% [3] Zha, H., He, X., Ding, C., Gu, M. and Simon, H., 2001. Spectral
%% relaxation for k-means clustering. Advances in neural information
%% processing systems, 14. 
%%              
%% INPUT:
%%
%%    A is n x n adjacency matrix of a community structure network
%%    M: upperbound on the number of communities
%%
%% OUTPUT
%%    label returns the cluster assignment vectors, associated to each of
%%    the n nodes.
%%
%%
%%_Remarks : None
%%
%%_Author :  Francois G. Meyer
%%
%%_Revisions History: 2025 Initial keying
%%
%%==============================================================================

function [label] = multiwaycluster (A, M)
  
  [V]   = eiglaplace (A,M);
  label = clusterassign(V);
  
end
