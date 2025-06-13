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
%%_Module_Name : eiglaplace.m
%%
%%_Description : compute the eigenvectors of the normalized
%%               Laplacian. We compute only the first nvec vectors, associated with 
%%               the nvec smallest eigenvalues. 
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


function [V] = eiglaplace (A, nvec)
  
  %%   We have a symmetric weight matrix A that describes the
  %%   graph topology and weights. Our next step consists in  
  %%   normalizing this matrix. We use the Laplace-Beltrami
  %%   normalization

  N = size (A,1);

  %% compute the degree of each node, sum across the rows

  deg = sum (A,2);                  

  %% make the degrees along the diagonal

  degree = speye (length (deg));

   %% set the diagonal of degree to the (deg)^{-1/2} we do this by
   %% looking at the matrix as a linear array, and advancing by exactly
   %% nColumns +1, so that at each time we hit the diagonal.
   % Normalized Beltrami Laplacian

  degree (1:N+1:end) = 1./sqrt(deg(deg ~= 0));

  %% Calculate the symmetric normalized adjacency matrix. The matrix
  %% corresponds to a kernel version of the probality transition
  %% K =  D^{1/2} A ^{-1/2} 

  K = degree * A * degree;

  clear degree;

  %% in principle, it should be symmetric a this point of the game,
  %% but we can enforce this numerically.

  K = (K + K')*.5;

  %%
  %%   solve eigenvalue problem
  %%
  
  %%  disp ('computing eigenvectors...');
  %%   timex = cputime;

  %% the nvec eigenvectors are stored in columns in eigvec
  %% columns of eigvec are the eigenvectors
  
  [eigvec, eigval]= eigs (K, nvec);

  %%  fprintf(1,'eigenvectors computation took %g s\n', cputime -timex);

  %%
  %%  eigenvalues of the normalized Laplaciance
  %%

  eigval = 1 - diag(eigval);

  %% Sort the eigenvalues from small to big, do the same on the eigenvectors 

  [lTemp lIdx] = sort(eigval);
  eigvec = eigvec(:,lIdx);

  V = eigvec (:,1:nvec);
  
  return
end
