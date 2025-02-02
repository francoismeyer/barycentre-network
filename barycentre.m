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
%%_Module_Name :
%%
%%_Description : Compute the barycentre graph from a sample of T graphs, each of which of size n
%%               The distance is l2 norm of the vector of eigenvalues of the normalized graph Laplacian.
%%               G is a tensor of T adjacency matrices of size nb x n: G(1:n,1:n,1:T)
%%               M is the number of communities in the SBM.
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

function [Ahat, Epd] = barycentre (G,M,n,T)
  %%
  %%

  %% set to 1 if you want some nice plots
  DISPLAY = 0;

  %%
  %% Sample mean adjacency matrix; we need it for the Soules basis
  %%

  mA = mean (G,3);

  %% if (DISPLAY)
  %%   figure;imagesc(mA); title ('sample mean A'); axis square, axis off; colormap ('jet'); colorbar;
  %% end

  oneshot = (T == 1);

  if (oneshot)

    %% compute the degree vector, and the normalization matrix = 1./sqrt(degree matrix)

    d = sum(mA,2);
    D = diag(d);

    D2 = diag(sqrt(d));
    iD2 = diag(1./sqrt(d));

    %% make it symmetric, for good measure

    An = iD2 * mA * iD2;
    An = 0.5.*(An + An');

    %% Compute the eigenvalues of \hat{A}

    lambda = 1 - eig(An);								

    %% sort the eigennvalues by increasing order, and get rid of the bulk eigenvalues

    eigL  = sort(lambda,'ascend');
    eigL2 = eigL (1:M);

  else

    %% compute the  eigenvalues of each normalized adjacency matrix
    
    lamb = zeros (n,T);
    Kn   = zeros (n,n);
    K    = zeros (n,n);
    
    for t=1:T
      K = G(:,:,t);
      
      dK   = sum(K,2);
      iDK2 = diag(1./sqrt(dK));
      
      Kn = iDK2 * K * iDK2;
      Kn = 0.5.*(Kn + Kn');
      
      lamb (:,t) = eig(Kn);
      
      fprintf ('\r %d eigensolve completed',t);
    end

    %% this is the mean spectrum
    
    lambda = 1 - mean (lamb,2);
    
    %% sort the eigennvalues by increasing order, and get rid of the bulk eigenvalues

    eigL  = sort(lambda,'ascend');
    eigL2 = eigL (1:M);

  end

  %%==================================================================================================
  %%
  %% Compute the best Soules basis
  %%
  %%==================================================================================================

  %% first eigenvector of the normalized graph Laplacian

  x = (1/sqrt(n)).*ones (n,1);

  %%
  %% allocate memory
  %%

  Qn = zeros (n, n);

  %%
  %% build the Soules basis by exploring the vector of cuts
  %%

  Qn = topdown (mA,x);

  %% display the Soules basis

  if (DISPLAY)
    
    f = figure;
    hold on;
    plot (Qn(:,1), 'color','cyan',   'linewidth',3)
    plot (Qn(:,2), 'color','blue',  'linewidth',3)
    plot (Qn(:,3), 'color','green', 'linewidth',4)
    plot (Qn(:,4), 'color','red',  'linewidth',2)
    
    lalegende = legend('\boldmath{$\psi$}$_1$','\boldmath{$\psi$}$_2$','\boldmath{$\psi$}$_3$','\boldmath{$\psi$}$_4$');
    set(lalegende,'Interpreter','latex');
    set(lalegende,'FontSize',22);
    set (lalegende,'FontName','Times New Roman');

    set (gca,'FontSize',20);
    set (gca,'FontName','Times New Roman');

    axis tight

    print -depsc 'les soules.eps';
    drawnow('update');
  end

  %%==================================================================================================
  %%
  %% remove the columns of the Soules basis associated to these eigenvalues  and compute EM (see paper)
  %% We also threshold EM to get an estimate of the blocks geometry, EP.
  %% Finally, we remove the diagonal in EP, and rerurn Epd.
  %%==================================================================================================

  Q = Qn (:,1:M);
  EM  = Q*Q';			
  EP  = EM > 1e-6;
  Epd = EP - diag (diag(EP));

  %% if (DISPLAY)
  %%   figure;imagesc(EM);title ('EM'); colormap('jet'); colorbar;
  %% end

  %%==================================================================================================
  %%
  %% estimate of the degree
  %% the key idea is to compute the average of the degrees in each block
  %%
  %%==================================================================================================

  d1     = sum (mA,2);
  blcksz = sum(Epd,2);
  dbar   = Epd*d1;
  dbar   = dbar./blcksz;

  Dbar   = diag(dbar);
  D2bar  = diag(sqrt(dbar));

  %%==================================================================================================
  %%
  %%   test the reconstruction of the normalized Laplacian using the first M Soules eigenvectors
  %%   where M is the number of communities
  %%   reconstruct an adjacency matrix using the expected degree matrix
  %%
  %%==================================================================================================

  L2 =  (Q * diag(eigL2) *Q');
  A2 =  D2bar * (EM - L2) * D2bar;

  A2 = 0.5.*(A2 + A2');
  A2 = A2 - diag(diag(A2));
  A2 = A2.* Epd;

  if (DISPLAY)
    figure;imagesc(A2); colormap('jet');colorbar;
    axis square, axis off;
    print -depsc 'truncated reconstruction.eps'
  end

  Ahat = A2;

end
