%%=========================================================================================
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
%%_Description : This is to test the barycentre function.
%%               The core function is barycentre. We generate some random realisations of SBMs,
%%               and we compute the mean absolute error (l1 norm normalized)
%%  INPUT:
%%
%%     n: size of the graphs
%%
%%     M: number of eigenvalues that are used to reconstruct the graph using the truncated Soules basis.
%%
%%     T: number of graph in the sample
%%
%%     equalsize: if TRUE then all the blocks have the same size (n/M) in the SBM
%%
%%     equalchance: if TRUE then all the blocks have the same edge probability density
%%                  a balanced SBM is equivalent to equalsize = equalchance = TRUE
%%
%%     DISPLAY: set it to TRUE to get nice figures
%%
%%
%%  OUTPUT:
%%
%%     mse: mean squared reconstruction error
%%
%%     code: if FALSE, then the alignment failed, and none of the output is sensible.
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
%%=========================================================================================

function [mse,code] = barycentreSBM (n,M, equalSize, equalChance, T)

  addpath ('./generate graph');
  addpath ('./soules');
  addpath ('./util');


  DISPLAY = 1;
  ONESHOT = (T == 1);

  %%==================================================================================================
  %%
  %% generate a sample from a hard-wired SBM
  %%
  %%==================================================================================================

  In = eye(n,n);
  J  = ones (n,n);

  [P,miniP,  sizeofBlock] = genSBM (n, equalSize, equalChance, M);

  %%
  %% build the vector of classes
  %%

  classes = zeros (n,1);

  deb = 1;
  for (m=1:M)
    fin = deb + sizeofBlock(m) - 1;
    idx = [deb:fin]; 
    classes(idx) = m;
    deb = fin + 1;
  end

  %%
  %% display the chance matrix for the SBM
  %%

  if (DISPLAY)
    figure;imagesc(P); colormap ('jet');colorbar;
    axis square, axis off;
    fig2pdf (gcf,'P');
  end

  miniPdec = sort(diag(miniP), 'descend');

  p = miniPdec(1);
  q = miniP(end);

  G = zeros (n, n ,T);

  A = zeros (n,n);
  B = zeros (n,n);
  U = zeros (n,n);
  Q = zeros (n,n);

  if (ONESHOT)
    %% one realization
    %%
    A = IHfast (P);

    U = ranPerm(n);		

    B = U'*A*U;

    A = round (0.5*(B + B'));

    if (DISPLAY)
      figure;imagesc(ones(size(A)) - A); colormap ('gray');
      axis square, axis off;
      fig2pdf (gcf,'A');
    end
    
  else  
    %%
    %% generate T realisations of the SBM; the node labels are  randomly
    %% permuted independently for each realization 
    %%

    for t=1:T

      A = IHfast (P);		%% one SBM with the edge probabilty P

      U = ranPerm(n);		%% randomly permute the nodes

      B = U'*A*U;
      A = round(0.5*(B + B'));	%% make it symmetric, and round to 0/1
      
      G(:,:,t) = A;		%% save the realization
    end
  end

  %%
  %%   we first perform a coarse realignment to assign the nodes to their respective classes
  %% 

  if (ONESHOT)
    [AG, reussi] = arrange (A,M, miniP);
  else
    [AG, reussi] = arrange (G,M, miniP);
  end

  %%
  %% we compute the barycentre
  %%

  if (reussi)
    [Qn, Atrunc,Epd] = barycentre (AG,M,n,T);
    code = true;
  else
    fprint('\n Warning: aligned failed in barycentreSBM failed');
    code = false;
    mse = 0;
    return;
  end

  clear A B U;

  %%
  %% relabel the classes to have the largest class in the top left corner
  %%

  Q = aligne (classes,M);

  P2 = Q*P*Q';

  del = abs((P2 - Atrunc).* Epd);

  mse = norm (del,"fro");	        %% frobenius norm

  mse = (mse*mse)/numel (del);		%% normalize

  %%
  %% display the Soules basis
  %%

  if (DISPLAY)
    
    S = Q'* Qn;
    f = figure;
    hold on;

    f.Theme = "light"; 		       %% white background

    plot (S(:,1), 'color','cyan',   'linewidth',3)
    plot (S(:,2), 'color','blue',  'linewidth',3)
    plot (S(:,3), 'color','green', 'linewidth',4)
    plot (S(:,4), 'color','red',  'linewidth',2)
    
    xline ([1 63 210 315 512] , 'Color','k', 'linestyle','-.','linewidth', 1);

    lalegende = legend('\boldmath{$\psi$}$_1$',...
		       '\boldmath{$\psi$}$_2$',...
		       '\boldmath{$\psi$}$_3$',...
		       '\boldmath{$\psi$}$_4$');
    
    set(lalegende,'Interpreter','latex');
    set(lalegende,'FontSize',22);
    set (lalegende,'FontName','Times New Roman');

    set (gca,'FontSize',20);
    set (gca,'FontName','Times New Roman');

    axis tight

    fig2pdf (f,'lesSoules');
    drawnow('update');

  end

  if (DISPLAY)

    %% the next lines are for display purposes only. We take the
    %% permutation associated with rearanging P with larges classes in the
    %% top left corner, and transpose it to re-order Atrunc, with the hope
    %% that we get something that looks like the original P.
    
    A2  = Q'*Atrunc*Q;
    A2  = 0.5*(A2+A2');

    figure;imagesc(A2); colormap('jet');colorbar;
    axis square, axis off;
    fig2pdf (gcf,'truncatedReconstruction');
    
    Epd2 = Q'*Epd*Q;

    del2 = abs(P - A2);
    del2 = del2.* Epd2./numel (del2);		%% clipp and normalize
    
    figure;imagesc (del2); colormap('jet'); colorbar;
    axis square, axis off;
    
    fig2pdf (gcf,'residualError');
  end

end
