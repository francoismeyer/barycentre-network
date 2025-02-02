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
%%
%%              
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
%%=========================================================================================
function [mae] = barycentreSBM (n,M, equalSize, equalChance, T)

%%==========================================================================================

DISPLAY = 0;
ONESHOT = (T == 1);

%%==========================================================================================
%%
%% generate a sample from a hard-wired SBM
%%
%%==========================================================================================

In = eye(n,n);
J  = ones (n,n);

[P,miniP,  sizeofBlock] = genSBM (n, equalSize, equalChance, M);
%%
%% display the chance matrix for the SBM
%%

if (DISPLAY)
  figure;imagesc(P); title ('P'); colormap ('jet'); colorbar;
  axis square, axis off;
  print -depsc 'P.eps'
end

miniPdec = sort(diag(miniP), 'descend');

p = miniPdec(1);
q = miniP(end);

G = zeros (n, n ,T);

if (ONESHOT)
  %% only one realization
  %%
  A = IHfast (P);

  if (DISPLAY)
    figure;imagesc(ones(size(A)) - A); colormap ('gray');
    axis square, axis off;
    print -depsc 'A.eps'
  end

else  
  %%
  %% generate T realisations of the SBM
  %%

  for t=1:T
    G(:,:,t) = IHfast (P);
  end
end

if (ONESHOT)
  [Ah, Epd] = barycentre (A,M,n,1);
else
  [Ah, Epd] = barycentre (G,M,n,T);
end

del = abs((P - Ah).* Epd);

mae = norm (del,1)/numel (P);

if (DISPLAY)
  figure;imagesc (del); colormap('jet'); colorbar;
  axis square, axis off;
  print -depsc 'residual error.eps'
end

end
