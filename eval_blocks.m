%% this script evaluates the error between the estimated Frechet mean
%% using the Soules basis, and the edge probability used to generate
%% the SBMs
%% We vary the  block size

labarre='–\|/';		% for display; ignore

n = 512;			     %% graph size

equalSize   = 1;
equalChance = 1;

nValues = 5;                        %% how many values of the sample size within the range of sample sizes
                                     %% in plain English: number of ticks on the x-axis
T = 1;				     %% number of samples

maxTrials = 64;                      %% number of experiments for a fixed sample size

fprintf ('\n We will process %d (x %d) graphs', nValues, maxTrials);

%% generate the block sizes, and number of blocks

dl1    = zeros (maxTrials, 1);       %% l1 distances for a trial, for a fixed sample size
histl1 = zeros (maxTrials, nValues); %% l1 distances for all trials, for a fixed sample size
meanl1 = zeros (1, nValues);         %% mean l1 distances (computed over maxTrials) for all sample sizes

nblcks = zeros(nValues,1);
blckSz = zeros(nValues,1);

nblcks(1)  = 2;
blckSz (1) = floor(n/2);
for b = 2:nValues 		     %% for increasing block size
  nblcks (b) = 2*nblcks (b-1);
  blckSz (b) = floor(blckSz (b-1)/2);
end

fprintf ('\n');
lastsize = 0;

for b = 1:nValues 		     %% for increasing block size
  
  M = nblcks (b);		    %% number of blocks

  ix = 1;  			     %% for display
  for ntrials=1:maxTrials            %% we run maxTrials independent experiments for the sample size N

    %%  compute the difference between the barycentre and the population mean, P.

    dl1 (ntrials) =  barycentreSBM (n, M, equalSize, equalChance, T);
    
    %% next two lines for display
    fprintf(1,'\b%c', labarre(ix));
    ix = mod(ix,4)+1;
  end
  
  histl1 (:,b) = dl1;
  meanl1 (:,b) = mean (dl1);
  
  fprintf(repmat('\b', 1, lastsize));
  lastsize = fprintf('\n%2d\n\n', b);
end

   % = hamming distance (population graph median, sample graph median)

coeff3   = polyfit (log(nblcks),log (meanl1),3)
linereg3 = polyval(coeff3, log(nblcks))';

%
% figure 2: hamming distance (population graph median, sample graph median)
%

save ("evalblock.mat","coeff3", "linereg3", "nblcks", "nValues", "meanl1", "histl1");

f = figure;f.Position = [2000 1000 1000 1000]; hold on;
p1 = boxplot(log(abs(histl1)), log(nblcks),'Labels', {'2','4','8','16','32'}, 'colors',[1 0 0]);
set(p1,'LineWidth',4);
p2      = plot ([1:nValues],linereg3,'blue');
set(p2,'LineWidth',4);
set (gca,'FontSize', 36);
set (gca,'FontName','Times New Roman');
set(get(gca, 'XLabel'),'String','M','FontName','Times New Roman','FontSize',48, 'LineWidth', 4);
set(get(gca, 'YLabel'),...
    'String','$\log_{10}\big(n^{-2}\|{\bf P}-\widehat{{\mu\mkern -10.2mu\mu}}_N[{\mathrm I}\mkern -4mu{\mathrm P}]\|_1\big)$',...
    'Interpreter','latex', 'FontName','Times New Roman','FontSize', 44, 'LineWidth',1.25);
axis tight;
set (get(gca,'XAxis'), 'LineWidth', 4)
set (get(gca,'YAxis'), 'LineWidth', 4)

print -depsc 'dl1_barycentre_block.eps'

keyboard;


