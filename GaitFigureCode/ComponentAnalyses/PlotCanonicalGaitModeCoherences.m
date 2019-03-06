function PlotCanonicalGaitModeCoherences( newData, cmp, corder )

if nargin < 2
    cmp = flipud(gray(2^16));
end

if nargin < 3
    corder = cbrewer('qual', 'Dark2', 4);
end

%% Set preferences for computation

% Select whether to use smoothed data
smoothing = true;

% Define lists of variables to extract
limbList = {'L1','L2','L3','R1','R2','R3'};
phaseVarList = strcat('InstantaneousPhase_', limbList, 'y');
if smoothing
    phaseVarList = strcat('smooth_', phaseVarList);
end

% Set the minimum forward speed
vfMin = 1;

% Define list of gaits
gaitList = {'tripod','left tetrapod','right tetrapod','wave'};

% Define offsets for each canonical gait
psi = [
    0, 1/2, 0, 1/2, 0, 1/2;
    1/3, 2/3, 0, 0, 1/3, 2/3;
    2/3, 0, 1/3, 0, 1/3, 2/3;
    1/6, 1/3, 1/2, 2/3, 5/6, 0;
    ];

% Set the intervals for conditioning on forward speed
vfBinEdges = (1:1:32);

% Set the query points for the coherences
xq = 0:0.001:1;

%% Compute the coherences for each mode

% Extract needed data from the table
vf = newData.smooth_forwardSpeed_mmPerSec;
Phi = newData{vf>vfMin, phaseVarList};
id = newData.videoID(vf>vfMin);
vf = vf(vf>vfMin);

% Take phases modulo 2 pi
Phi = mod(Phi, 2*pi);

% Find the number of bins
nBin = length(vfBinEdges)-1;

% Discretize the forward speeds to facilitate conditioning
vfDisc = discretize(vf, vfBinEdges);

% Eliminate NaN values
notNan = ~any(isnan(Phi),2) & ~isnan(vfDisc);
Phi = Phi(notNan,:);
id = id(notNan);
vf = vf(notNan);
vfDisc = vfDisc(notNan);

nf = length(unique(vfDisc));
ni = length(unique(id));

% Compute resultants for projections into each mode
numPsi = size(psi,1);
phiMean = nan(size(Phi,1), numPsi);
for ind = 1:numPsi
    phiMean(:,ind) = nanmean(exp(1i * (Phi + 2*pi * psi(ind,:))),2);
end

% Compute the modulus and argument of each resultant
r = abs(phiMean);
theta = angle(phiMean);

%% Plot marginal PDFs of coherences

fPDFr = nan(length(xq),numPsi);
for ind = 1:numPsi
    fPDFr(:,ind) = ksdensity(r(:,ind), xq, 'kernel','epanechnikov');
end

figure('Position',[200,500,500,700],'WindowStyle','docked');
hold on;
set(gca, 'colororder', corder);
plot(xq, fPDFr, 'linewidth', 2);
axis('square');
legend(gaitList, 'location','eastoutside');
xlabel('coherence');
ylabel('pdf');
ConfAxis('fontSize', 16);

%% Plot marginal PDFs of coherences, conditioned on forward speed

% Compute conditional PDFs
fPDFgivenVf = nan(nBin, length(xq), numPsi);
for vInd = 1:nBin
    idx = vfDisc == vInd;
    if any(idx)
        for psiInd = 1:numPsi
            fPDFgivenVf(vInd, :, psiInd) = ksdensity(r(idx, psiInd), xq, 'kernel','epanechnikov');
        end
    end
end

% Plot the PDFs
maxP = round(max(fPDFgivenVf(:)));
for ind = 1:numPsi
    figure('Position',[200,500,500,700],'WindowStyle','docked');
    imagesc(vfBinEdges(1:end-1), xq, fPDFgivenVf(:,:,ind)');
    
    axis('xy','square','tight');
    xlabel('v_{||} (mm/s)');
    ylabel('coherence');
    title(gaitList{ind});
    cbar = colorbar;
    ylabel(cbar, 'pdf given v_{||}');
    caxis([0 maxP]);
    colormap(cmp);
    ConfAxis('fontSize', 16);
end

%% Plot joint PDFs of mode coherences

xq = 0:0.01:1;
% [xx,yy] = meshgrid(xq);
% xx = xx(:);
% yy = yy(:);

indRef = 1;
fPDFjoint = nan(length(xq)-1,length(xq)-1, numPsi-1);
for ind = 2:numPsi
    %     fPDFjoint(:,:,ind) = reshape(ksdensity([r(:,indRef), r(:,ind)], [xx(:), yy(:)], 'kernel','epanechnikov'),length(xq),length(xq));
    fPDFjoint(:,:,ind) = histcounts2(r(:,indRef), r(:,ind), xq, xq, 'normalization','count')';
end

maxP = max(abs(fPDFjoint(:)));
for ind = 2:numPsi
    figure('Position',[200,500,500,700],'WindowStyle','docked');
    %     imagesc(xq,xq, log10(fPDFjoint(:,:,ind)));
    imagesc(xq,xq, fPDFjoint(:,:,ind));
    axis('xy','equal','tight');
    xlabel(sprintf('%s mode coherence', gaitList{indRef}));
    ylabel(sprintf('%s mode coherence', gaitList{ind}));
    cbar = colorbar;
    ylabel(cbar, 'pdf');
    %     caxis([0 maxP]);
    xlim([0 1]);
    ylim([0 1]);
    colormap(cmp);
    ConfAxis('fontSize', 16);
end

%% Plot average coherence in each mode as a function of forward velocity

mu = nan(nBin,numPsi);
ciUpper = nan(nBin, numPsi);
ciLower = nan(nBin, numPsi);

for ind = 1:numPsi
    s = accumarray([id,vfDisc], r(:,ind), [ni,nf],  @nanmean);
    ci = bootci(1000, {@(x) squeeze(nanmean(x,1)), s});
    
    mu(:,ind) = nanmean(s,1);
    ciLower(:,ind) = ci(1,:);
    ciUpper(:,ind) = ci(2,:);
end

figure('Position',[200,500,500,700],'WindowStyle','docked');
PlotAsymmetricErrorPatch(vfBinEdges(1:end-1)', mu, ciLower, ciUpper, corder);
axis('square');
ylim([0 1]);
legend(gaitList, 'location','eastoutside');
xlabel('v_{||} (mm/s)');
ylabel('mean coherence');
ConfAxis('fontSize', 16);

%% Plot percentage of samples instantaneously best described by each mode as a function of forward velocity

mu = nan(nBin,numPsi);
ciUpper = nan(nBin, numPsi);
ciLower = nan(nBin, numPsi);

for ind = 1:numPsi
    maxCoherence = double(all(r(:,ind) > r(:,(1:4)~=ind),2));
    s = accumarray([id,vfDisc], maxCoherence, [ni,nf],  @nanmean);
    ci = bootci(1000, {@(x) squeeze(nanmean(x,1)), s});
    
    mu(:,ind) = nanmean(s,1);
    ciLower(:,ind) = ci(1,:);
    ciUpper(:,ind) = ci(2,:);
end

figure('Position',[200,500,500,700],'WindowStyle','docked');
PlotAsymmetricErrorPatch(vfBinEdges(1:end-1)', mu, ciLower, ciUpper, corder);
axis('square');
legend(gaitList, 'location','eastoutside');
xlabel('v_{||} (mm/s)');
ylabel('instantaneous fraction of points best described by mode');
ConfAxis('fontSize', 16);
ylim([0 1]);

% Plot cumulative percentages
cumulativePct = cumsum(mu,2);
figure('Position',[200,500,500,700],'WindowStyle','docked');
hold on;
set(gca, 'colororder', corder);
% Plot cumulative percentages as lines
plot(vfBinEdges(1:end-1)',cumulativePct, 'linewidth',2);

% Fill in patches
x = [vfBinEdges(1:end-1)'; flipud(vfBinEdges(1:end-1)')];
cumulativePct = [zeros(size(cumulativePct,1),1),cumulativePct];
for ind = 2:numPsi+1
    y = [cumulativePct(:, ind-1); flipud( cumulativePct(:, ind) )];
    patch(x, y, 1, 'FaceColor', corder(ind-1,:), 'FaceAlpha', 1/2, 'EdgeColor', 'None');
end
axis('square','tight');
legend(gaitList, 'location','eastoutside');
xlabel('v_{||} (mm/s)');
ylabel('instantaneous fraction of points best described by mode');
ConfAxis('fontSize', 16);
ylim([0 1]);

%% Plot the distribution of coherences for those points instantaneously best described by a given mode

% Set the query points for the coherences
xq = 0:0.001:1;

fPDF = nan(length(xq),numPsi);
for ind = 1:numPsi
    maxCoherence = all(r(:,ind) > r(:,(1:4)~=ind),2);
   fPDF(:,ind) = ksdensity(r(maxCoherence,ind), xq, 'kernel','epanechnikov');
end

figure('Position',[200,500,500,700],'WindowStyle','docked');
hold on;
set(gca, 'colororder', corder);
plot(xq, fPDF, 'linewidth', 2);
axis('square');
legend(gaitList, 'location','eastoutside');
xlabel('distribution of coherence for points best described by given mode');
ylabel('pdf');
ConfAxis('fontSize', 16);

end
