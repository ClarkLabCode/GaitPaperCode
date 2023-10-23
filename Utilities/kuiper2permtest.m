function [p, pci, V, V0] = kuiper2permtest(countx,county, nperm, alpha_level)
% kuiper2permtest: two-sample Monte Carlo permutation test using the Kuiper
% test statistic

if nargin < 3
    nperm = 1000;
end

if nargin < 4
    alpha_level = 0.05;
end

countxy = [countx, county];

n = size(countx,2);
m = size(county,2);
nm = n + m;

%% Compute the value of the test statistic for the actual samples

a = sum(countx,2);
b = sum(county,2);

cdfx = cumsum(a) / sum(a);
cdfy = cumsum(b) / sum(b);

V = max(cdfx - cdfy) + max(cdfy - cdfx);

%% Monte Carlo permutation test

V0 = nan(nperm,1);

for ind = 1:nperm
    
    p = randperm(nm);
    
    a = sum(countxy(:,p(1:n)),2);
    b = sum(countxy(:,p(n+1:end)),2);
    
    cdfa = cumsum(a) / sum(a);
    cdfb = cumsum(b) / sum(b);
    
    V0(ind) = max(cdfa - cdfb) + max(cdfb - cdfa);
    
end

%% Compute the p-value

% Compute the p-value
x = nnz(V0>V);
p = x / nperm;

% Compute the Clopper–Pearson confidence interval for the p-value
pci = [betainv(alpha_level/2,x, nperm-x+1); betainv(1-alpha_level/2,x+1,nperm-x)];
pci(isnan(pci)) = 0;

end