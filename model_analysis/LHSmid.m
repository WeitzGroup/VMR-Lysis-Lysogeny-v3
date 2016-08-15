% Part of the code used in:
% Weitz et al. Lysis, Lysogeny, and Virus-Microbe Ratios
% 
% From https://github.com/WeitzGroup/VMR-Lysis-Lysogeny-v2
% MIT License

function LHSample = LHSmid(nS,pmin,pmax,r)
%function LHSample = LHSmid(nS,pmin,pmax,r)
%
%This function performs Latin Hypercube Sampling with sampling the midpoint
%of the selected hypercubes.  
%LHSample is a r*nS x nPar array of sampled points where nPar is the number
%of parameters.
%pmin is a vector of minima for the sampled parameter range.  pmax is a vector
%of maxima for the sampled parameter range.  r is on optional input that is
%the number of iterations of the LHS.  Set >1 to perform simple resampling.
%--------------------------------------------------------------------------
assert(numel(pmin) == numel(pmin), ...
    'The number of elements in pmin and pmax must be the same.' );
if nargin == 3
    r=1;
end
final = gensample(nS,pmin,pmax);
    LHSample=zeros(r.*nS,length(pmin));
for k=0:r-1

    LHSample2 = zeros(size(final));
    for m = 1:length(pmin)
        LHSample2(m,:) = final(m,randperm(size(final,2)));
    end
    LHSample((1+k*nS):(k+1)*nS,:) = LHSample2';
end
end
function final = gensample(nS,pmin,pmax)

nV = length(pmin);

dX =(pmax-pmin)/nS;
test = 0.5:nS-.5;
final = repmat(reshape(dX,nV,1),1,nS).*repmat(test,nV,1)...
    +repmat(reshape(pmin,nV,1),1,nS);
end
