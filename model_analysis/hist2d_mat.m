% Part of the code used in:
% Weitz et al. Lysis, Lysogeny, and Virus-Microbe Ratios
% 
% From https://github.com/WeitzGroup/VMR-Lysis-Lysogeny-v3
% MIT License


function mHist = hist2d_mat (mX, vYEdge, vXEdge)
nCol = size(mX, 2);
if nCol < 2
    error ('mX has less than two columns')
end

nRow = length (vYEdge)-1;
nCol = length (vXEdge)-1;

vRow = mX(:,1);
vCol = mX(:,2);

mHist = zeros(nRow,nCol);

for iRow = 1:nRow
    rRowLB = vYEdge(iRow);
    rRowUB = vYEdge(iRow+1);
    
    vColFound = vCol((vRow > rRowLB) & (vRow <= rRowUB));
    
    if (~isempty(vColFound))
        
        
        vFound = histc (vColFound, vXEdge);
        
        nFound = (length(vFound)-1);
        
        if (nFound ~= nCol)
            disp([nFound nCol])
            error ('hist2d error: Size Error')
        end
        
        [nRowFound, nColFound] = size (vFound);
        
        nRowFound = nRowFound - 1;
        nColFound = nColFound - 1;
        
        if nRowFound == nCol
            mHist(iRow, :)= vFound(1:nFound)';
        elseif nColFound == nCol
            mHist(iRow, :)= vFound(1:nFound);
        else
            error ('hist2d error: Size Error')
        end
    end
    
end


