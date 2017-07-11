
% select K patches from blocks with the most unsimilarity each other

function D = selectDic ( blocks, K )

D = zeros ( size(blocks,1), K);

engeryBlocks = sum ( blocks );
index = find (~engeryBlocks);
blocks(:,index) = [];
sim = blocks' * blocks;
minSim = min (sim);
[div, idx] = sort ( minSim );

D(:,1:K) = blocks(:,idx(1:K)); 
