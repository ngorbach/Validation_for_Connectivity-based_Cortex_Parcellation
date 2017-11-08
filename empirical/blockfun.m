function X=blockfun(X,S,fun)
%BLOCKFUN applies a function on blocks of an array
% Y=BLOCKFUN(X,S,funHandle) applies the function funHandle on blocks of size 
% S on the array X.  Y is of size ceil(size(X)./S)
% 
% For instance, if X is of size [9 9] and S is [3 3], then 
% X is partitionned in 9 [3 3] blocks as follow :   
%  | B1 | B4 | B7 |
%  | B2 | B5 | B8 | where Bi's are [3 by 3] matrices
%  | B3 | B6 | B9 | 
% Y is then a [3 by 3] matrix with 
%  Y(i) = fun(Bi(:))
% 
% This function is just a simple wrap for accumarray !
% See Also : accumarray
%
% EXAMPLES  :
% rand('twister',12);
% x=floor(2*rand(9,9))
%   x =
%      0     0     1     1     1     0     0     0     0
%      1     0     0     0     1     0     0     1     1
%      0     1     1     1     0     0     0     0     0
%      1     1     0     1     0     0     1     1     1
%      0     1     0     0     0     0     0     0     0
%      1     0     0     0     1     0     0     1     1
%      1     1     0     1     0     1     1     1     0
%      0     1     1     1     1     1     1     1     1
%      1     0     0     1     0     0     0     0     1
%
% y=blockfun(x,[3 3],@sum)
% y =
%      4     4     2
%      4     2     5
%      5     6     6
% 
% y=blockfun(x,[3 3],@median)
% y =
%      0     0     0
%      0     0     1
%      1     1     1
%
% y=blockfun(x,[5 5],@sum)
% y =
%     12     5
%     11    10
% 
% x=floor(2*rand(6,6,6))
% y=blockfun(x,[3 2 3],@(x) length(find(x>0))>4)
%
% x=floor(2*rand(512,512,50));
% tic; y=blockfun(x,[5 5 5],@sum);  toc

assert(nargin<4,'Too many arguments');
assert((exist('X','var') & ~isempty(X)),'X is empty or is missing');
assert((exist('S','var') & ~isempty(S)),'S is empty or is missing');
assert((exist('fun','var') & ~isempty(fun)),'F is empty or is missing');
assert(isvector(S),'S should be a vector');
assert(isa(fun,'function_handle'),'fun should be a handle');

assert(ndims(X)==length(S),'S should have ndims(X) elements');
assert(numel(X)<=intmax('uint64'),'X is too big');


% Precompute stuff we will reuse
sizeX       = size(X);
nDim        = ndims(X);

% resulting size 
sizeY = ceil(sizeX./S);

% number of blocks
nBlock      = prod(sizeY);

% select adapted class for element adressing
classI = uintSelect(nBlock);

% Compute block index for each element of X
I_block = cell(1,nDim);

% sub indexes
for iDim = 1:nDim
    I_block{iDim} = classI(ceil((1:sizeX(iDim))/S(iDim)));
end  

I = classI(reshape(1:nBlock,sizeY));
I = I(I_block{:});

% And now perform the accumulation
I = colshape(I);
X = colshape(X);

X = accumarray(I,X,[nBlock 1],fun);
X = reshape(X,sizeY);


function x=colshape(x)
%COLSHAPE returns a column from an array
% X=COLSHAPE(X) is the same than x=x(:);

x=x(:);

function C=uintSelect(I)
%uintSelect returns the uint type adapted to a specific adressing 
% C=uintSelect(I) returns a casting function for selecting the type of
% uinteger depending on the max I to desribe.
%  
% Examples : 
%  uintSelect(255)
%   ans = 
%     @uint8
%  uintSelect(256)
%   ans = 
%     @uint16
%  uintSelect(2^16-1)
%   ans = 
%     @uint16
%  uintSelect(2^16)
%   ans = 
%     @uint32
    
list = {'uint64','uint32','uint16','uint8'};

N       = cell2mat(cellfun(@(x) uint64(intmax(x)),list,'UniformOutput',false));
iClass  = find(I<=N,1,'last');

assert(~isempty(iClass),'You are asking for too much');

C = str2func(list{iClass});
