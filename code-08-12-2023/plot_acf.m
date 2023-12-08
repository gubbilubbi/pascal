function acfVals=plot_acf(meg)    

    nROI=size(meg,2);
    nSubjs=size(meg,3);

    maxTime=500;

    acfVals=zeros(maxTime,nSubjs);

    for s=1:length(nSubjs)

        for r=1:nROI

            acfVals(:,s)=autocorr(meg(:,r,s),maxTime);

        end
    end
end

function x = autocorr(A,maxLag)

[row,col] = size(A);
if (row ~= 1 && col ~= 1)
    error('The input should be a vector, not a matrix!');
end
if row == 1
    A = A';
end

N = length(A);

% If no maxLag is defined then set it to the length of the timeseries
if(~maxLag)
    maxLag=N;
end
x = zeros(maxLag,1);
x(1) = sum(A.*A);

for ii = 2:maxLag
    B = circshift(A,-(ii-1));
    B = B(1:(N-ii+1));
    x(ii) = sum(B.*(A(1:(N-ii+1))));
end
x = x/x(1);
end