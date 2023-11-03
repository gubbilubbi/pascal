%% 0. Settings
clear all

rng("default")

% set data path
mypath='/Users/Mortimtosh/Desktop/kth-MSc/data';

fs=200; % Sampling freq.

% The row indicies of the unordered data
order=[19,119,11,111,3,103,22,122,10,110,2,102,21,121,17,117,12,112,...
    14,114,13,113,5,105,20,120,7,107,18,118,9,109,8,108,1,101,6,106,...
    16,116,15,115,4,104];

%Order so that first 22 rows is left hemisphere, last 22 right hemisphere
[orderVals, orderInds]=sort(order); 

% visual: blue, sensorimotor: red, auditory: yellow, Temporal: black, 
% Posterior: cyan, anterior: magenta
colors=["b","r","y","k","c","m"];

% Each ROI belongs to a region
hemisphereRegions=[repmat("b",5,1);repmat("r",4,1);repmat("y",3,1);...
               repmat("k",2,1);repmat("c",4,1);repmat("m",4,1)];
cortexRegions=[hemisphereRegions;hemisphereRegions];

normMethod="zscore"; %'eucNorm', 'zscore', 'none'
boolCoh=0; % Calculate Coherence 
FCtoAdj="geq0"; % Value to determine FC to Adjacency conversion method

% Start and end of data timepoints to use
startInd=1;
endInd=84000;%94000;

%Strings for plotting
plotLegends=["PD ses1","PD ses2","HC ses1","HC ses2"];
subplotsInds={1,3,2,4};

% 1. Loading required data
L=endInd-startInd+1;


%IDs for PD and HC subjects
%idsPD=[1,2,3,5,7,8,9,10,11,13,14,15,21,28,31,33,37];
%idsHC=[4,6,12,16,17,18,19,20,22,23,24,25,26,27,29,30,32,34,35,36]; 
idsPD=[1,3];
idsHC=[2,4];
nPD=length(idsPD);
nHC=length(idsHC);

nROI=44;

data=load_data('scilife/meg/',nROI,startInd:endInd,idsHC);% time x regions x subjects
% Reorder data so in ascending order based on ID
data=data(:,orderInds,:);



nSubjs=size(data,3);
nROI=size(data,2);

startSes1HC=nPD*2+1;
endSes1HC=nPD*2+nHC;
startSes2HC=endSes1HC+1;
endSes2HC=endSes1HC+nHC;
% Start and end indicies for PD and HC in the dataset for ses1 and ses2
indsPDHC={1:nPD,nPD+1:nPD*2,startSes1HC:endSes1HC,startSes2HC:endSes2HC};
nGroups=length(indsPDHC); % 4 different groups in total

% Normalise data
switch normMethod

    case 'eucNorm' % Norm each region by its Euclidean norm 
        eNorm=vecnorm(data,2,1);
        data=data./eNorm;
        disp("Euclidean normalisation")

    case 'zscore' % Normalize each column (region) to mean 0 & std 1
        data=zscore(data); 
        disp("zscore normalisation")

    case 'none'
        disp("No normalisation")
    otherwise
        error(['Invalid Normalisation method: ' normMethod])
end

% Correlation
FC=zeros(nROI,nROI,nSubjs);
for s=1:nSubjs
    corrTemp=corr(data(:,:,s));
    FC(:,:,s)=corrTemp-diag(diag(corrTemp));
end

%The adjacency matrix is the FC but w_ij >= 0
A=FC;
switch FCtoAdj
    case 'abs'
        A=abs(FC);
    case 'geq0' % All values below 0 set to 0
        A(FC<0)=0;
    otherwise
        error(['Invalid FC to Adj value: ' FCtoAdj])
end

dirName="./scilife/data/Norm_"+normMethod+"_FCtoAdj_"+FCtoAdj;

if ~exist(dirName,'dir')
   mkdir(dirName)
end

save(dirName+"/dataNorm.mat","data")
save(dirName+"/A.mat","A");

