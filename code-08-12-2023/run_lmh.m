
%% MAIN 2022

%%  Graph Signal Processing (GSP) for Neuroimaging - Giulia Preti 17.12.2022

% Author: Valter Lundegårdh

% Inspiration from Giulia Preti's Graph Signal Processing
% https://github.com/gpreti/GSP_StructuralDecouplingIndex/tree/master
% DOI: 10.1109/ICASSP49357.2023.10095285

%% 0. Settings
clear all

set(0,'DefaultFigureVisible','on')

rng("default")

fs=200; % Sampling freq.

% The row indicies of the unordered data
order=[19,119,11,111,3,103,22,122,10,110,2,102,21,121,17,117,12,112,...
    14,114,13,113,5,105,20,120,7,107,18,118,9,109,8,108,1,101,6,106,...
    16,116,15,115,4,104];
%Order so that first 22 rows is left hemisphere, last 22 right hemisphere
[orderVals, orderInds]=sort(order); 

% Visual: blue, Sensorimotor: red, Auditory: yellow, Temporal: green, 
% Posterior: cyan, Anterior: magenta
areaColors=["b";"r";"y";"g";"c";"m";"b";"r";"y";"g";"c";"m"];

% Each ROI belongs to a region
hemisphereColors=[repmat("b",5,1);repmat("r",4,1);repmat("y",3,1);...
               repmat("g",2,1);repmat("c",4,1);repmat("m",4,1)];
cortexColors=[hemisphereColors;hemisphereColors];

% Settings
methodFC="Avgcorr";
boolFDR=1; %Apply FDR on the FC matrix and the correlation values
boolPrec=0;

boolParCorr=0;
normMethod="zscore"; %'eucNorm','zscore','none'
FCtoAdj="abs"; % Value to determine FC to Adjacency conversion method
boolGroupFC=0; %Whether each group or subject has one FC

capOutliers=0; % To cap outliers in the data to a value (4) if zscored
boolFC=1; % If we have FC or SC

boolSavePlots=0;

%%% CHOICE OF CUTOFF
% choose here between 'fixed' (arbitrary threshold) 
% 'fixed' arbitrary selected cut-off
% 'median': median for each subject and for each group compared
% Should be the same though
cutoff_option='fixed';

% Start and end of data timepoints to use
startInd=1;
endInd=94000;

%IDs for PD and HC subjects
idsPD=[1,2,3,5,7,8,9,10,12,13,14,20,27,30,32,36];
idsHC=[4,6,11,15,16,17,18,19,21,22,23,24,25,26,28,29,31,33,34,35];
%idsPD=[1,3];
%idsHC=[2,4];
nPD=length(idsPD);
nHC=length(idsHC);

%Strings for plotting
plotLegends=["PD ses1","PD ses2","HC ses1","HC ses2"];
subplotsInds={1,3,2,4};

% 1. Loading required data
L=endInd-startInd+1;

% Length of cortexColors = nROI
%
data=load_data('/archive/21098_pd_dimensionality/MEG/',length(cortexColors),startInd:endInd,idsHC);% time x regions x subjects
%data=load_data('scilife/MEG/',length(cortexColors),startInd:endInd,idsHC);% time x regions x subjects

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

        if(capOutliers)
            %thr=3; %0.00135 => ~1% of data for every region
            thr=4; %0.00003 => ~0.3% of data for every region
            data(data>thr)=thr;
            data(data<-thr)=-thr;
        end

    case 'none'
        disp("No normalisation")
    
    otherwise
        error(['Invalid Normalisation method: ' normMethod])
end

%% ONLY RUN THIS ONCE PER methodFC
saveDotFolder="./pascal/data/dot-product";
%savePlotsFolder="./scilife/testData/"+folderName;
if ~exist(saveDotFolder, 'dir')
       mkdir(saveDotFolder)
end
create_corr_timeseries(data,methodFC,indsPDHC,saveDotFolder);

%%


[FC,A,pvalsMat,corrMeans,corrStds]=create_connectivity_matrix(data,methodFC,indsPDHC,boolGroupFC,FCtoAdj,boolFDR,boolParCorr,boolPrec);

% normalize adjacency matrix
DegMat=zeros(nROI,nROI,nSubjs); % Cleared afterwards
Anorm=zeros(nROI,nROI,nSubjs); % Cleared afterwards
Lnorm=zeros(nROI,nROI,nSubjs); % Cleared afterwards
U=zeros(nROI,nROI,nSubjs);
LambdaL=zeros(nROI,nSubjs);

for s=1:nSubjs
    DegMat(:,:,s)=diag(sum(A(:,:,s),2)); %% Degree
    %Symmetrically normalized Adjacency Matrix 
    Anorm(:,:,s)=DegMat(:,:,s)^(-1/2)*A(:,:,s)*DegMat(:,:,s)^(-1/2); 
    
    % STRUCTURAL CONNECTOME DECOMPOSITION INTO STRUCTURAL HARMONICS
    % Compute symmatrically normalized Laplacian
    Lnorm(:,:,s)=eye(nROI)-Anorm(:,:,s); 
    % Laplacian eigendecomposition (eigenvectors and values)
    [tempU,tempLambdaL] = eig(Lnorm(:,:,s));
    % eigenvalues sorted ascendingly = spatial frequencies
    [LambdaL(:,s), Ind]=sort(diag(tempLambdaL)); 
    % eigenvectors = structural harmonics 
    U(:,:,s)=tempU(:,Ind); 

    clear tempU; clear tempLambdaL; 
    clear Ind; clear D; clear Anorm; clear Lnorm;
end 

%%% 2. GRAPH FOURIER TRANSFORM
X_hat=zeros(nROI,L,nSubjs);
for s=1:nSubjs %graph fourier transform
    %nROIxnROI*nROIxtime
    % U Complex Transpose
    X_hat(:,:,s)=U(:,:,s)'*data(:,:,s).';
end

% Plot eigenvalues (i.e. variation of eigenvectors) and select cutoff
% plot_subjects_and_mean(LambdaL,"Laplacian eigenvalues/Variation of eigenvectors" ...
%     ,"Spectral index","Variation",nGroups,indsPDHC,subplotsInds,plotLegends, ...
%     NaN,0)

folderName=normMethod+"_FC="+methodFC+"_FCtoAdj="+FCtoAdj+"_boolFDR="+boolFDR+...
    "boolParCorr_+"+boolParCorr+"_boolPrec="+boolPrec+...
    "_capOutliers="+capOutliers+"_start="+startInd+"_end="+endInd;
savePlotsFolder="./pascal/data/"+folderName;
%savePlotsFolder="./scilife/testData/"+folderName;
if ~exist(savePlotsFolder, 'dir')
       mkdir(savePlotsFolder)
end

% Choose cutoffs
for low=4:22 %4:22p
    disp(low)
    for up=23:41 %23:41

        NN=[low,low,low,low;up,up,up,up];
        strCutoffs=strjoin(reshape((string([low;up])).',1,[]),"_");
                   cutoffFolder=savePlotsFolder+"/cutoff_"+cutoff_option+"_"+strCutoffs+"/";
                   mkdir(cutoffFolder)
        
        %%% Create 2 versions of the U matrix, containing only low- or high- frequency Laplacian eigenvectors
        Ulow=zeros(size(U));
        Umid=zeros(size(U));
        Uhigh=zeros(size(U));
        for g=1:nGroups
            Ulow(:,1:NN(1,g),indsPDHC{g})=U(:,1:NN(1,g),indsPDHC{g}); %low spatial frequencies  
            Umid(:,NN(1,g)+1:NN(2,g)-1,indsPDHC{g})=U(:,NN(1,g)+1:NN(2,g)-1,indsPDHC{g});
            Uhigh(:,NN(2,g):end,indsPDHC{g})=U(:,NN(2,g):end,indsPDHC{g}); %high spatial frequencies
        end
        
        %%% reconstruct functional signals containing only low / high frequencies
        
        % This both filters and inverts the data back to time domain 
        % because V^-1 * H^ (H binary diag matrix) takes only the specified values
        % from V^-1, i.e. low or high components.
        
        Xc=zeros(L,nROI,nSubjs);
        Xm=zeros(L,nROI,nSubjs);
        Xd=zeros(L,nROI,nSubjs);
        for s=1:nSubjs
            Xc(:,:,s)=(Ulow(:,:,s)*X_hat(:,:,s)).'; %reconstruction of coupled signal (only low frequencies)
            Xm(:,:,s)=(Umid(:,:,s)*X_hat(:,:,s)).';
            Xd(:,:,s)=(Uhigh(:,:,s)*X_hat(:,:,s)).'; %reconstruction of decoupled signal (only high frequencies)
        end
        
        % Energy (l2-norm) over time for each region and subject
        energyX=zeros(nROI,nSubjs,3);
        energyX(:,:,1)=squeeze(vecnorm(Xc,2,1));
        energyX(:,:,2)=squeeze(vecnorm(Xm,2,1));
        energyX(:,:,3)=squeeze(vecnorm(Xd,2,1));
        
        meanX=zeros(nROI,nGroups,3); %low,mid,high
        stdX=zeros(nROI,nGroups,3); %low,mid,high
        
        for g=1:nGroups
            % Mean and std of reconstructed signals over time and session group 
            meanX(:,g,1)=mean(abs(Xc(:,:,indsPDHC{g})),[1,3]);
            meanX(:,g,2)=mean(abs(Xm(:,:,indsPDHC{g})),[1,3]);
            meanX(:,g,3)=mean(abs(Xd(:,:,indsPDHC{g})),[1,3]);
            stdX(:,g,1)=std(abs(Xc(:,:,indsPDHC{g})),0,[1,3]);
            stdX(:,g,2)=std(abs(Xm(:,:,indsPDHC{g})),0,[1,3]);
            stdX(:,g,3)=std(abs(Xd(:,:,indsPDHC{g})),0,[1,3]);
        end
        
        save(cutoffFolder+"/meanX.mat",'meanX')
        save(cutoffFolder+"/stdX.mat",'stdX')
        save(cutoffFolder+"/energyX.mat",'energyX')
    
    end
end

% Save FC, pvals from correlation calc as well as means and stds from
% the dot-product timeseries
save(savePlotsFolder+"/FC.mat","FC")
save(savePlotsFolder+"/pvalsMat.mat","pvalsMat")

% params for exponential fitted to acf
acfVals=plot_acf(data);
save(savePlotsFolder+"/acfVals.mat","acfVals")


%% Show figures yes or no
show_figures=0;
if(~show_figures)
    set(0,'DefaultFigureVisible','off')
end
%% Plot on brain
%TODO: Uncomment if the codebook is not created yet
%create_codebook 

% Atlas + Selection
switch nROI
    case 44
        disp("44 regions")
        CodeBookpath=which('myCodeBook.mat');
    case 360
        disp("360 regions")
        CodeBookpath=which('Glasser360_2mm_codebook.mat');
    otherwise
        error(['No available atlas for the number of regions being: ' nROI])
end

CodeBook=load(CodeBookpath);
CodeBook=CodeBook.myCodeBook;

% Pre-allocate space
meanLowU=zeros(nROI,nGroups);meanMidU=zeros(nROI,nGroups);meanHighU=zeros(nROI,nGroups);
stdLowU=zeros(nROI,nGroups);stdMidU=zeros(nROI,nGroups);stdHighU=zeros(nROI,nGroups);
meanLowX=zeros(nROI,nGroups);meanMidX=zeros(nROI,nGroups);meanHighX=zeros(nROI,nGroups);
stdLowX=zeros(nROI,nGroups);stdMidX=zeros(nROI,nGroups);stdHighX=zeros(nROI,nGroups);

meanEigVecs=zeros(nROI,nROI,nGroups);


% Energy (l2-norm) over time for each region and subject
XcE=squeeze(vecnorm(Xc,2,1));
XmE=squeeze(vecnorm(Xm,2,1));
XdE=squeeze(vecnorm(Xd,2,1));
% Brain plots
for g=1:nGroups
    % Mean and std magnitude of eigvecs over freqs and sessions
    % meanLowU(:,g)=mean(abs(Ulow(:,1:NN(1,g),indsPDHC{g})),[2,3]);
    % meanMidU(:,g)=mean(abs(Uhigh(:,NN(1,g)+1:NN(2,g)-1,indsPDHC{g})),[2,3]);
    % meanHighU(:,g)=mean(abs(Uhigh(:,NN(2,g):end,indsPDHC{g})),[2,3]);
    % stdLowU(:,g)=std(abs(Ulow(:,1:NN(1,g),indsPDHC{g})),0,[2,3]);
    % stdMidU(:,g)=std(abs(Uhigh(:,NN(1,g)+1:NN(2,g)-1,indsPDHC{g})),0,[2,3]);
    % stdHighU(:,g)=std(abs(Uhigh(:,NN(2,g):end,indsPDHC{g})),0,[2,3]);

    % Mean and std of reconstructed signals over time and session group 
    meanLowX(:,g)=mean(abs(Xc(:,:,indsPDHC{g})),[1,3]);
    meanMidX(:,g)=mean(abs(Xm(:,:,indsPDHC{g})),[1,3]);
    meanHighX(:,g)=mean(abs(Xd(:,:,indsPDHC{g})),[1,3]);
    stdLowX(:,g)=std(abs(Xc(:,:,indsPDHC{g})),0,[1,3]);
    stdMidX(:,g)=std(abs(Xm(:,:,indsPDHC{g})),0,[1,3]);
    stdHighX(:,g)=std(abs(Xd(:,:,indsPDHC{g})),0,[1,3]);
end

%brain_plots

%% Plot timeseries

% Indicies of the data used
times=startInd:endInd;

% Only need the indicies for the first session for PD and HC
plot_timeseries(data,Xc,Xd,times,indsPDHC{1},nPD,"PD-patient-", ...
                savePlotsFolder,boolSavePlots)
plot_timeseries(data,Xc,Xd,times,indsPDHC{3},nHC,"HC-patient-", ...
                savePlotsFolder,boolSavePlots)

%% Plot energy concentration

% When the data has been EucNorm i.e. 1=||x||_2^2 =
% = sqrt(sum(x_l)_i^2+sum(x_m)_j^2+sum(x_h)_k^2)
% then we know everything under tahe sqrt = 1 => sum(x_l)_i^2 <= 1
% so it is sort of a measure of how much of the total signal that is 
% present in the lower frequencies. This holds for other normalisation
% schemes too but they might not equal to 1.

% This is not true, the sum of two added vectors can have a shorter length
% than the vectors e.g. [1,3]+[1,-2] or [-1,1]+[1,-2]

% Energy (l2-norm) over time for each region and subject
XcE=squeeze(vecnorm(Xc,2,1));
XdE=squeeze(vecnorm(Xd,2,1));

% Mean energy for each session
meanXcE=zeros(44,nGroups);
meanXdE=zeros(44,nGroups);
stdXcE=zeros(44,nGroups);
stdXdE=zeros(44,nGroups);


for g=1:nGroups
    meanXcE(:,g)=mean(XcE(:,indsPDHC{g}),2);
    meanXdE(:,g)=mean(XdE(:,indsPDHC{g}),2);
    stdXcE(:,g)=std(XcE(:,indsPDHC{g}),0,2);
    stdXdE(:,g)=std(XdE(:,indsPDHC{g}),0,2);
end

nROIPerArea=[5,4,3,2,4,4,5,4,3,2,4,4];

meanAreaXcE=zeros(length(nROIPerArea),nGroups);
meanAreaXdE=zeros(length(nROIPerArea),nGroups);
stdAreaXcE=zeros(length(nROIPerArea),nGroups);
stdAreaXdE=zeros(length(nROIPerArea),nGroups);

ind=0;
for i=1:length(nROIPerArea)
    meanAreaXcE(i,:)=mean(meanXcE(ind+1:nROIPerArea(i)+ind,:),1);
    meanAreaXdE(i,:)=mean(meanXdE(ind+1:nROIPerArea(i)+ind,:),1);
    stdAreaXcE(i,:)=std(meanXcE(ind+1:nROIPerArea(i)+ind,:),0,1);
    stdAreaXdE(i,:)=std(meanXdE(ind+1:nROIPerArea(i)+ind,:),0,1);
    ind=nROIPerArea(i);
end

%For areas (i.e. regions -> areas)
plot_energy_conc(meanAreaXcE,stdAreaXcE,areaColors, ...
    plotLegends,subplotsInds,"Low freq. average energy per area (Top half LH, bottom half RH)", ...
    savePlotsFolder+"/energyLfreqArea.png",boolSavePlots)
plot_energy_conc(meanAreaXdE,stdAreaXdE,areaColors, ...
    plotLegends,subplotsInds,"High freq. average energy per area (Top half LH, bottom half RH)", ...
    savePlotsFolder+"/energyHfreqArea.png",boolSavePlots)

%For regions 
plot_energy_conc(meanXcE,stdXcE,cortexColors, ...
    plotLegends,subplotsInds,"Low freq. average energy per region (Top half LH, bottom half RH)", ...
    savePlotsFolder+"/energyLfreqRegion.png",boolSavePlots)
plot_energy_conc(meanXdE,stdXdE,cortexColors, ...
    plotLegends,subplotsInds,"High freq. average energy per region (Top half LH, bottom half RH)", ...
    savePlotsFolder+"/energyHfreqRegion.png",boolSavePlots)


%%% Difficult to see things here... %%%
% % Sort energy in ascending order
% [sortedMeanXcE, indXc]=sort(meanXcE);
% [sortedMeanXdE, indXd]=sort(meanXdE);
% colorXc=cortexColors(indXc);
% colorXd=cortexColors(indXd);
% 
% % Sort left hemisphere in ascending order
% [sortedLHMeanXcE, indXcLH]=sort(meanXcE(1:nROI/2,:));
% [sortedLHMeanXdE, indXdLH]=sort(meanXdE(1:nROI/2,:));
% colorXcLH=hemisphereColors(indXcLH);
% colorXdLH=hemisphereColors(indXdLH);
% 
% % Sort right hemisphere in ascending order
% [sortedRHMeanXcE, indXcRH]=sort(meanXcE(nROI/2+1:end,:));
% [sortedRHMeanXdE, indXdRH]=sort(meanXdE(nROI/2+1:end,:));
% colorXcRH=hemisphereColors(indXcRH);
% colorXdRH=hemisphereColors(indXdRH);

% plot_energy_conc(sortedMeanXcE,colorXc,plotLegends,subplotsInds,"Low freq.")
% plot_energy_conc(sortedMeanXdE,colorXd,plotLegends,subplotsInds,"High freq.")
% 
% plot_energy_conc(sortedLHMeanXcE,colorXcLH,plotLegends,subplotsInds,"LH Low freq.")
% plot_energy_conc(sortedLHMeanXdE,colorXdLH,plotLegends,subplotsInds,"LH High freq.")
% 
% plot_energy_conc(sortedRHMeanXcE,colorXcRH,plotLegends,subplotsInds,"RH Low freq.")
% plot_energy_conc(sortedRHMeanXdE,colorXdRH,plotLegends,subplotsInds,"RH High freq.")

%% Visualization of the Connectivity Matrix's Mean and Variance %%
%%% PLOTTING %%%
stringMat='FC'; % Which connectivity matrix to plot
% Only for Adj as log(-FC) does not work
boolLog=0; % Log-plots (The log scale allows you to better see contrast)

% Connectivity matrix to be plotted
switch stringMat
    case 'Adj'
        cMat=A; 
        cStr=" Adj ";
    case 'FC'
        cMat=FC;
        cStr=" FC ";
    otherwise
end

maxMatVar=max(var(cMat,0,3),[],"all");
plotType=["Mean","Var"];
% Plot Mean connectivity 
for i=0:1
    f=figure;
    for g=1:nGroups
        if(~i)
            toPlot=mean(cMat(:,:,indsPDHC{g}),3);
        else
            toPlot=var(cMat(:,:,indsPDHC{g}),0,3);
        end
        logStr="";
        if(boolLog)
            logStr="Log ";
            meanPlot=log(meanPlot);
            varPlot=log(varPlot);
        end
        subplot(3,2,subplotsInds{g});imagesc(toPlot);
        title(logStr+plotType(i+1)+cStr+plotLegends(g));
        xlabel('regions');ylabel('regions');colorbar; 
        if(~i) %meanFC limited to [0,1]
            clim([0,1]);
        else %Var limited from 0 to maxVar + 5%
            clim([0,maxMatVar+maxMatVar*0.05]);
        end
    end
    if(~i)
        subplot(3,2,5);
        imagesc(abs(mean(cMat(:,:,indsPDHC{1}),3)-mean(cMat(:,:,indsPDHC{2}),3)));
        colorbar;
        clim([0,1])
        subplot(3,2,6);
        imagesc(abs(mean(cMat(:,:,indsPDHC{3}),3)-mean(cMat(:,:,indsPDHC{4}),3)));
        colorbar;
        clim([0,1])
    end
    if(boolSavePlots)
        saveas(f,savePlotsFolder+"/"+plotType(i+1)+"FC.png")
    end
end

%% Plot eigenvalues and structural harmonics and wZC

% Plot distribution of eigenvalues
eigFig=figure;
sgtitle("Eigenvalue distributions")
for g=1:nGroups
    subplot(2,2,subplotsInds{g})
    histfit(reshape(LambdaL(:,indsPDHC{g}),1,[]),44);
    title(plotLegends{g})
    xlabel("Eigenvalue");ylabel("Count")
    %ylim([0,])
end

saveas(eigFig,savePlotsFolder+"/eigDistr.png")

% Plot eigenvalues (i.e. variation of eigenvectors)
plot_subjects_and_mean(LambdaL,"Laplacian eigvals/Variation of eigvecs" ...
    ,"Spectral index","eigVal",nGroups,indsPDHC,subplotsInds,plotLegends, ...
    savePlotsFolder,boolSavePlots)

% Plot all of the eigenvectors values, mean over each session
meanEigVec=figure;
sgtitle('Mean Laplacian Eigvecs');
minU=min(U,[],"all");
maxU=max(U,[],"all");
for g=1:nGroups
    subplot(2,2,subplotsInds{g});imagesc(mean(U(:,:,indsPDHC{g}),3));
    title(plotLegends(g))
    xlabel('Spectral index');ylabel('regions')  
    colorbar;
    clim([minU,maxU]);
end

saveas(meanEigVec,savePlotsFolder+"/meanEigVec.png")

% Plot num zero-corssings
plot_subjects_and_mean(ZC,"Zero Crossings (ZC)" ...
    ,"Connectome harmonics","ZC",nGroups,indsPDHC,subplotsInds,plotLegends, ...
    savePlotsFolder,boolSavePlots)
% Plot weighted zero-crossings
plot_subjects_and_mean(wZC,"Weigthed Zero Crossings (wZC)" ...
    ,"Connectome harmonics","wZC",nGroups,indsPDHC,subplotsInds,plotLegends, ...
    savePlotsFolder,boolSavePlots)
%% Spectral X
%for each timepoint, you have one coefficient for each eigenvector, i.e., the contribution of each eigenvector to the functional signal

%figure;imagesc(X_hat(:,:,4));title('Spectral Coefficients - Subject 1');colorbar;xlabel('time');ylabel('spectral index')

% These are simply the values of the signal in the freq domain. It makes
% sense that we get +- large values for low spectral indices because those
% are the indices belonging to aligned eigenvectors => we would expect the energy to be
% the largest among these if the signal is overall fairly aligned with the
% network. But for the PD RS data we don't see this...

%% PLOTTING
for g=1:nGroups
    % If we want the same colorbar we can add ,[-0.2,0.2] after 
    % Uhigh in imagesc
    figure;
    sgtitle(["Mean high and low eigvecs " plotLegends{g}])
    subplot(3,1,1);
    imagesc(mean(Uhigh(:,:,indsPDHC{g}),3));colorbar;title('High Frequency Structural Harmonics');
    xlabel('spectral index');ylabel('regions')
    subplot(3,1,2);
    imagesc(mean(Umid(:,:,indsPDHC{g}),3));colorbar;title('Mid Frequency Structural Harmonics');
    xlabel('spectral index');ylabel('regions')
    subplot(3,1,3);
    imagesc(mean(Ulow(:,:,indsPDHC{g}),3));colorbar;title('Low Frequency Structural Harmonics');
    xlabel('spectral index');ylabel('regions')
end
%% 
function plot_timeseries(meg,coupled,decoupled,times,ids,nPats,titStr,folder,boolSave)
    % meg: 3d Matrix to plot  time x regions x subjects
    % coupled: coupled timeseries regions x time x subjects
    % decoupled: decoupled timeseries regions x time x subjects
    % times: plotStart:plotEnd range of ints Int
    % ids: ids for patients
    % nPats: number of PD or HC patients in total
    % titStr: string for the title ("PD id " or "HC id ")

    % Plots the MEG,coupled and decoupled timeseries for the given patients
    % left column ses1, right column ses 2.

    sess=["(ses 1)","(ses 2)"];
    SPinds={[1,3,5],[2,4,6]}; %indices subplot
    
    for i=1:nPats % Patients
        f = figure;
        f.Position = [100 100 1400 700];
        
        %Find min and max over sessions for the subject
        subj=ids(i);
        subjs=[subj,subj+nPats];
        megYmin=min(meg(times,:,subjs),[],"all");
        megYmax=max(meg(times,:,subjs),[],"all");
        cYmin=min(coupled(times,:,subjs),[],"all");
        cYmax=max(coupled(times,:,subjs),[],"all");
        dYmin=min(decoupled(times,:,subjs),[],"all");
        dYmax=max(decoupled(times,:,subjs),[],"all");
        
        % Column
        for j=1:2
            subplot(3,2,SPinds{j}(1));
            plot(meg(times,:,subj+(j-1)*nPats));
            title("MEG timecourses "+sess(j));xlabel('time');ylabel('amplitude')
            ylim([megYmin,megYmax]);

            subplot(3,2,SPinds{j}(2));
            plot(coupled(times,:,subj+(j-1)*nPats))
            title("Coupled signal timecourses "+sess(j));xlabel('time');ylabel('amplitude');
            ylim([megYmin,megYmax]);
            %ylim([cYmin,cYmax]);

            subplot(3,2,SPinds{j}(3));
            plot(decoupled(times,:,subj+(j-1)*nPats))
            title("Decoupled signal timecourses "+sess(j));xlabel('time');ylabel('amplitude');
            ylim([megYmin,megYmax]);
            %ylim([dYmin,dYmax]);
        end
        sgtitle(titStr+string(i))
        timeFolder=folder+"/time/";
        if ~exist(timeFolder, 'dir')
           mkdir(timeFolder)
        end
        if(boolSave)
            saveas(f,timeFolder+titStr+string(i)+"_timeseries.png")
        end
        clf(f)
    end
end

function plot_subjects_and_mean(toBePlotted,sgtit,xlab,ylab,nGroups,indsPDHC, ...
    subplotsInds,plotLegends,folder,savePlot)

indvFig=figure;
sgtitle("Individual "+sgtit);
% Plot individual subjects and group mean
meanG=zeros(size(toBePlotted,1),nGroups);
stdG=zeros(size(toBePlotted,1),nGroups);

%upper=meanG(:,g)+stdG(:,g);
%lower=meanG(:,g)-stdG(:,g);
yMax=max(toBePlotted,[],"all");
for g=1:nGroups 
    subplot(2,2,subplotsInds{g});plot(toBePlotted(:,indsPDHC{g}));
    hold on;
    meanG(:,g)=mean(toBePlotted(:,indsPDHC{g}),2);
    stdG(:,g)=std(toBePlotted(:,indsPDHC{g}),0,2); 
    subplot(2,2,subplotsInds{g});plot(meanG(:,g),'k');
    title(plotLegends(g));xlabel(xlab);ylabel(ylab);
    ylim([0 yMax+0.05*yMax]);
end
if(savePlot)
    saveas(indvFig,folder+"/indv_"+ylab+".png")
end

% Compare mean between groups
meanFig=figure;
for g=1:nGroups
    errorbar(1:size(meanG, 1), meanG(:,g), stdG(:,g))
    xlabel(xlab);ylabel(ylab);
    hold on;
end
title("Mean "+sgtit);
legend(plotLegends,'Location','southeast');

if(savePlot)
    saveas(meanFig,folder+"/mean_"+ylab+".png")
end

end

