%First and last eigenvectors

brainFolder=savePlotsFolder+"/brainPlots/";
brainSubFolders={brainFolder+"eigVecs/",... 
                brainFolder+"lowEigVec/",...
                brainFolder+"highEigVec/",...
                brainFolder+"lowSignal/",...
                brainFolder+"highSignal/"};

if ~exist(brainFolder, 'dir')
   mkdir(brainFolder)
   mkdir(brainSubFolders{1})
   mkdir(brainSubFolders{2})
   mkdir(brainSubFolders{3})
   mkdir(brainSubFolders{4})
   mkdir(brainSubFolders{5})
end


disp("Plotting Single Group Mean Eigvecs")
vecs=[1,2,3,4,nROI-3,nROI-2,nROI-1,nROI];
CA=[0,max(mean(abs(U(:,vecs,:)),3),[],"all")];
for v=vecs
    disp("Vector: "+v)
    f=figure;
    f.Position(3:4)=[700 600];
    sgtitle("Mean abs. eigVec"+v,'horizontalAlignment','right');
    meanPD1=mean(squeeze(abs(U(:,v,indsPDHC{1}))),2);
    meanPD2=mean(squeeze(abs(U(:,v,indsPDHC{2}))),2);
    meanHC1=mean(squeeze(abs(U(:,v,indsPDHC{3}))),2);
    meanHC2=mean(squeeze(abs(U(:,v,indsPDHC{4}))),2);
    plot_on_brain([meanPD1,meanPD2,meanHC1,meanHC2],0,CodeBook,CA, ...
        ["PD1","PD2","HC1","HC2"],0);%always 0 threshold because we want to 
    % compare eigenvectors within subjects and between. Setting a thr
    % will make this difficut because for low eig vecs everything will be
    % close to 0.
    if(boolSavePlots)
        saveas(f,brainSubFolders{1}+"eigVec"+v+".png")
    end
    clf(f);
end
disp("Single Group Mean Eigvecs Done")

dataStrs=["Low","High"];


% for plotting with and without threshold
for boolThr=[0,1]
    disp("BoolThr: "+boolThr)
    disp("Plotting eigvecs")
    vecData={meanLowU,meanHighU};
    signalData={meanLowX,meanHighX};

    for i=1:length(vecData)
        data=vecData{i};
        datastr=dataStrs(i);

        f=figure;
        f.Position(3:4)=[700 600];
        sgtitle("Mean mag. "+datastr+" freq. eigvecs",'horizontalAlignment','right');
        plot_on_brain(data,0,CodeBook,0,["PD1","PD2","HC1","HC2"],boolThr)
        if(boolSavePlots)
            saveas(f,brainSubFolders{i+1}+"mean"+datastr+"EigVecs_thr"+boolThr+".png")
        end
        clf(f);

        % Only plot diff in eigvecs if we have FC and not SC
        % because SC will be same for diff PD and HC sessions
        if(boolFC) 

            f=figure;
            sgtitle("Abs Diff PD "+datastr+" eigenvectors")
            plot_on_brain([data(:,1),data(:,2),abs(data(:,1)-data(:,2))], ...
                0,CodeBook,0,["PD1","PD2","abs diff"],boolThr)
            if(boolSavePlots)
                saveas(f,brainSubFolders{i+1}+"mean"+datastr+"EigVecDiffPD_thr"+boolThr+".png")
            end
            clf(f);

            f=figure;
            sgtitle("Abs Diff HC "+datastr+" eigenvectors")
            plot_on_brain([data(:,3),data(:,4),abs(data(:,3)-data(:,4))], ...
                0,CodeBook,0,["HC1","HC2","abs diff"],boolThr)
            if(boolSavePlots)
                saveas(f,brainSubFolders{i+1}+"mean"+datastr+"EigVecDiffHC_thr"+boolThr+".png")
            end
            clf(f);
        end
        disp("Done " + i)
    end
    
    disp("Plotting signal")
    for i=1:length(signalData)
        data=signalData{i};
        datastr=dataStrs(i);
    
        f=figure;
        f.Position(3:4)=[700 600];
        sgtitle("Mean abs. "+datastr+" freq. signal",'horizontalAlignment','right');
        plot_on_brain(data,0,CodeBook,0,["PD1","PD2","HC1","HC2"],boolThr)
        if(boolSavePlots)
            saveas(f,brainSubFolders{i+3}+"mean"+datastr+"X_thr"+boolThr+".png")
        end
        clf(f);
    
        f=figure;
        sgtitle("Abs Diff PD mean "+datastr+" signal");
        plot_on_brain([data(:,1),data(:,2),abs(data(:,1)-data(:,2))], ...
            0,CodeBook,0,["PD1","PD2","abs diff"],boolThr)
        if(boolSavePlots)
            saveas(f,brainSubFolders{i+3}+"mean"+datastr+"XDiffPD_thr"+boolThr+".png");
        end
        clf(f);
    
        f=figure;
        sgtitle("Abs Diff HC mean "+datastr+" signal");
        plot_on_brain([data(:,3),data(:,4),abs(data(:,3)-data(:,4))], ...
            0,CodeBook,0,["HC1","HC2","abs diff"],boolThr)
        if(boolSavePlots)
            saveas(f,brainSubFolders{i+3}+"mean"+datastr+"XDiffHC_thr"+boolThr+".png");
        end
        clf(f);
        disp("Done " + i)
    end
    

end



% % % %%% Mean Eigenvectors %%%
% disp("Plotting low U")
% f=figure;
% f.Position(3:4)=[700 600];
% sgtitle("Mean mag. low freq. eigvecs",'horizontalAlignment','right');
% plot_on_brain(meanLowU,0,CodeBook,0,["PD1","PD2","HC1","HC2"],boolTHR)
% if(boolSavePlots)
%     saveas(f,brainSubFolders{2}+"meanLowEigVecs.png")
% end
% clf(f);
% 
% disp("Plotting high U")
% f=figure;
% f.Position(3:4)=[700 600];
% sgtitle("Mean mag. high freq. eigvecs",'horizontalAlignment','right');
% plot_on_brain(meanHighU,0,CodeBook,0,["PD1","PD2","HC1","HC2"],boolTHR)
% if(boolSavePlots)
%     saveas(f,brainSubFolders{3}+"meanHighEigVecs.png")
% end
% clf(f);
% 
% % There will only be a difference in the Laplacian Matrix if we use FC
% if(boolFC)
% 
%     disp("Plotting diff PD low eigvecs")
%     f=figure;
%     sgtitle("Abs Diff PD low eigenvectors")
%     plot_on_brain([meanLowU(:,1),meanLowU(:,2),abs(meanLowU(:,1)-meanLowU(:,2))], ...
%         0,CodeBook,0,["PD1","PD2","abs diff"],boolTHR)
%     if(boolSavePlots)
%         saveas(f,brainSubFolders{2}+"meanEigVecDiffLowPD.png")
%     end
%     clf(f)
% 
%     disp("Plotting diff HC low eigvecs")
%     f=figure;
%     sgtitle("Diff HC ses1 and HC ses 2 mean abs low eigenvectors")
%     plot_on_brain([meanLowU(:,3),meanLowU(:,4),abs(meanLowU(:,3)-meanLowU(:,4))], ...
%         0,CodeBook,0,["HC1","HC2","abs diff"],boolTHR)
%     if(boolSavePlots)
%         saveas(f,brainSubFolders{2}+"meanEigVecDiffLowHC.png")
%     end
%     clf(f);
% 
%     disp("Plotting diff PD high eigvecs")
%     f=figure;
%     sgtitle("Abs Diff PD high eigenvectors")
%     plot_on_brain([meanHighU(:,1),meanHighU(:,2),abs(meanHighU(:,1)-meanHighU(:,2))], ...
%         0,CodeBook,0,["PD1","PD2","abs diff"],boolTHR)
%     if(boolSavePlots)
%         saveas(f,brainSubFolders{3}+"meanEigVecDiffHighPD.png")
%     end
%     clf(f);
% 
%     disp("Plotting diff HC high eigvecs")
%     f=figure;
%     sgtitle("Diff HC ses1 and HC ses 2 mean abs high eigenvectors")
%     plot_on_brain([meanHighU(:,3),meanHighU(:,4),abs(meanHighU(:,3)-meanHighU(:,4))], ...
%         0,CodeBook,0,["HC1","HC2","abs diff"],boolTHR)
%     if(boolSavePlots)
%         saveas(f,brainSubFolders{3}+"meanEigVecDiffHighHC.png")
%     end
%     clf(f);
% end









% %%% Low signal plots %%%
% disp("Plotting low X")
% f=figure;
% f.Position(3:4)=[700 600];
% sgtitle("Mean abs. low freq. signal",'horizontalAlignment','right');
% plot_on_brain(meanLowX,0,CodeBook,0,["PD1","PD2","HC1","HC2"],boolTHR)
% if(boolSavePlots)
%     saveas(f,brainSubFolders{4}+"meanLowX.png")
% end
% clf(f);
% 
% f=figure;
% sgtitle("Abs Diff PD mean low signal");
% plot_on_brain([meanLowX(:,1),meanLowX(:,2),abs(meanLowX(:,1)-meanLowX(:,2))], ...
%     0,CodeBook,0,["PD1","PD2","abs diff"],boolTHR)
% if(boolSavePlots)
%     saveas(f,brainSubFolders{4}+"meanXDiffLowPD.png");
% end
% clf(f);
% 
% f=figure;
% sgtitle("Abs Diff HC mean low signal");
% plot_on_brain([meanLowX(:,3),meanLowX(:,4),abs(meanLowX(:,3)-meanLowX(:,4))], ...
%     0,CodeBook,0,["HC1","HC2","abs diff"],boolTHR)
% if(boolSavePlots)
%     saveas(f,brainSubFolders{4}+"meanXDiffLowHC.png");
% end
% clf(f);
% 
% %%% High signal plots %%%
% disp("Plotting high X")
% f=figure;
% f.Position(3:4)=[700 600];
% sgtitle("Mean abs. high freq. signal",'horizontalAlignment','right');
% plot_on_brain(meanHighX,0,CodeBook,0,["PD1","PD2","HC1","HC2"],boolTHR)
% if(boolSavePlots)
%     saveas(f,brainSubFolders{5}+"meanHighX.png")
% end
% clf(f);
% 
% f=figure;
% sgtitle("Abs Diff PD mean high signal");
% plot_on_brain([meanHighX(:,1),meanHighX(:,2),abs(meanHighX(:,1)-meanHighX(:,2))], ...
%     0,CodeBook,0,["PD1","PD2","abs diff"],boolTHR);
% if(boolSavePlots)
%     saveas(f,brainSubFolders{5}+"meanXDiffHighPD.png");
% end
% clf(f);
% 
% f=figure;
% sgtitle("Abs Diff HC mean high signal");
% plot_on_brain([meanHighX(:,3),meanHighX(:,4),abs(meanHighX(:,3)-meanHighX(:,4))], ...
%     0,CodeBook,0,["HC1","HC2","abs diff"],boolTHR);
% if(boolSavePlots)
%     saveas(f,brainSubFolders{5}+"meanXDiffHighHC.png");
% end
% clf(f);
% 
% 