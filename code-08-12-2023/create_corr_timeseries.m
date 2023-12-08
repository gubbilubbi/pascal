function create_corr_timeseries(data,connType,indsPDHC,saveFolder)
    % Calculates a new timeseries from the existing one
    % that at every time point are the values of two regions
    % multiplied and then returns the mean and std of these new timeseries
    % data: matrix of size nTimepoints x nROI

    nTimes=size(data,1);
    nROI=size(data,2);
    nSubjs=size(data,3);

    corrMeans=zeros(nROI,nROI,nSubjs);
    corrStds=zeros(nROI,nROI,nSubjs);

    switch connType
        case 'corr'
            for s=1:nSubjs
                [corrMeans(:,:,s),corrStds(:,:,s)]=create_corr_timeseries(subData);
            end
        case 'Avgcorr'
            for i=1:length(indsPDHC)
                groupData=data(:,:,indsPDHC{i});

                % Stack all of the timeseries in each ROI (column) 
                % for each group
                groupROIData=reshape(permute(groupData,[1,3,2]),size(groupData,1)*size(groupData,3),size(groupData,2));
                [groupCorrMeans,groupCorrStds]=create_corr_timeseries(groupROIData);
                for s=indsPDHC{i}
                    corrMeans(:,:,s)=groupCorrMeans;
                    corrStds(:,:,s)=groupCorrStds;
                end
            end
    end
        
    corrTimeseries=zeros(nTimes,nROI,nROI);

    for regi=1:nROI
        for regj=regi+1:nROI
            corrTimeseries(:,regi,regj)=data(:,regi).*data(:,regj);
        end
    end
    corrMeans=squeeze(mean(corrTimeseries,1));
    corrStds=squeeze(std(corrTimeseries,0,1));

    save(saveFolder+"/meansCorr.mat","corrMeans")
    save(saveFolder+"/stdsCorr.mat","corrStds")
end