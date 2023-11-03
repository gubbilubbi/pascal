%% More Functions


function data=load_data(mypath,nROI,times,idsHC)
    
    % mypath: 'scilife/meg/';
    % times: range of start to end indx for the data

    %mypath='scilife/meg/';
    %nROI=44;
    %times=1:5;
    %idsHC=[2,4]; % Test example, %Which of subjects 1-37 are HC
    % Real case
    %idsPD=[1,2,3,5,7,8,9,10,11,13,14,15,21,28,31,33,37];
    %idsHC=[4,6,12,16,17,18,19,20,22,23,24,25,26,27,29,30,32,34,35,36]; 

    % subjects in data are ordered as:
    % PD first sessions, PD second sessions
    % HC first sessions, HC second sessions

    % data: Matrix time x regions x subject 

    dataFolders=dir([mypath '0*']); % Only consider folders that start with 0
    idx = cellfun(@(x)contains(x,'0354'),{dataFolders.name}, 'Unif',0);  % EXclude 0354 because no .npy        
    dnew = dataFolders(~[idx{:}]);                                         
    dataFolderNames = {dnew.name}; 

    nSubjs=length(dataFolderNames);
    boolPD=create_bool_PD(nSubjs,idsHC);

    nPD=sum(boolPD);
    nHC=sum(~boolPD);

    indPD=1;
    indHC=1;

    data=zeros(length(times),nROI,nSubjs*2);
    pdData=zeros(length(times),nROI,nPD*2);
    hcData=zeros(length(times),nROI,nHC*2);

    %vvv=load('bool.mat');
    bool=cell(nSubjs,2);
    strs=["-ses1-","-ses2-"];
    for i=1:length(dataFolderNames)
        for j=1:length(strs)
            dataName=dataFolderNames{i}+strs(j)+"labeltc.npy";
            subj=readNPY(dataName)'; %time x regions    
            %subj=vvv.bool{i,j};
            % if PD
            if(boolPD(i))
                % if ses 1
                if(mod(j,2))
                    pdData(:,:,indPD)=subj(times,:);
                % if ses 2
                else
                    pdData(:,:,indPD+nPD)=subj(times,:);
                    indPD=indPD+1;
                end
            % if HC
            else
                % if ses 1
                if(mod(j,2))
                    hcData(:,:,indHC)=subj(times,:);
                % if ses 2
                else
                    hcData(:,:,indHC+nHC)=subj(times,:);
                    indHC=indHC+1;
                end
            end
            bool{i,j}=subj;
        end
    end
    % Stack the data PDses1, PDses2, HCses1, HCses2
    data(:,:,1:2*nPD)=pdData;
    data(:,:,2*nPD+1:end)=hcData;

    %save('./scilife/data/misc/bool.mat','bool')
    disp("Data extracted")
end

function boolPD=create_bool_PD(nSubjs,idsHC)
    boolPD=ones(nSubjs,1); % PD indices set to 1
    % HC indicies set to 0
    for row=idsHC
        boolPD(row)=0;
    end
    boolPD=logical(boolPD);
end

