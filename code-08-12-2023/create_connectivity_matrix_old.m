function [FC,A]=create_connectivity_matrix(data,connType,indsPDHC,boolGroupFC, ...
                                       FCtoAdj)
    % data: time x nROI x nSubjs 
    % connType: String with text for type of FC to be calculated

    nROI=size(data,2);
    nSubjs=size(data,3);

    FC=zeros(nROI,nROI,nSubjs);
    switch connType
        case 'corr'
            % Correlation
            for s=1:nSubjs
                corrTemp=corr(data(:,:,s));
                FC(:,:,s)=corrTemp-diag(diag(corrTemp)); %set diag to 0
            end
            clear corrTemp;
        otherwise
    
    end


    % All subject in a group receive the same group mean connectivity matrix 
    if(boolGroupFC)
        for g=1:length(indsPDHC) % for each group
              groupFC=mean(FC(:,:,indsPDHC{g}),3);
              for sub=indsPDHC{g}
                    FC(:,:,sub)=groupFC;
              end
        end
    end

    
    A=FC;
    %If FC !>= 0 then the adjacency matrix is the FC but w_ij >= 0
    switch FCtoAdj
        case 'abs'
            A=abs(FC);
        case 'geq0' % All values below 0 set to 0
            A(FC<0)=0;
        otherwise
            error('No FC to Adj conversion')
    end
end