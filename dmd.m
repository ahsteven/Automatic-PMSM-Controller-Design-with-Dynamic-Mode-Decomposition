function [ DMD_eigs, DMD_modes,Admd,S,x0] = dmd2(X,method,ndelays,r)
    % this code differs from Exact DMD Tu et al. presented in Kutz book because
    % the SVD is taken before the time shift. 

    %%%%%%%
    % This code is from Scott Dawson February 2020 
    % Script to implement various versions of DMD
    %
    % Inputs
    % Data: Snapshots of data, assumed to have uniform timestep
    % Method: string specifying method to use. Must be one of:
    %           'DMD': standard DMD
    %           'tls' total least-squares DMD
    %           'fb' forward-backward DMD
    % r: Truncaton level for SVD of Data
    %
    % Outputs
    % Admd: DMD propagation matrix
    % DMD_eigs: (discrete time) DMD eigenvalues
    % DMD_modes: DMD modes
    %
    % For more information, see Dawson, Hemati, Williams & Rowley,
    % "Characterizing and correcting for the effect of sensor noise in the
    % dynamic mode decomposition", Experiments in Fluids, 2016.
    %%%%%%%
    nt = length(X(1,:));

    % stack delay embeddings
    Data = zeros((ndelays+1),nt-ndelays); % noisy

    current = NaN((ndelays+1)*1,1);

    for j = 1:nt
        current = [X(:,j);current(1:end-1)];
        if j>ndelays
            Data(:,j-ndelays)=current;
        end
    end
    
    x0 = Data(:,1);
    m = size(Data,1);
    n = size(Data,2);
    if nargin < 3
        r = min(m,n);
        if nargin  < 2
            method = 'dmd';
        end   
    end

    [U,S,V] = svd(Data,'econ');

    Ur = U(:,1:r);
    Sr = S(1:r,1:r);
    Vr = V(:,1:r);

    % Project data onto first r POD modes
    DataProjected = Ur'*Data;

    Z1 = DataProjected(:,1:end-1);
    Z2 = DataProjected(:,2:end);

    if strcmp(method,'dmd')
        Admd = Z2/Z1;

    elseif strcmp(method,'tls')
        [U,~,~] = svd([Z1;Z2],'econ');
        U11 = U(1:r,1:r);
        U21 = U((r+1):end,1:r);
        if rank(U11)<r
            error('TLS solution does not exist')
        end
        Admd  = U21/U11; %note that this is still in projected space

    elseif strcmp(method,'fb')
        Af = Z2/Z1;
        Ab = Z1/Z2;
        Admd = (Af/Ab)^0.5;
%         Admd = (Af*pinv(Ab))^0.5;
    else
        disp('Invalid method specified')
    end

    [Evecs, DMD_eigs] = eig(Admd,'vector');
    DMD_modes = Ur*Evecs;

end