function main(file,gamma,save_name)
% File should include multi-omics data as well as network structure 
% Example:file='Metabric_breast_hprd_all_two_n.mat'
% 
% Gamma is a parameter for across channels interactions, 
% suggested value is 0.1 
%
% save_name is the file name of calculated distance matrix

% Jiening Zhu, 06/19/2021
%% Data loading and Markov invariant meeasures calculation
load(file);
[n,l]=size(RNA_s);
if exist('methyl_s','var')	
    % TCGA has three kinds of omic data: gene expressioan, methylation and copy number alteration
	[W_ge,W_me,W_cn]  = invariant_3( Adj, gene_s, RNA_s, methyl_s, CNA_s);
    Node_weights=zeros(n,l,3);
    Node_weights(:,:,1)=W_ge./sum(W_ge);
    Node_weights(:,:,2)=W_me./sum(W_me);
    Node_weights(:,:,3)=W_cn./sum(W_cn); 
else
    % METABRIC has only two kinds of omic data: gene expressioan and copy number alteration
	[W_ge,~,W_cn]  = invariant_3( Adj, gene_s, RNA_s, RNA_s, CNA_s);
    Node_weights=zeros(n,l,2);
    Node_weights(:,:,1)=W_ge./sum(W_ge);
    Node_weights(:,:,2)=W_cn./sum(W_cn);     
end

%% Graph structures
% graph G
m=sum(sum(Adj))/2;
D1=zeros(n,m); 
D2=zeros(n,m); %D=Divergence matrix
count=1;
for i=1:n-1
    for j=i+1:n
        if Adj(i,j)&&i~=j
            D1(i,count)=1;
            D2(j,count)=1;
            count=count+1;
        end
    end
end
m=count-1;
D1=D1(:,1:m);
D2=D2(:,1:m);

% graph F
nc=size(Node_weights,3);
if nc==3
    F2     = [1 0 0;1 0 0]';
    F1     = [0 1 0;0 0 1]';
elseif nc==2
    F2     = [1 0]';
    F1     = [0 1]';
end
[c d]=size(F1);

%% Pairwise distances calculation
distance=zeros(l,l);
parfor i=1:l-1
    i
    Distance_matrix_row= zeros(1,l);
    for j=i+1:l
        rho0=reshape(Node_weights(:,i,:),[n,nc]);
        rho1=reshape(Node_weights(:,j,:),[n,nc]);
        drho=rho0-rho1;
        Distance_matrix_row(j)=dist_cvx(drho(:),kron(speye(c),D1-D2),kron(F1-F2,speye(n)),gamma,m,n,c,d);
    end
    distance(i,:)= Distance_matrix_row;
end

distance=distance+distance';

%% Save results
save(save_name, 'distance')
end

