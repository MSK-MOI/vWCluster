function  [W_ge,W_me,W_cn]  = invariant_3( Adj, gene_s, RNA_s, methyl_s, CNA_s )
% Calculate Markov chain invariant measure
% Jiening Zhu, 06/19/2021

methyl_F=1-methyl_s; 

n=size(RNA_s,2);
m=length(gene_s);
W_ge=zeros(size(RNA_s));
W_me=zeros(size(RNA_s));
W_cn=zeros(size(RNA_s));
for i=1:n
    R_sum=zeros(m,1);
    R=repmat(RNA_s(:,i),1,m);
    R_sum=sum(R.*Adj)'; 
    
    C_sum=zeros(m,1);
    C=repmat(CNA_s(:,i),1,m);
    C_sum=sum(C.*Adj)';   
    
    M_sum=zeros(m,1);
    M=repmat(methyl_F(:,i),1,m);
    M_sum=sum(M.*Adj)';

    W_ge(:,i)=RNA_s(:,i).*R_sum;
    W_me(:,i)=methyl_F(:,i).* M_sum;
    W_cn(:,i)= CNA_s(:,i).* C_sum;
    
end

end

