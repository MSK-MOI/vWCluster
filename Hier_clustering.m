%% Hierarchical clustering
% load('./Data/TCGA_breast_hprd_all_three')
% load('./Data/distance_matrix1_01')
load('data/Metabric_breast_hprd_all_two_n')
load('data/Mdistance_matrix1_01_n')


D= linkage(distance, 'ward');


NumCluster=4;

C = cluster(D,'Maxclust',NumCluster);

 
figure(1)
leafOrder = optimalleaforder(D, distance);
%cases_s=Cases(leafOrder);

%C=C(leafOrder);
 S_C=C(leafOrder);

color = D(end-NumCluster+2,3)-eps;
H = dendrogram(D, 0,'ReOrder', leafOrder,'ColorThreshold', color);
set(H,'LineWidth',1.5)
h = gca();

%// Changing the colours
lineColours = cell2mat(get(H,'Color'));
colourList = unique(lineColours, 'rows');

myColours = [0,0,255;
             247,53,99;
             129,43,146;
             253,154,82;
            255,0,255]/255;
     
%// Replace each colour (colour by colour). Start from 2 because the first colour are the "unclustered" black lines             
for colour = 2:size(colourList,1)
    %// Find which lines match this colour
    idx = ismember(lineColours, colourList(colour,:), 'rows');
    %// Replace the colour for those lines
    lineColours(idx, :) = repmat(myColours(colour-1,:),sum(idx),1);
end
%// Apply the new colours to the chart's line objects (line by line)
for line = 1:size(H,1)
    set(H(line), 'Color', lineColours(line,:));
end
 
%% Kaplan-Meier
%[num,txt,raw]=xlsread('.\data\TCGA_breast_60.xlsx');
[num,txt,raw]=xlsread('.\data\Metabric_breast.xlsx');
%C=num(:,2);
num=num(:,1);
e=txt(2:end,3);
e = strrep(e,'''','');
[p,fh,stats]=MatSurv(num,e,C,'GroupsToUse',[1 2 3 4]);