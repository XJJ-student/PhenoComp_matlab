function [populationDEG]=PhenoComp(expnormal,stablepair_freq,expE,individualfdr,method,populationfdr)
%stable_pairs
Allpair=[];
expN=cell2mat(expnormal(:,1));
gid=expN(:,1);
exp=expN(:,2:end);
[stable_pairs]=stablepair(gid,exp,stablepair_freq);
Allpair=stable_pairs;
if size(expnormal,2)==1
Allpair=Allpair;   
end   
if size(expnormal,2)>1
for j=2:size(expnormal,2)
	expN=cell2mat(expnormal(:,j));
	gid=expN(:,1);
    exp=expN(:,2:end);
	[stable_pairs]=stablepair(gid,exp,stablepair_freq);
	Allpair=intersect(Allpair,stable_pairs,'rows');
end
end
%individual_DEGs
stable_pair1=Allpair;
gidP1=unique([stable_pair1(:,1);stable_pair1(:,2)]);
gid1=expE(:,1); 
con_gid=intersect(gidP1,gid1);
[index,~]=ismember(gid1,con_gid);
expE=expE(index,:);
gid=expE(:,1); 
[result,~,~]=individual_def_gene(Allpair,expE,individualfdr);
L=[];
H=[];
for j=1:(size(expE,2)-1)
    a1=cell2mat(result(j,2));
    [Lia,~] = ismember(gid,a1);
    H(:,j)=Lia;
    b1=cell2mat(result(j,3));
    [Lia,~] = ismember(gid,b1);
    L(:,j)=Lia;
end
clear i j
%population_DEGs
[up_DEG,down_DEG]=ind2pop(H,L,gid,method,populationfdr);
populationDEG={up_DEG,down_DEG};
end










