function [result,p_value,index]=individual_def_gene(stable_pair,expE,fdr_cutoff)
stable_pair1=stable_pair;
gidP1=unique([stable_pair1(:,1);stable_pair1(:,2)]);
gid1=expE(:,1); 
con_gid=intersect(gidP1,gid1);
[index,~]=ismember(stable_pair1,con_gid);
stable_pair=stable_pair1(index(:,1)&index(:,2),:);
gidP=unique([stable_pair(:,1);stable_pair(:,2)]);
[index,~]=ismember(expE(:,1),con_gid);
expE=expE(index,:);
gid=expE(:,1); 
exp_tumor=expE(:,2:end);
for i = 1 : length(gidP)
    length(gidP)-i
   loc=find(stable_pair(:,1)==gidP(i));
    pair1=stable_pair(loc,:);
    a=size(pair1,1);
    loc=find(stable_pair(:,2)==gidP(i));
    pair2=stable_pair(loc,[2,1]);
    b=size(pair2,1);
    [~,loc1]=ismember(pair1,gid);
    indexx1=exp_tumor(loc1(:,1),:)-exp_tumor(loc1(:,2),:)<0;
    x=sum(indexx1,1);
    [~,loc2]=ismember(pair2,gid);
    indexx2=exp_tumor(loc2(:,1),:)-exp_tumor(loc2(:,2),:)>0;
    y=sum(indexx2,1);
    c=a-x+y;
    d=b-y+x;
    clear pair1 pair2 indexx1 indexx2  loc  loc1  loc2
    for j = 1 : size(exp_tumor,2)
       table=[a,b;c(j),d(j)];
       [~,p,~]=fishertest(table);
       p_value(i,j)=p;
       index(i,j)=(x(j))<y(j);   %1 up; 0 down
    end
    
end
    for k =1 : size(exp_tumor,2)
        time2=size(exp_tumor,2)-k
        [def_all,def_up,def_down]=fdr_adjust_bh(gidP,p_value(:,k),index(:,k),fdr_cutoff);
        result(k,1:3)={def_all,def_up,def_down};
    end
end