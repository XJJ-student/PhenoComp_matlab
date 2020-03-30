function [up_DEG,down_DEG]=ind2pop(H,L,gene,method,populationfdr)
p0_result=[];      
for i=1:size(H,2)
  p_result1=(sum(H(:,i)==1)+sum(L(:,i)==1))/size(H,1);
  p0_result=[p0_result;p_result1];
end
pup_result=[];      
for i=1:size(H,2)
  p_result1=sum(H(:,i)==1)/(sum(H(:,i)==1)+sum(L(:,i)==1));
  pup_result=[pup_result;p_result1];
end
pdown_result=[];      
for i=1:size(L,2)
  p_result1=sum(L(:,i)==1)/(sum(H(:,i)==1)+sum(L(:,i)==1));
  pdown_result=[pdown_result;p_result1];
end
if method==1
   p0_up=median(p0_result)*median(pup_result);
   p0_down=median(p0_result)*median(pdown_result);
end
if method==2
   p0_up=mean(p0_result)*mean(pup_result);
   p0_down=mean(p0_result)*mean(pdown_result);
end
up_result=[];
for i = 1:size(H,1)
    k=sum(H(i,:)==1);
    s=size(H,2);
    up_result(i,2)=k/s;
    up_result(i,3)=1-binocdf(k-1,s,p0_up);
end      
up_result(:,1)=gene;
up_result(:,4)=mafdr(up_result(:,3),'BHFDR',true);
    
up_DEG=up_result(up_result(:,4)<populationfdr,1);
down_result=[];
for i = 1:size(L,1)
    k=sum(L(i,:)==1);
    s=size(L,2);
    down_result(i,2)=k/s;
    down_result(i,3)=1-binocdf(k-1,s,p0_down);
end 
down_result(:,1)=gene;
down_result(:,4)=mafdr(down_result(:,3),'BHFDR',true);   
down_DEG=down_result(down_result(:,4)<populationfdr,1);
overlap_gene=intersect(up_DEG,down_DEG)  
	overlap_gene_num=length(overlap_gene)
	if overlap_gene_num==0
	  up_DEG=up_DEG;
	  down_DEG=down_DEG
	end
	if overlap_gene_num>0
	  [index,~]=ismember(up_DEG,overlap_gene);
      up_DEG=up_DEG(~index);
	  [index,~]=ismember(down_DEG,overlap_gene);
      down_DEG=down_DEG(~index);    
    end
end



