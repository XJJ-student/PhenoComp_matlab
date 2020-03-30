function [stable_pairs]=stablepair(gid,exp,freq)
len=length(gid)-1;  
m=length(exp(1,:)); 
pair=nan(nchoosek(length(gid),2),3); 
j=1;
for i = 1 : len
    tic
    time1=len-i 
    gid_pair=[gid(i).*ones(len-i+1,1),gid(i+1:end,:)]; 
    RE=exp(i*ones(len-i+1,1),:)-exp(i+1:end,:)>0; 
    index=sum(RE,2)/size(RE,2)>0.5;
    f=sum(RE,2)/size(RE,2);
    pair1=[gid_pair(index,:),f(index,:)];
    pair2=[gid_pair(~index,[2,1]),1-f(~index,:)];
    pair(j:j+length(RE(:,1))-1,:)=[pair1;pair2];
    j=j+length(RE(:,1));
    toc
end
stable_pairs=pair(pair(:,3)>freq,1:2);  
end

