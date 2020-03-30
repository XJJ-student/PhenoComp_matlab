function [def_all,def_up,def_down]=fdr_adjust_bh(gene,p,index,fdr_cutoff)
m=length(p);
[p_sort,loc]=sort(p);
indexx(:,1)=(1:m).*(p_sort'<=(1:m)*fdr_cutoff/m);
n_sig=max(indexx);
if (n_sig)
    sig_loc=loc(1:n_sig);
    def_all=gene(sig_loc,:);
    index=index(sig_loc);
    def_up=def_all(index==1);
    def_down=def_all(index==0);
else
    def_all=[];
    def_up=[];
    def_down=[];
end
end