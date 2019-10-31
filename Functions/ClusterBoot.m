% function [s_boot_corrected;clusters] = ClusterBoot(s_boot,statlimit,timeaxis)
my_array = s_boot;
my_num = statlimit;
clusters = [];numclu = 0;
consec = 1;

my_array(1:(find(timeaxis==tw_bc_PDR(2))-1))=NaN;
i = find(timeaxis==tw_bc_PDR(2));
while i<numel(my_array)
    if my_array(i) == NaN
        i = i+1;
    elseif ( my_array(i+1) == (my_array(i)) )
        numclu = numclu + 1;
        g = i;
        consec = 1;
        clusters(numclu,1) = timeaxis(g);
        while g<numel(my_array) && my_array(g+1)==my_array(g)
            consec = consec+1;
            g=g+1;
        end
        if consec<=statlimit
            my_array(1,i:g)=NaN;
            clusters(numclu,2) = NaN;
            i=g+1;
        else
            i=g+1;
            clusters(numclu,2) = timeaxis(g);
        end
    else
        i=i+1;
    end
end
s_boot_corrected=my_array;
