function statlimit = FindStatLimit(tw_bc,s_boot,timeaxis)
my_array=s_boot(find(timeaxis==tw_bc(1)):find(timeaxis==tw_bc(2)));
k = 1;
thiscluster=1;
while k<numel(my_array)
    if my_array(k) == NaN
    k = k+1;
    elseif (my_array(k+1) == (my_array(k)))
        g = k;
        consec = 1;
        while g<numel(my_array) && my_array(g+1)==my_array(g) 
        consec = consec+1;
        g=g+1;
        end
        num(thiscluster) = consec;
        thiscluster=thiscluster+1;
        k=g+1;
    else
    k=k+1;
    end
end
if ~exist('num','var')
    statlimit = 0;
else
statlimit = max(num);
end
end