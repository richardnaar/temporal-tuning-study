function [allns] = colFun(n, doSorting)
allns = [];
while n ~= 1
    if  mod(n, 2) == 1 %
        n = (n*3+1);
        allns = [allns, n];
    else
        n = n/2;
        allns = [allns, n];
    end
end
if doSorting == 1
    allns = sort(allns, 'ascend');
end
end