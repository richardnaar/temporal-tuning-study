function data = wcond(data, avgDim) % standardize within conditions over all bins and subjects
    datMean = mean(data,avgDim);
    numer = bsxfun(@minus, data, min(datMean));
    data = bsxfun(@rdivide, numer, max(datMean) - min(datMean));
end