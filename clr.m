function y = clr(data)

y = log(data./repmat(geomean(data,2),1,size(data,2)));

end