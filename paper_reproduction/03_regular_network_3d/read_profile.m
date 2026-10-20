function d=read_profile(file)
raw=strrep(fileread(file),'D','E');
c=textscan(raw,'%f%f','Delimiter',',');d=[c{1} c{2}];d=d(all(isfinite(d),2),:);
[~,idx]=sort(d(:,1));d=d(idx,:);
end
