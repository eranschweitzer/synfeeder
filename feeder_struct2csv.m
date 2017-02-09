function feeder_struct2csv(n,e,fname)
extention = strfind(fname,'.csv');
if ~isempty(extention)
   fname = fname(1:extention-1); 
end
% convert node and edge feeder structures to csv
writetable(struct2table(n),[fname '_node.csv'])
writetable(struct2table(e),[fname '_edge.csv'])
