
%Script to determine sequence name from any mrd file
tmp = fread(fid)';
tmp = char(tmp);
tmp((tmp==44)|(tmp==10)|(tmp==13)|(tmp==0)) = 32;    % erase special characters
tmp = strtrim(regexp(tmp,' :','split'));  % clean up whitespace
tmp = regexp(tmp,'\s+','split');

 for j=2:size(tmp,2)-1,
                if isequal(tmp{1,j}{1},'VAR') || isequal(P{1,j}{1},'VAR_ARRAY'), P{1,j}(1) = []; end %builds a params structure
                eval(['params.' genvarname(P{1,j}{1}) ' = P{1,j}(2:end);']);
            end
            obj.param =  params;
frewind(fid);


