function direc=mydir(inputdir,searchstrings,addtasksflg)
% this function provides an easy way to get the data inside a directory
%inputdir= string input of the searched folder
%searchstrings=cell string of all data
%addtasksflg: 1=

% case 1) if we have no input data we just want an overview over the
% current directory

if nargin==0 
    inputdir=cd;    
end

% (implicit) case 2) we have just an input image and want the list - this is
% the "standard" output
direc=struct2cell(dir(inputdir));
direc=direc([1:2,5],3:end)';

% case 4) we have specific output requirements
if nargin>2
    switch addtasksflg
        case 1
            for a=1:size(direc,1)                
                direc{a,1}=fullfile(direc{a,2},direc{a,1});
            end
            direc=direc(:,1);
        case 2
            direc=direc(:,1);
    end
end

% case 3): we have a search string
if nargin>1 && ~isempty(searchstrings)
    % we have a search string
    if ischar(searchstrings) %if we have a character string we just run that
        direc=direc(contains(direc(:,1),searchstrings),:);
    elseif iscell(searchstrings)
        for a=1:length(searchstrings) %if we have more strings, we loop them
            direc=direc(contains(direc(:,1),searchstrings{a}),:);
        end
    elseif islogical(searchstrings)
        direc=direc(searchstrings,:);
    end
end

end