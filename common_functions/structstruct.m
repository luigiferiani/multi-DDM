% structstruct(S) takes in a structure variable and displays its structure.
% 
% INPUTS:
% 
% Recursive function 'structstruct.m' accepts a single input of any class.
% For non-structure input, structstruct displays the class and size of the
% input and then exits.  For structure input, structstruct displays the
% fields and sub-fields of the input in an ASCII graphical printout in the
% command window.  The order of structure fields is preserved.
% 
% OUTPUTS:
% 
% (none yet!)

% Copyright (c) 2011, Andrew J. Joslin
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% Modified to improve handling of N-dimensional data by Luigi Feriani, 2018, luigi.feriani@gmail.com

function structstruct(S)

% Figure the type and class of the input
whosout = whos('S');
sizes = whosout.size;

% Luigi: modified the follwing 3 lines to display N-dimensional data (N>2)
sizestr = int2str(sizes(1));
for sc = 2:numel(sizes)
    sizestr = [sizestr,'x',int2str(sizes(sc))];
end

endstr = [':  [' sizestr '] ' whosout.class];

% Print out the properties of the input variable
disp(' ');
disp([inputname(1) endstr]);

% Check if S is a structure, then call the recursive function
% if isstruct(S)
    recursor(S,0,'');
% end

% Print out a blank line
disp(' ');

end



function recursor(S,level,recstr)

recstr = [recstr '  |'];

fnames = fieldnames(S);

for i = 1:length(fnames)
    
    %% Print out the current fieldname
    
    % Take out the i'th field
    tmpstruct = S.(fnames{i});
    
    % Figure the type and class of the current field
    whosout = whos('tmpstruct');
    sizes = whosout.size;
    
    %     sizestr = [int2str(sizes(1)),'x',int2str(sizes(2))];
    % Luigi: modified the following 3 lines to display N-dimensional data (N>2)
    sizestr = int2str(sizes(1));
    for sc = 2:numel(sizes)
        sizestr = [sizestr,'x',int2str(sizes(sc))];
    end

    endstr = [':  [' sizestr '] ' whosout.class];
    
    % Create the strings
    if i == length(fnames) % Last field in the current level
        str = [recstr(1:(end-1)) '''--' fnames{i} endstr];
        recstr(end) = ' ';
    else % Not the last field in the current level
        str = [recstr '--' fnames{i} endstr];
    end
    
    % Print the output string to the command line
    disp(str);
    
    %% Determine if each field is a struct
    
    % Check if the i'th field of S is a struct
    if isstruct(tmpstruct) % If tmpstruct is a struct, recursive function call
        recursor(tmpstruct,level+1,recstr); % Call self
    end
    
end

end