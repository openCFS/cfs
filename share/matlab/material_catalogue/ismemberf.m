function varargout = ismemberf(A, S, varargin)
 % function tf = ismemberf(A, S)
 % [tf loc] = ismemberf(A, S)
 %
 % Purpose: Floating-point ISMEMBER (i.e., with round-off tolerance)
 % See help ISMEMBER (Matlab buit-in function) for more detail about calling
 % Note that Matlab ISMEMBER uses strict exact comparison between floats.
 %
 % ismemberf(A, S, 'row') operates on rows of A and S (2D arrays)
 % and returns true (1) if they match, false (0) otherwise.
 % ismemberf(..., 'tol', tol) select a suitable the tolerance
 % - tol can be scalar or vector (using with 'row')
 % - When tol is vector, each element is applied to a specific 
 % column of A and S
 % - If not provided, or NaN, tol is 1e-10 relative to variations
 % of S values (separation for each column of S)
 % - When tol is provided as zero, ISMEMBERF calls Matlab ISMEMBER
 %
 % Examples:
 %
 % [tf, loc]=ismemberf(0.3, 0:0.1:1) % <- This return 0 by ISMEMBER
 %
 % [X Y]=meshgrid(0:0.1:10,0:0.1:10);
 % S=[3 1;
 % 3 3;
 % 5 6;
 % 5.2 5.5;
 % 3 3];
 % A = [X(:) Y(:)];
 % [tf loc]=ismemberf(A, S, 'row', 'tol', 0.5);
 % imagesc(reshape(loc,size(X)));
 %
 % Author: Bruno Luong <brunoluong@yahoo.com>
 % Original: 15/March/2009


 % Call the native Matlab ismember for strings
 if ~isnumeric(A) || ~isnumeric(S)
     out = cell(1,max(nargout,1));
     [out{:}] = ismember(A, S, varargin{:});
     varargout = out;
     return
 end

 % Preprocess the optional inputs (for parsing later on)
 vin = varargin;
 for k=1:length(vin)
     if ischar(vin{k})
         vin{k} = strtrim(lower(vin{k}));
     else
         vin{k}='';
     end
 end

 % parsing the varargin to set tol
 tol = NaN; % <- automatic tolerancce, see ismember1 bellow
 tolloc = strmatch('tol', vin);
 if ~isempty(tolloc)
     tolloc = tolloc(end);
     if length(vin)>=tolloc+1
         tol = varargin{tolloc+1};
         if ~isnumeric(tol)
             error('tolerance must be a number');
         end
     else
         error('tolerance value is not provided');
     end
 end

 % No tolerance, call Matlab ismember for array
 if all(tol==0)
     out = cell(1,max(nargout,1));
     % Remove tolerance parameters, not supported by Matlab
     v=varargin;
     v(tolloc:tolloc+1)=[];
     [out{:}] = ismember(A, S, v{:});
     varargout = out;
     return 
 end

 % one dimensional ismember
 if isempty(strmatch('row', vin))
     [tf loc] = ismember1(A, S, tol);
     
 else % row option
     % Check fo compatible dimension
     if ndims(A) ~= 2 || ndims(A) ~= 2 
         error('A and S must be 2-dimensional arrays');
     end
     if size(A,2) ~= size(S,2)
         error('A and S must have the same number of columns');
     end
     if isempty(S) % Few exception cases
         tf = false(size(A,1),1);
         loc = zeros(size(tf));
         if size(S,2)==0 % Hah, compare a 0-dimensional set is always true
             tf(:) = true;
             loc(:) = 1;
         end
     else % S must not be empty
         % duplicate tol if necessary
         if isempty(tol)
             tol = nan(size(A,2),1);
         elseif isscalar(tol)
             tol(1:size(A,2)) = tol;
         end
         
         % Loop over column (dimension)
         B = true;
         for j=1:size(A,2)
             % Get the binary matrix
             [tfj locj Bj] = ismember1(A(:,j), S(:,j), tol(j));
             B = B & Bj;
         end
         [iA iS] = find(B);
         tf = false(size(A,1),1);
         tf(iA) = true;
         if isempty(iS)
             loc = zeros(size(tf));
         else
             % This seems to be faster than max inside accumarray
             B = accumarray([iA(:),iS(:)], iS(:), size(B));
             loc = max(B,[],2);
             %loc = accumarray(iA(:),iS(:),size(tf),@(j) max(j));
         end
     end % if empty(S)
 end % row option

 % Process output
 [out{1:2}] = deal(tf, loc);
 nout = min(max(nargout,1), length(out));
 varargout(1:nout) = out(1:nout);

 end % ismemberf


 % Nested function, working linearly 
 function [tf loc B] = ismember1(A, S, tol)

     % Work only on subset of S, Su is sorted
     [Su I J] = unique(S(:),'last');
     
     % Set the tolerance automatically
     if isempty(tol) || isnan(tol)
         maxS = max(Su);
         minS = min(Su);
         if maxS == minS
             tol = 1e-10;
         else
             tol = (maxS - minS)*1e-10;
         end
     else
         tol=tol(1);
     end
     
     % Take a braket [-tol,tol] round each element
     xlo = Su - tol;
     xhi = Su + tol;
     % Look where points in A are located
     [dummy ilo] = histc(A, [xlo; Inf]);
     [dummy ihi] = histc(A, [-Inf; xhi]);
     % Test if they belong to the bracket
     tf = ilo & ihi & (ilo >= ihi);
     
     % ilo is the last indice
     loc = zeros(size(tf));
     loc(tf) = I(ilo(tf));
     %
     % Building a logical matrix of boolean B of size (m x n)
     % where m = numel(A), n = numel(S)
     % B(m,n) is true if two elements in A(m) and S(n) is "identical"
     %
     if nargout>=3
         % index in S
         % Group all the index of S when they map to the same Su
         left = ihi(tf);
         right = ilo(tf);

         % Find the index in S, duplicated by number of elements in A
         % belong to it
         [iS nele] = getiS(left(:), right(:), J);
     
         % index in A
         % This is a trick to generate the same vector long as iS
         % with a ids (1, 2, 3...) repeated for each cell elements (of
         % length nele)
         iA = cumsum(accumarray(cumsum([1; nele]),1));
         
         inonly = find(tf);
         iA = inonly(iA(1:end-1));
         
         % Logical matrix, in sparse
         B = sparse(iA(:),iS(:),true,numel(A),numel(S));
     end % if nargout>=3
     
 end % ismember1

 function [Icat lengthII] = Jsubset(J)
 % J is an array from the third argument of UNIQUE (mapping index)
 % Group the mapping J by subset (in Icat), each subset has the length
 % stored in lengthII. The subset are sorted.
 % In other word, perform equivalent to the following:
 % (but optimized for speed)
 % II = accumarray(J(:), (1:numel(J)).', [max(J) 1], @(x) {x});
 % Icat = cat(1,II{:});
 % lengthII = cellfun(@length, II)

 [Js Icat]=sort(J(:));
 n=Js(end);
 m=length(J);

 jump=diff([0; Js(:)])>0;
 last=zeros(n,1);
 last(Js(jump)) = diff([find(jump); m+1]);
 lengthII=diff([0; cumsum(last)]);

 end % Jsubset

 function [v csl] = catbraket(l, r)
 % Concatenate I1:=(l(1):r(1))', I2=(l(2):r(2)', etc ... 
 % in v = [I1; I2; ... ]
 % Note: at the entrance r(i) must be l(i)-1 for empty braket

     if isempty(l)
         v = [];
         csl = 1;
         return
     end
     l=l(:);
     r=r(:);
     csl=cumsum([0; r-l]+1);

     v = accumarray(csl(1:end-1), (l-[0; r(1:end-1)]-1));
     % Adjust the size of v
     sv = csl(end)-1; % wanted length 
     if size(v,1)>sv
         v(sv+1:end)=[];
     elseif size(v,1)<sv % pad zero
         v(sv)=0;
     end
     v = cumsum(v+1);
     %[l r v]
 end % catbraket

 function [iS nele] = getiS(left, right, J)
 % Do this (but avoid cell to improve speed)
 % iS = arrayfun(@(l,r) cat(1,II{l:r}).', left, right, ...
 % 'UniformOutput', false);
 % nele = cellfun(@length, iS);
 % iS = [iS{:}]; % concatenate in a long row vector
 % This is awfully hard to read, because of the optimization
     
     [Icat lengthII] = Jsubset(J);
     % Do the following
     % is1 = arrayfun(@(l,r) (l:r).', left, right, ...
     % 'UniformOutput', false);
     % is1 = cat(1,is1{:});
     [is1 csl] = catbraket(left, right); 
     
     % Compute the length of each subset in iS
     % nele(k) will be length of cat(1,II{left(k):right(k)})
     csIIis = cumsum([0; lengthII(is1)]);
     nele = csIIis(csl(2:end))-csIIis(csl(1:end-1));
     
     % Build the left/right brakets when II cells will be expanded
     ss=cumsum([0; lengthII])+1;
     l=ss(is1);
     r=ss(is1+1)-1; 
     
     % Last step, concatenate II and retrieve data in the braket
     % 
     iS = Icat(catbraket(l, r));

 end % getiS