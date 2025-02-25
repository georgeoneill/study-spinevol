function mesh = go_readGmsh(file)

% GMESH reader, based on the specifications of the GMESH 2.0 file format as
% laid out in:
% https://bacula.nti.tul.cz/~jan.brezina/flow123d_doc/old_doc/msh_20.pdf

% Determine the number of nodes, which is written in the 5th line, and then
% the number of elements which is the 6+nNth line. (But remember things are
% zero indexed in the dlmread world!
nN          = dlmread(file,'',[5-1 1-1 5-1 1-1]);
nE          = dlmread(file,'',[7+nN 0 7+nN 0]);

% nodes are not necessarily stored in numerical order (!!!) so get the ids,
% and coordinate seperately. 
nodeID      = dlmread(file,'',[5 0 4+nN 0]); % These are 1-inedxed
nodes       = dlmread(file,'',[5 1 4+nN 3]);

% So in this instance we are expecting all the elements to be the same type,
% otherwise we may run into trouble later, lets check what type they are..
elementType      = dlmread(file,'',[8+nN 1 7+nN+nE 1]);
if numel(unique(elementType)) > 1
    error('There are multiple element classes in this file');
end

% For our use, we need there to be only one tag per element, 
% which represents the compartment ID, lets check that this is the case;
nT          = dlmread(file,'',[8+nN 2 7+nN+nE 2]);
if numel(unique(nT)) > 1
    error('Different elements have different number of tags');
elseif unique(nT) ~= 1
    error('Elements need to have only one tag');
end

% Okay, now lets work out how many nodes per element there are;
switch unique(elementType)
    case 4
        % Tetrahedron
        nPE = 4;
    case 5
        % Hexahedron
        nPE = 8;
    otherwise
        error('Unsupported element type, expecting either tetra- or hexahedrons');
end

% Element IDs are meant to be stored in numerical order, but lets be
% explicit incase this suddenly changes...
elementID   = dlmread(file,'',[8+nN 0 7+nN+nE 0]);
% Get element tags, and conver them from 0-index to 1-index
elementTag  = dlmread(file,'',[8+nN 3 7+nN+nE 3])+1;
% Get the nodes in each element
elementNodes = dlmread(file,'',[8+nN 4 7+nN+nE 3+nPE]); % 1-indexed

% Consolodate everything into a structure and return.
[~, loc] = sort(nodeID,'ascend');
mesh.pos    = nodes(loc,:); % puts nodes into numerical order. :)
[~, loc] = sort(elementID,'ascend');
switch unique(elementType)
    case 4
       mesh.tet = elementNodes(loc,:);
    case 5
       mesh.hex = elementNodes(loc,:);
end
mesh.tissue = elementTag(loc);

end

