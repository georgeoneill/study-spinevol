function h = fem_plot_mesh(node,elem,axis,pm,val)

% George O'Neill, 2024


% Inputs
% node - vertices within the mesh
% elem - which nodes are connected and tissue ID
% axis - which axis do we want to slice along (takes a value between 1-3)
% pm - which direction are we going to omit (takes a '>' or '<' string)
% val - which point along axis do we make the cut?

if nargin == 2
    axis = 1;
    pm = '<';
    val = 0;
end

types = [];

switch size(elem,2)
    case 4
        % just tetrahedra
        els = elem;
    case 5
        % tetrahedra with specific tissue types
        types = elem(:,5);
        els = elem(:,1:4);
    otherwise
        error('unsupported element type!')
end


% for this example, want to look at everything in the positive X direction

switch pm
    case '>'
bad_node = find(node(:,axis) < val);
    case '<'
  bad_node = find(node(:,axis) > val);
end

bad_el = find(sum(ismember(els,bad_node),2)>0);

good_el = els;
good_el(bad_el,:) = [];
good_types = types;
good_types(bad_el) = [];

F=[good_el(:,[2 1 3]);... %face 1 2 3
    good_el(:,[1 2 4]);... %face 1 2 4
    good_el(:,[2 3 4]);... %face 2 3 4
    good_el(:,[3 1 4])];   %face 1 3 4
C=repmat(good_types,4,1);
% [faces, c] = element2patch(good_el,good_types);
cmap = hsv(numel(unique(good_types)));
% cmap = flipud(brewermap(3,'BuPu'));

figure
patch('vertices',node,'faces',F,'facevertexcdata',cmap(C,:),'facecolor','flat');

end

% function [varargout]=element2patch(varargin)

