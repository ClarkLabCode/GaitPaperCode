function s = scattern(x,sz,c,varargin)
%SCATTERN: N-dimensional scatter plots - JZV

switch size(x,2)
    case 1
        s = scatter( x,rand(length(x),1)-0.5, sz, c, varargin{:});
    case 2
        s = scatter(x(:,1), x(:,2), sz, c, varargin{:});
    case 3
        s = scatter3(x(:,1), x(:,2), x(:,3), sz, c, varargin{:});
    otherwise
        error('Invalid data dimensionality: %d.', size(x,2));
end

end

