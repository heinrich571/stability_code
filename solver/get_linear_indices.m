function linear_inds = get_linear_indices(mat, row_inds, column_inds)

% DESCRIPTION
%   This function returns the linear indices for the specified matrix based on specified row and column indices.
% INPUT
%   mat             matrix          [matrix]
%   row_inds        row indices     [vector]
%   column_inds     column indices  [vector]
% OUTPUT
%   linear_inds     linear indices corresponding to 'row_inds' and 'column_inds' in the input matrix

row_inds    = row_inds(:);
column_inds = column_inds(:);
linear_inds = sub2ind(size(mat), row_inds, column_inds);

end