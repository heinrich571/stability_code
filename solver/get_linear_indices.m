function Linear_Indices = get_linear_indices(Matrix, Row_Indices, Column_Indices)

% DESCRIPTION
%   This function returns the linear indices for the specified matrix
%   based on specified row and column indices.
% 
% INPUT
%   Matrix          matrix          [matrix]
%   Row_Indices     row indices     [vector]
%   Column_Indices  column indices  [vector]
% 
% OUTPUT
%   Linear_Indices  linear indices corresponding to "Row_Indices" and
%   "Column_Indices" in the "Matrix".

% Input validation
if ~isequal(size(Row_Indices), size(Column_Indices))
    error('get_linear_indices:SizeMismatch', ...
          'Inputs for "Row_Indices" and "Column_Indices" must be the same size.')
end

% Calculate linear indices of matrix entries
row_inds    = Row_Indices(:);
column_inds = Column_Indices(:);

Linear_Indices = sub2ind(size(Matrix), row_inds, column_inds);

end
