function A = qs_mergestruct(A, B)
%QS_MERGESTRUCT  Recursive struct merge: fields of B override fields of A.
%
%   Used to build the resolved configuration as
%       defaults  <-  case.json  <-  study.json  <-  command-line overrides
%   Scalar structs recurse; everything else (numbers, char, cell, arrays,
%   struct arrays) is replaced wholesale by B's value.
%
%   Replacing rather than merging non-structs is deliberate: if a study says
%   "pc_kPa": {...}, it means THAT grid, not that grid unioned with whatever
%   the defaults happened to contain.

    if isempty(B), return; end
    assert(isstruct(A) && isstruct(B), 'qs_mergestruct: both inputs must be structs');
    assert(isscalar(A) && isscalar(B), 'qs_mergestruct: struct arrays are not merged');

    f = fieldnames(B);
    for k = 1:numel(f)
        key = f{k};
        vB  = B.(key);
        if isfield(A, key) && isstruct(A.(key)) && isstruct(vB) ...
                           && isscalar(A.(key)) && isscalar(vB)
            A.(key) = qs_mergestruct(A.(key), vB);
        else
            A.(key) = vB;
        end
    end
end
