function [] = general_print( f_id, varargin )
assert(isempty(varargin) == false, 'Not enought inputs for function')

if (length(varargin) == 1)
    fprintf(varargin{1})
    fprintf(f_id, varargin{1});
else
    fprintf(varargin{1}, varargin{2:end})
    fprintf(f_id, varargin{1}, varargin{2:end});
end


end

