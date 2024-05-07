function dispstatus(process_name, status)

separator = ' - ';
emphasizer_str = '';

if nargin < 2
    status = -1;
end

if nargin == 0
    status = -2;
    process_name = '';
end

switch status
    case -2
        status_str     = '';
        separator      = '';
        emphasizer_str = '';
    case -1
        status_str     = '';
        separator      = '';
        emphasizer_str = '----> ';
    case 0
        status_str = 'started';
    case 1
        status_str = 'completed';
    otherwise
        status_str = 'N/A';
end

message = [emphasizer_str status_str separator process_name '\n'];

fprintf(message)

end
