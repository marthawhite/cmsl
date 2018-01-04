function string = struct2str(s)
% Takes the fields in the struct and prints them
% in one line.
  
  fn=fieldnames(s);
  sc=struct2cell(s(1));
  
  string = [];
  if length(fn) > 0
    string = [fn{1} '=' getString(sc{1})];
  end 
  for i = 2:length(fn)
    string = [string ', ' fn{i} '=' getString(sc{i})];        
  end
  
  function s = getString(field_value)
    s = '';
    if isnumeric(field_value) || islogical(field_value)
      s = num2str(field_value);    
    elseif isstr(field_value)
      s = field_value;
    elseif isa(field_value, 'function_handle')
      s = func2str(field_value);
    else
      s = 'unknown';
    end    
  end
  
end

