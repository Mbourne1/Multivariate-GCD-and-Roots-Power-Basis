function [] = verify_yes_no_input(input)

if strcmp(input,'y') || strcmp(input,'n')
else
    error('Input must be either y or n')
end

end
