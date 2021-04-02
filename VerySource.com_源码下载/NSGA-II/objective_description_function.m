function [number_of_objectives, number_of_decision_variables, min_range_of_decesion_variable, max_range_of_decesion_variable] = objective_description_function()

%% function [number_of_objectives, number_of_decision_variables, min_range_of_decesion_variable, max_range_of_decesion_variable] = objective_description_function()
% This function is used to completely describe the objective functions and
% the range for the decision variable space etc. The user is prompted for
% inputing the number of objectives, numebr of decision variables, the
% maximum and minimum range for each decision variable and finally the
% function waits for the user to modify the evaluate_objective function to
% suit their need.

g = sprintf('Input the number of objective: ');
% Obtain the number of objective function
number_of_objectives = input(g);
if number_of_objectives < 2
    error('This is a multi-objective optimization function hence the minimum number of objectives is two');
end
g = sprintf('\nInput the number of decision variables: ');
% Obtain the number of decision variables
number_of_decision_variables = input(g);
clc
for i = 1 : number_of_decision_variables
    clc
    g = sprintf('\nInput the minimum value for decision variable %d : ', i);
    % Obtain the minimum possible value for each decision variable
    min_range_of_decesion_variable(i) = input(g);
    g = sprintf('\nInput the maximum value for decision variable %d : ', i);
    % Obtain the maximum possible value for each decision variable
    max_range_of_decesion_variable(i) = input(g);
    clc
end
g = sprintf('\n Now edit the function named "evaluate_objective" appropriately to match your needs.\n Make sure that the number of objective functions and decision variables match your numerical input. \n Make each objective function as a corresponding array element. \n After editing do not forget to save. \n Press "c" and enter to continue... ');
% Prompt the user to edit the evaluate_objective function and wait until
% 'c' is pressed.
x = input(g, 's');
if isempty(x)
    x = 'x';
end
while x ~= 'c'
    clc
    x = input(g, 's');
    if isempty(x)
        x = 'x';
    end
end    