clear all
clc
data = create_data_quadratic;
[~, ~, ~] = BEM_collocation_quadratic(data);