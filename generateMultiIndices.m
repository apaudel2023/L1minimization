function multi_indices = generateMultiIndices(num_vars, order)
    % Function to generate multi-indices for a specified order
    % Includes zero-order terms if order is 1, otherwise includes only terms of the specified order

    % Initialize the multi_indices matrix
    if order == 1
        multi_indices = zeros(1, num_vars); % Include zero-order for first order
    else
        multi_indices = []; % Do not include zero-order for higher orders
    end

    % Recursive function to fill the multi_indices matrix
    function iter(indices, pos, remaining_order)
        if pos == num_vars
            indices(end) = remaining_order;
            if sum(indices) == order
                multi_indices = [multi_indices; indices];
            end
        else
            for i = 0:remaining_order
                new_indices = indices;
                new_indices(pos) = i;
                iter(new_indices, pos + 1, remaining_order - i);
            end
        end
    end

    % Start the recursion
    iter(zeros(1, num_vars), 1, order);
end
