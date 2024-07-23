function uniqueArray = removeDuplicates(inputArray)
    % Round the input array to three decimal places
    roundedArray = round(inputArray, 3);
    
    % Sort the rounded array
    sortedArray = sort(roundedArray);
    
    % Initialize an empty array to store unique elements
    uniqueArray = [];
    
    % Iterate through the sorted array
    for i = 1:length(sortedArray)
        % If the current element is not equal to the previous element,
        % add it to the unique array
        if i == 1 || sortedArray(i) ~= sortedArray(i-1)
            uniqueArray = [uniqueArray, sortedArray(i)];
        end
    end
end
