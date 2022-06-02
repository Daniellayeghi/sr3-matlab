function arr = remove_mid_elems(arr, times)
    for sparse = 1:times
        arr(1:2:end-1, :) = []; 
    end
end