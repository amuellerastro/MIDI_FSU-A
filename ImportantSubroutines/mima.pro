function mima, arr

arr_new = arr[sort(arr)]

minval = min(arr_new)
maxval = max(arr_new)

return, [minval, maxval]


end