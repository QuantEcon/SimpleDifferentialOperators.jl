# find i such that xs[i] is nearest to x
# ASSUMES xs is sorted (from least to largest) and xs[1] <= x
function findnearestindex(xs, x)
    n = length(xs)
    # using binary search algorithm O(log n), search the last element found
    # i.e., xs[i_first] <= x < xs[min(i_first + 1, n)] and i_first >= 1
    i_first = max(1, searchsortedlast(xs, x))
    i_next = min(i_first + 1, n)
    # suffices to compare xs[i_first] and xs[i_next]
    return ((x - xs[i_first]) < (xs[i_next] - x)) ? i_first : i_next
end