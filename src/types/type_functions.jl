"""
    get_parametric_type(i)

Return the types `T1`,..,`Tn` of a parametric type `i::T{T1,..,Tn}`.

!!! note
    - If the parametric type has only one parameter, use `get_parametric_type(i)[1]`.
    - If no parametric type, returns an array containing the type itself.
"""
function get_parametric_type(i)
    length(collect(typeof(i).parameters)) > 0 ? collect(typeof(i).parameters) : [typeof(i)]
end
