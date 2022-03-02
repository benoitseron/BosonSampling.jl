function get_parametric_type(i)

    """returns the types T1,... of a parametric type i::T{T1,...}

    if only one parameter, use
        get_parametric_type(i)[1]

    in case there is no parameter type, returns an array containing the type itself"""

    length(collect(typeof(i).parameters)) > 0 ? collect(typeof(i).parameters) : [typeof(i)]

end
