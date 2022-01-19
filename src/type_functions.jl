function get_parametric_type(i)

    """returns the types T1,... of a parametric type i::T{T1,...}

    if only one parameter, use
        get_parametric_type(i)[1]"""

    collect(typeof(i).parameters)

end
