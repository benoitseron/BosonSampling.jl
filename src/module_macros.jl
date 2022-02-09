#
# Copyright 2021 Clemens Cords
# Created on 26.12.2021 by clem (mail@clemens-cords.com)
#

"""
export all non-temporary, non-imported values (values not having a '#' at the
start of it's symbol when listed via Base.names) in a module

Example:
    module MyModule

        module SubModule
            _sub_variable
        end

        _variable
    end

    @make_public MyModule

    # exports:
    #   MyModule
    #   MyModule.SubModule
    #   MyModule._variable

See also: @make_public_rec
"""
macro make_public(module_name::Symbol)

    eval(Meta.parse("export " * string(module_name)))

    as_module = eval(module_name)
    @assert as_module isa Module

    for name in names(as_module; all = true)
        if (string(name)[1] != '#')
            #println("export " * string(name))
            as_module.eval(Meta.parse("export " * string(name)))
        end
    end

    return nothing
end
export make_public

"""
export all non-temporary, non-imported values (values not having a '#' at the
start of it's symbol when listed via Base.names) in a module and all
such values in any submodule, recursively

Example:
    module MyModule

        module SubModule
            _sub_variable
        end

        _variable
    end

    @make_public_rec MyModule

    # exports:
    #   MyModule
    #   MyModule.SubModule
    #   MyModule.SubModule._sub_variable
    #   MyModule._variable

See also: @make_public
"""
macro make_public_rec(module_name::Symbol)

    function make_public_aux(child_name::Symbol, super::Module) ::Nothing

        if (string(child_name)[1] != '#')
            #println("export " * string(child_name))
            super.eval(Meta.parse("export " * string(child_name)))
        end

        child = super.eval(child_name)
        if (child isa Module && child != super)
            for name in names(child; all = true)
                make_public_aux(name, child)
            end
        end

        return nothing
    end

    origin = eval(module_name)
    @assert origin isa Module

    for name in names(origin; all = true)
        make_public_aux(name, origin)
    end

    return nothing
end
export make_public_rec
