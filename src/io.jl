using CSV, DataFrames

function make_directory(path; delete_contents_if_existing = true)

    """makes a directory and goes inside, removes the previous contents of a
    similarly named directory if existing"""

    try
        if delete_contents_if_existing
            rm(path, recursive = true)
        else
            pass()
        end
    catch err

    finally
        mkpath(path)
        cd(path)
    end
end

function print_error(err)

	"""prints the error in a try/catch"""
	buf = IOBuffer()
    showerror(buf, err)
    message = String(take!(buf))
end

# to save julia variables : using JLD2

function save_matrix(mat, filename = "matrix.csv")

	CSV.write(filename,  DataFrame(mat), header=false)

end
