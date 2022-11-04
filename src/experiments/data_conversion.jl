"""
    convert_csv_to_samples(path_to_file::String, m::Int, input_type = ThresholdModeOccupation, input_format = (input_data) -> ModeList(input_data, m))

Reads a CSV file containing experimental data. Converts it into a collection of samples of type `input_type` following the type of detectors available. `input_format` holds how the data is encoded, for instance as a `ModeList` or `ModeOccupation`.

Defaults work for files written as follow: each line is of form

    # 49, 1, 4, 5, 6, 7, 8, 10
    # for finding 49 occurences of the mode list (1, 4, 5, 6, 7, 8, 10)

For this line, the function will output 49 identical samples. This is not the most efficient use of memory but allows for an easier conversion with existing functions.

Example usage:

    convert_csv_to_samples("data/loop_examples/test.csv", 10)

"""
function convert_csv_to_samples(path_to_file::String, m::Int, input_type = ThresholdModeOccupation, input_format = (input_data) -> ModeList(input_data, m))

    data = readdlm(path_to_file, ',', Int)

    samples = Vector{input_type}()

    for (output_number, output_pattern) in enumerate(eachrow(data[:,2:end]))
        for count in 1:data[output_number,1]
            push!(samples, input_type(input_format(convert(Array{Int},output_pattern))))
        end
    end

    samples

end
