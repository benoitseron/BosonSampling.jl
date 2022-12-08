"""
    convert_csv_to_samples(path_to_file::String, m::Int, input_type = ThresholdModeOccupation, input_format = (input_data) -> ModeList(input_data, m), samples_type = MultipleCounts, renormalize = false)

Reads a CSV file containing experimental data. Converts it into a collection of samples of type `input_type` following the type of detectors available. `input_format` holds how the data is encoded, for instance as a `ModeList` or `ModeOccupation`.

You can choose the encoding of the samples to be either a Vector of observed events (like detector readings) or a `MultipleCounts` if there are many events.

If using `ModeOccupation`, defaults work for files written as follow: each line is of form

    # 45, 1,1,0
    # for finding 45 occurences of the mode occupation 1,1,0

and you should use

    input_format = (input_data) -> ModeOccupation(input_data)

If using `ModeList`, defaults work for files written as follow: each line is of form

    # 49, 1, 4, 5, 6, 7, 8, 10
    # for finding 49 occurences of the mode list (1, 4, 5, 6, 7, 8, 10)

and you should use

    input_format = (input_data) -> ModeList(input_data, m)

For this line, the function will output 49 identical samples, or a `MultipleCounts`. This is not the most efficient use of memory but allows for an easier conversion with existing functions.

Example usage:

    convert_csv_to_samples("data/loop_examples/test.csv", 10)

Setting `renormalize = false` outputs the number of times an event was observed, while true gives a probability of this particular event.

"""
function convert_csv_to_samples(path_to_file::String, m::Int, input_type = ThresholdModeOccupation, input_format = (input_data) -> ModeOccupation(input_data); samples_type = MultipleCounts, renormalize = false)

    data = readdlm(path_to_file, ',', Int)

    samples = Vector{input_type}()

    if samples_type == Vector

        for (output_number, output_pattern) in enumerate(eachrow(data[:,2:end]))
            for count in 1:data[output_number,1]
                push!(samples, input_type(input_format(convert(Array{Int},output_pattern))))
            end
        end

        return samples

    elseif samples_type == MultipleCounts

        counts = renormalize ? Vector{Real}() : Vector{Int}()

        for (output_number, output_pattern) in enumerate(eachrow(data[:,2:end]))

                push!(samples, input_type(input_format(convert(Array{Int},output_pattern))))
                push!(counts, data[output_number,1])

        end

        return MultipleCounts(samples, renormalize ? counts ./ sum(counts) : counts) 
    else
        error("not implemented")
    end

    return nothing

end
