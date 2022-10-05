# convert CSV data into a vector of events

# conventions:
# the CSV file needs to have the same number of columns - make different files if using a different m

using DelimitedFiles
data = readdlm("data/test_file_experimental_data.csv", ',', Int)

for line in 1:size(data)[1]
    
