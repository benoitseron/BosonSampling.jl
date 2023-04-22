x_range = 0:1:10
f(x) = (x)^2 - 2x + 10

f_array = f.(x_range)

plot(x_range, f_array)

ft_array = fft(f_array)

plot(real.(fft(f_array)))
plot!(imag.(fft(f_array)))

# pad with zeros ft_array for fourier transform

ft_array_padded = similar(ft_array, 2 .* length(ft_array))
ft_array_padded[1:length(ft_array)] = ft_array
ft_array_padded[length(ft_array)+1:end] = zeros(length(ft_array))

ft_array_padded

plot!(real.(ft_array_padded))
plot!(imag.(ft_array_padded))

f_array_recovered = real.(ifft(ft_array))
f_array_recovered_padded = real.(ifft(ft_array_padded))



plot(f_array)
plot!(f_array_recovered)
plot!(f_array_recovered_padded)

f_array_recovered


### padding with zeros the signal ###

f_array_padded = zeros(eltype(f_array), 2 .* length(f_array))
f_array_padded[1:length(f_array)] = f_array

plt_real = plot(f_array_padded, label="padded")

ft_array_padded = fft(f_array_padded)

plt_fourier = plot(real.(fft(ft_array_padded)))
plot!(plt_fourier, imag.(fft(ft_array_padded)))

f_array_recovered_padded = real.(ifft(ft_array_padded))

plot!(plt_real, f_array_recovered_padded, label="recovered")

### padding with zeros only the fourier transform ###

# we obtain the fourier coefficients in our theoretical calculations

