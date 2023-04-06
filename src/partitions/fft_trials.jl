x_range = 0:1:10
f(x) = sech(x)

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