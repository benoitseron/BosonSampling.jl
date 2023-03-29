### fft ###

# using FFTW

# probas_fourier

# probas_fourier_matrix = reshape(probas_fourier, ((n_max + 1), div(length(probas_fourier), (n_max + 1))))

# pdf_matrix = ifft(probas_fourier_matrix)

# pdf = reshape(pdf_matrix, (length(probas_fourier),))

# clean_pdf(pdf)

# bar(real.(pdf), alpha = 0.5)

# a = randn(20)

# ifft(fft(a))

