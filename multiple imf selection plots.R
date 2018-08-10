library(pracma)
library(seewave)
library(signal)
library(Rlibeemd)
#sample_rate <- 20480
sample_rate <- 1365 #CHANGED  FROM 1344
nq_freq <- sample_rate/2
#rpm = 2000
rpm = 987
rps = rpm/60
Path1="E:/nanoprecise/test2/2nd_test/2nd_test"
Files = list.files(Path1,full.names = 1)
#data_raw <- read.table(Files[530], sep="\t")
#selected_imf <- data_raw[,1]
data_raw <- read.csv("E:/nanoprecise/s3-files_download/pump_de_comp539-2018-08-07T16-17-10.csv")
freq=seq(0,(sample_rate-sample_rate/nrow(data_raw)),sample_rate/nrow(data_raw))
f=freq[1:(length(freq)/2)]
selected_axis <- data_raw[,1]
###########removing mean offset from y-axis ##################
selected_axis <- selected_axis - mean(selected_axis)

#amms <- 9.80665*1000*data_raw[,1]
#a_time <- (0:(length(amms)-1))/sample_rate
#vmms <- cumtrapz(a_time, amms)
#bf <- butter(5, 0.2)
#v_filtered <- filter(bf, vmms)
#imf_data=eemd(data_raw[,1], num_imfs = 14, ensemble_size = 250L, noise_strength = 0.2,
#                 S_number = 4L, num_siftings = 50L, rng_seed = 0L, threads = 0L)
#imf_data=ceemdan(data_raw[,1], num_imfs = 14, ensemble_size = 250L, noise_strength = 0.2,
#             S_number = 4L, num_siftings = 50L, rng_seed = 0L, threads = 0L)
imf_data=ceemdan(selected_axis, num_imfs = 14, ensemble_size = 250L, noise_strength = 0.2,
                              S_number = 4L, num_siftings = 50L, rng_seed = 0L, threads = 0L)

########### imf selection through Correlation coefficients#################
#selected_imf <- imf_kurtosis_selection(imf_data)
#no_dc_imf <- rowSums(imf_data[,1:10], dims = 1)
selected_imf <- imf_selection_sum(imf_data, selected_axis)
#selected_imf <- imf_local_kurtosis_selection(imf_data, sample_rate, rps)
#princ_comp <- princomp(imf_data[,-14])
#selected_imf <- princ_comp$scores[,1]
#selected_imf <- imf_data[,5]
selected_imf_hilb=abs(hilbert(selected_imf, sample_rate, fftw= FALSE))

imf_data_fft=(2/nrow(data_raw))*abs(fft(selected_imf[,1]))
imf_data_fft=imf_data_fft[1:length(f)]
#plot.frequency.spectrum(imf_data_fft)

imf_data_fft_hilb=(2/nrow(data_raw))*abs(fft(selected_imf_hilb[,1]))
imf_data_fft_hilb=imf_data_fft_hilb[1:length(f)]
p1 <- plot.frequency.spectrum(imf_data_fft_hilb, sample_rate)
p1 <- p1 + ggtitle("Correlation-coeff")

############ imf selection through global kurtosis ################
selected_imf <- imf_kurtosis_selection(imf_data)
#selected_imf <- imf_selection_sum(imf_data, data_raw[,1])
#selected_imf <- imf_local_kurtosis_selection(imf_data, sample_rate, rps)
#princ_comp <- princomp(imf_data[,-14])
#selected_imf <- princ_comp$scores[,1]
#selected_imf <- imf_data[,5]
selected_imf_hilb=abs(hilbert(selected_imf, sample_rate, fftw= FALSE))

imf_data_fft=(2/nrow(data_raw))*abs(fft(selected_imf[,1]))
imf_data_fft=imf_data_fft[1:length(f)]
#plot.frequency.spectrum(imf_data_fft)

imf_data_fft_hilb=(2/nrow(data_raw))*abs(fft(selected_imf_hilb[,1]))
imf_data_fft_hilb=imf_data_fft_hilb[1:length(f)]
p2 <- plot.frequency.spectrum(imf_data_fft_hilb, sample_rate)
p2 <- p2 + ggtitle("global kurtosis")

############ imf selection through local kurtosis ################
#selected_imf <- imf_kurtosis_selection(imf_data)
#selected_imf <- imf_selection_sum(imf_data, data_raw[,1])
selected_imf <- imf_local_kurtosis_selection(imf_data, sample_rate, rps)
#princ_comp <- princomp(imf_data[,-14])
#selected_imf <- princ_comp$scores[,1]
#selected_imf <- imf_data[,5]
selected_imf_hilb=abs(hilbert(selected_imf, sample_rate, fftw= FALSE))

imf_data_fft=(2/nrow(data_raw))*abs(fft(selected_imf[,1]))
imf_data_fft=imf_data_fft[1:length(f)]
#plot.frequency.spectrum(imf_data_fft)

imf_data_fft_hilb=(2/nrow(data_raw))*abs(fft(selected_imf_hilb[,1]))
imf_data_fft_hilb=imf_data_fft_hilb[1:length(f)]
p3 <- plot.frequency.spectrum(imf_data_fft_hilb, sample_rate)
p3 <- p3 + ggtitle("local kurtosis")

############ imf selection through princ components ################
#selected_imf <- imf_kurtosis_selection(imf_data)
#selected_imf <- imf_selection_sum(imf_data, data_raw[,1])
#selected_imf <- imf_local_kurtosis_selection(imf_data, sample_rate, rps)
#princ_comp <- princomp(imf_data[,-14])
#selected_imf <- princ_comp$scores[,1]
selected_imf <- imf_data[,1]

selected_imf_hilb=abs(hilbert(selected_imf, sample_rate, fftw= FALSE))

imf_data_fft=(2/nrow(data_raw))*abs(fft(selected_imf[,1]))
imf_data_fft=imf_data_fft[1:length(f)]
#plot.frequency.spectrum(imf_data_fft)

imf_data_fft_hilb=(2/nrow(data_raw))*abs(fft(selected_imf_hilb[,1]))
imf_data_fft_hilb=imf_data_fft_hilb[1:length(f)]
p4 <- plot.frequency.spectrum(imf_data_fft_hilb, sample_rate)
p4 <- p4 + ggtitle("first IMF")

multiplot(p1,p2,p3,p4,cols = 2)
