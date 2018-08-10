#function to calculate correlation coefficient
corrcoef=function(x,y){
  corr_coef=sum(x*y)/(sqrt(sum(x*x))*sqrt(sum(y*y)))
  return(corr_coef)
}

##function for imf selection and recombination based on correlation coefficient
imf_selection_sum=function(imf_data,data_raw1){
  cc <- apply(imf_data, 2, function(x) corrcoef(data_raw1, x))
  #threshold selection method
  #A CRITERION FOR SELECTING RELEVANT INTRINSIC MODE FUNCTIONS IN EMPIRICAL MODE DECOMPOSITION
  cc_threshold=mean(abs(cc))*0.95
  imf_data1 <- imf_data[, cc > cc_threshold]
  #princ_comp <- princomp(imf_data1)
  imf_select=rowSums(imf_data1, na.rm = FALSE, dims = 1)#reconstruction by recombing the selected IMFs
  #imf_select=as.data.frame(princ_comp$scores[,1])
  return(imf_select)
}

plot.frequency.spectrum <- function(X.k, fs, xlimits = c(0, 1000)) {
  plot.data  <- data.frame(cbind((0:(length(X.k)-1))*fs/(2*length(X.k)), Mod(X.k)))
  names(plot.data) <- c("freq", "strength")
  plot.data[1,] <- 0
  # TODO: why this scaling is necessary?
  #plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  plot.data$top_str <- NA
  top_str_pos <- order(plot.data$strength, decreasing = T)[1:20]
  plot.data$top_str[top_str_pos] <- round(plot.data$freq[top_str_pos])
  plot.data <- plot.data[1:1000, ]
  plt <- ggplot(plot.data, aes(x = freq, y = strength)) + geom_line() + 
    geom_text(aes(label = top_str)) + 
    scale_x_discrete(breaks = seq(0, nrow(plot.data), 50)) 
  return (plt)
  #   plot(plot.data, t="h", lwd=2, main="", 
  #        xlab="Frequency (Hz)", ylab="Strength", 
  #        xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}

library(moments)
imf_kurtosis_selection <- function(imf_data){
  kurt <- apply(imf_data, 2, kurtosis)
  zero_crosses <- apply(imf_data, 2, function(x) sum(diff(sign(x))!=0))
  high_kurt_pos <- which(kurt > 2.9 & kurt < 15 & zero_crosses > 6)
  imf_select <- data.frame(rowSums(imf_data[,high_kurt_pos], dims = 1))
  return(imf_select)
}

imf_local_kurtosis_selection <- function(imf_data, fs, rps){
  wl <- floor(fs/rps)
  local_kurt <- matrix(NA, nrow = floor(nrow(imf_data)/wl)-1, ncol = ncol(imf_data))
  for(i in 1:ncol(imf_data)){  
    for(j in 1:floor(nrow(imf_data)/wl)-1){
      local_kurt[j,i] <- kurtosis(imf_data[((j-1)*wl+1):(j*wl),i])
    }
  }
  mean_kurt <- apply(local_kurt, 2, mean)
  imf_select <- data.frame(rowSums(imf_data[,mean_kurt >= 2.9], dims = 1))
}

axis_select <- function(raw_data){
  princ_comp <- princomp(raw_data[,1:3])
  prim_comp <- princ_comp$scores[,1]
  corr_coef <- apply(raw_data[,1:3], 2, function(x) corrcoef(x, prim_comp))
  return(raw_data[,which(corr_coef == max(corr_coef))])
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}