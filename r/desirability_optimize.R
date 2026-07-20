library(desirability)

range(design_coded$time)
range(design_coded$peak_rate)  
range(design_coded$conversion)

#desirability functions 
d_time <- dMin(low = 106.5533, high = 811.9060)
d_peak_rate<-dMin(low=0.8516844,high=4.0315215)
d_conv <- dMax(low = 0.98, high = 0.99)

predict(d_time, 811.91)  
predict(d_time, 106.55)  
predict(d_conv, 0.8267)  
predict(d_conv, 1.0)      
predict(d_peak_rate, 4.0315215) #we dont want it high because it needs more cooling


#the contatiner it keeps the three of them 
overall<-dOverall(d_time,d_conv,d_peak_rate)