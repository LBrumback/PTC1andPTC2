## 8/22/19
## Example R script for estimating PTC1 and PTC2
## Brumback LC, Jacobs Jr DR, Duprez DA. PTC1 and PTC2: New indices of blood pressure waveforms and
## cardiovascular disease. Am. J. Epidemiology (under review) 
## LynB@uw.edu

## Example beat-specific pressure waveform
pressure<-c(83.5,83.6,84.6,86.1,87.7,89.8,92.3,95.5,99.1,103.0,107.4,112.1,117.0,122.2,127.1,131.8,136.4,
            140.1,143.4,146.5,148.9,151.2,152.5,153.8,154.8,155.7,156.2,156.9,156.9,157.2,157.4,157.7,157.7,
            157.5,157.7,157.5,157.5,157.4,157.7,156.9,156.7,156.4,155.9,155.1,154.8,153.8,153.1,151.8,151.2,
            150.0,149.2,148.1,146.6,145.5,143.7,142.2,140.6,138.7,136.5,134.9,133.0,130.7,129.2,127.6,126.1,
            124.7,123.3,121.9,121.4,120.4,119.9,119.1,118.5,118.1,118.1,117.7,117.7,117.3,117.7,117.5,117.5, 
            117.5,117.2,117.0,117.2,116.8,116.8,116.8,116.3,116.0,116.2,115.4,115.0,114.7,114.1,113.7,113.3,
            112.4,112.3,111.5,110.7,110.2,109.7,108.9,108.5,107.6,106.7,106.3,105.4,104.8,104.1,103.8,103.0,
            102.7,102.2,101.5,101.0,100.7,100.1,99.7,99.4,99.3,98.8,98.3,97.8,97.6,97.3,96.5,96.2,96.0,95.5,
            95.3,94.4,94.2,93.7,93.7,93.1,92.6,92.3,91.8,91.6,91.3,91.0,90.8,90.5,90.3,90.3,90.3,89.8,90.1,
            89.8,89.7,89.5,89.5,89.2,89.2,89.2,88.8,88.5,88.5,88.0,87.9,87.5,87.2,86.9,86.9,86.6,86.4,86.1,
            85.7,85.7,85.4,85.4,85.4,84.9,84.4,84.4,84.3,84.1,83.6,83.6,83.3,83.1,83.1)
samplingRate<-200 # sampling rate=200Hz (in other examples it might be 250Hz or ...)
tim<-seq(0,by=1/samplingRate,length=length(pressure)) # times in seconds
waveform<-data.frame(pressure,tim)
## beat-specific pressure waveform decay only ...
maxp<-max(waveform$pressure)
tmaxp<-waveform$tim[waveform$pressure==maxp]
tmaxp<-tmaxp[length(tmaxp)] # if time of max is not unique, use the last time
decay<-waveform[waveform$tim>=tmaxp,] # waveform decay only
decay$tim<-decay$tim-tmaxp # means P(0) is maximum pressure
## the following 2 lines, which calibrate the decay to [0,1], are not necessary
# the resulting PTC1 and PTC2 are the same with and without calibration
# we calibrate so it is easier to check for poor fits, a1exp(-a2t)~0 or a3exp(-a4t)cos(a5t+a6)~0
calibration<-summary(lm(c(1,0)~c(maxp,min(decay$pressure))))$coef
decay$pressure<-calibration[1]+calibration[2]*decay$pressure
## Fit Windkessel model to pressure decay and obtain beat-specific PTC1 and PTC2
# use nonlinear least squares (nls function), take advantage of conditional linearity of model
# use a grid of starting values
startval<-expand.grid(a2=c(0.5,1.5,2.5,3.5), a4=c(1,3,5,7,9), a5=c(10,20,30), a6=c(pi/4,pi/2,pi*3/4) )
df<-data.frame(a2=numeric(0), a4=numeric(0), a5=numeric(0), sigma=numeric(0))
for (k in 1:length(startval$a2))  
{
  fit<-try(nls(pressure~cbind(1,exp(-a2*tim),exp(-a4*tim)*cos(a5*tim+a6)),data=decay,
               start=startval[k,],  
               algorithm="plinear", trace=F, nls.control(tol=1e-5)),silent=T)
  if (!(attr(fit,"class")=="try-error")) df[k,]<-c(coef(fit)[1:3], summary(fit)$sigma)
  # use the try function and above line to avoid when errors or kth set of starting values fails to converge
}
df<-df[!is.na(df$sigma),]
result<-df[df$sigma==min((df$sigma)),][1,] # if converges for more than 1 set of starting values, use best fit
while(abs(result$a5)>2*pi/0.005) result$a5<-abs(result$a5)-samplingRate*2*pi # |a5|<200*2*pi for identifiability
ptc2<-with(result, 1/(2*a4+a2)*1000) # multiplied by 1000 so PTC2 in milliseconds
ptc1<-with(result, 2*a4*( (a2+a4)^2 + a5^2)  * ptc2/( a2*(a4^2+a5^2) )) #PTC1 in milliseconds
c(ptc1,ptc2) # PTC1 and PTC2, in milliseconds
