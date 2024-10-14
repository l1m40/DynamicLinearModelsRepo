# File:   
#
# 
# 
# 
# Inspired on: 
# https://kevinkotze.github.io/ts-4-tut/
# https://hedibert.org/wp-content/uploads/2015/03/EconometriaAvancada-aula7.pdf
# 
# INSTALL AND LOAD PACKAGES ################################

# Installs pacman ("package manager") if needed
if (!require("pacman")) install.packages("pacman")

# Use pacman to load add-on packages as desired
# Packages I load every time; uses "pacman"
pacman::p_load(pacman,tidyverse,tidyquant,lubridate,gridExtra,dlm) 



# Current worklog ##########################################









# Price forecast with Kalman Filter ||| First approach #####
price_forecast_not_estimated <- function(Y_df,period_input,build_input){
  result_df <- NULL
  for(d in (period_input+5):nrow(Y_df)){
    train_df <- Y_df %>% slice((d-period_input-0):(d-1)) #%>% mutate(time_exp=log(300+row_number()))
    test_df  <- Y_df %>% slice(d) 
    gradient <- NA
    dV=NA;dW1=NA;dW2=NA;NVR1=NA;NVR2=NA
    # How to estimate variance of the observation noise?
    # fit <- dlmMLE(y=train_df$price,parm=rep(0,3),build=build,method="BFGS",control=list(trace=3,maxit=50))
    # gradient=as.numeric(regmatches(fit$counts[2],regexpr("\\d+",fit$counts[2])))
    # dV=exp(fit$par[1])
    # dW1=exp(fit$par[2])
    # dW2=exp(fit$par[3])
    # NVR1=exp(fit$par[2])/exp(fit$par[1])
    # NVR2=exp(fit$par[3])/exp(fit$par[1])
    # mod <- build_input(fit$par)
    # best fit for all? .. dV = 0.1 .. dW1 = 0.1 .. dW2 = 0.000001
    mod <- build_input(c(log(.1),log(.1),log(.00000000001)))
    filtered <- dlmFilter(train_df$price, mod)
    sd=residuals(filtered)$sd
    sd_residuals=sd[length(sd)]
    # sd_y_f=sd(filtered$y-filtered$f)
    m=filtered$m[period_input+1]
    ci=sqrt(dlmSvd2var(filtered$U.C[[period_input+1]],filtered$D.C[period_input+1,])[1,1])
    # paste("NVR[1]",diag(mod$W)[1]/mod$V,"  NVR[2]",diag(mod$W)[2]/mod$V) # just to double check if NVR of fit is the same
    # smoothed <- dlmSmooth(filtered)
    # smoothed$s
    forecast <- dlmForecast(filtered,n=1)$f[1,]
    result_df <- result_df %>% rbind(test_df %>% cbind(data.frame(gradient,dV,dW1,dW2,NVR1,NVR2,sd_residuals,ci,m,forecast))) %>% as_tibble()
  }
  return(result_df)
}

Y_df <- tq_get("ITUB4.SA",from=today()-300,to=today()+1,get="stock.prices") %>% mutate(price=adjusted,asset=str_remove(symbol,".SA")) %>% mutate(adjust_gap=close-adjusted,open=open-adjust_gap,high=high-adjust_gap,low=low-adjust_gap)
Y_df <- Y_df %>% arrange(date) %>% mutate(log_returns=log(price/lag(price,5))) %>% select(date,open,high,low,price,log_returns) %>% mutate(price_mean_2=rollapply(price,list(-seq(1:2)),mean,fill=NA)) %>% tail(300) %>% mutate(row_no=row_number())
build <- function(params) { return(dlmModPoly(dV=exp(params[1]),dW=c(exp(params[2:3])))) }
result_df <- price_forecast_not_estimated(Y_df,80,build)
result_df %>% summary

result_df %>% ggplot(aes(x=date))+#geom_line(aes(y=dW1))
  geom_ribbon(aes(ymax=(forecast+sd_residuals*qnorm(0.95)),ymin=(forecast-sd_residuals*qnorm(0.95))),alpha=.1,fill="cyan3")+
  #geom_ribbon(aes(ymax=(forecast+ci*qnorm(0.75)),ymin=(forecast-ci*qnorm(0.75))),alpha=.1)+
  geom_line(data=Y_df,aes(y=price),color="gray")+
  geom_line(aes(y=price   ),color="black")+geom_point(aes(y=price   ),color="black",size=.5)+
  geom_line(aes(y=forecast),color="cyan3")+geom_point(aes(y=forecast),color="cyan3",size=.5)+
  #geom_point(data=result_df %>% mutate(forecast=ifelse(NVR1<100,forecast,NA)),aes(y=forecast),size=3,color="green3",alpha=.5)+
  theme_minimal()





# Price forecast with Kalman Filter ||| MLE ################
price_forecast_mle <- function(Y_df,period_input,build_input){
  result_df <- NULL
  for(d in (period_input+5):nrow(Y_df)){
    train_df <- Y_df %>% slice((d-period_input-0):(d-1)) #%>% mutate(time_exp=log(300+row_number()))
    test_df  <- Y_df %>% slice(d) 
    gradient <- NA
    dV=NA;dW1=NA;dW2=NA;NVR1=NA;NVR2=NA
    # How to estimate variance of the observation noise?
    fit_calc <- tryCatch({dlmMLE(y=train_df$price,parm=rep(0,3),build=build,method="BFGS",control=list(trace=3,maxit=50))},warning=function(e){print(e)},error=function(e){print(e)},finally=NA)
    if(exists("par",fit_calc)) fit <- fit_calc
    gradient=as.numeric(regmatches(fit$counts[2],regexpr("\\d+",fit$counts[2])))
    dV=exp(fit$par[1])
    dW1=exp(fit$par[2])
    dW2=exp(fit$par[3])
    NVR1=exp(fit$par[2])/exp(fit$par[1])
    NVR2=exp(fit$par[3])/exp(fit$par[1])
    mod <- build_input(fit$par)
    filtered <- dlmFilter(train_df$price,mod)
    sd=residuals(filtered)$sd
    sd_residuals=sd[length(sd)]
    # sd_y_f=sd(filtered$y-filtered$f)
    m=filtered$m[period_input+1]
    ci=sqrt(dlmSvd2var(filtered$U.C[[period_input+1]],filtered$D.C[period_input+1,])[1,1])
    # paste("NVR[1]",diag(mod$W)[1]/mod$V,"  NVR[2]",diag(mod$W)[2]/mod$V) # just to double check if NVR of fit is the same
    # smoothed <- dlmSmooth(filtered)
    # smoothed$s
    forecast <- dlmForecast(filtered,n=1)$f[1,]
    result_df <- result_df %>% rbind(test_df %>% cbind(data.frame(gradient,dV,dW1,dW2,NVR1,NVR2,sd_residuals,ci,m,forecast))) %>% as_tibble()
  }
  return(result_df)
}

Y_df <- tq_get("ITUB4.SA",from=today()-300,to=today()+1,get="stock.prices") %>% mutate(price=adjusted,asset=str_remove(symbol,".SA")) %>% mutate(adjust_gap=close-adjusted,open=open-adjust_gap,high=high-adjust_gap,low=low-adjust_gap)
Y_df <- Y_df %>% arrange(date) %>% mutate(log_returns=log(price/lag(price,5))) %>% select(date,open,high,low,price,log_returns) %>% mutate(price_mean_2=rollapply(price,list(-seq(1:2)),mean,fill=NA)) %>% tail(300) %>% mutate(row_no=row_number())
build <- function(params) { return(dlmModPoly(dV=exp(params[1]),dW=c(exp(params[2:3])))) }
result_basic_df <- price_forecast_not_estimated(Y_df,80,build)
result_mle_df   <- price_forecast_mle(Y_df,80,build)

NULL %>% 
  rbind(result_basic_df %>% summarise(error=mean((forecast/price)-1,na.rm=T)) %>% mutate(test="basic estimate")) %>% 
  rbind(result_mle_df   %>% summarise(error=mean((forecast/price)-1,na.rm=T)) %>% mutate(test="maximum likelihood"))

result_mle_df %>% ggplot(aes(x=date))+#geom_line(aes(y=dW1))
  geom_ribbon(data=result_basic_df,aes(ymax=(forecast+sd_residuals*qnorm(0.95)),ymin=(forecast-sd_residuals*qnorm(0.95))),alpha=.1,fill="gray")+
  geom_ribbon(aes(ymax=(forecast+sd_residuals*qnorm(0.95)),ymin=(forecast-sd_residuals*qnorm(0.95))),alpha=.1,fill="cyan3")+
  #geom_ribbon(aes(ymax=(forecast+ci*qnorm(0.75)),ymin=(forecast-ci*qnorm(0.75))),alpha=.1)+
  geom_line(data=Y_df,aes(y=price),color="gray")+
  geom_line(aes(y=price   ),color="black")+geom_point(aes(y=price   ),color="black",size=.5)+
  geom_line(aes(y=forecast),color="cyan3")+geom_point(aes(y=forecast),color="cyan3",size=.5)+
  #geom_point(data=result_df %>% mutate(forecast=ifelse(NVR1<100,forecast,NA)),aes(y=forecast),size=3,color="green3",alpha=.5)+
  theme_minimal()












# Use typical price ########################################
Y_df <- tq_get("ITUB4.SA",from=today()-300,to=today()+1,get="stock.prices") %>% mutate(price=adjusted,asset=str_remove(symbol,".SA")) %>% mutate(adjust_gap=close-adjusted,open=open-adjust_gap,high=high-adjust_gap,low=low-adjust_gap)
Y_df <- Y_df %>% arrange(date) %>% mutate(log_returns=log(price/lag(price,5))) %>% select(date,open,high,low,price,log_returns) %>% mutate(price_mean_2=rollapply(price,list(-seq(1:2)),mean,fill=NA)) %>% tail(300) %>% mutate(row_no=row_number())
build <- function(params) { return(dlmModPoly(dV=exp(params[1]),dW=c(exp(params[2:3])))) }
result_basic_df <- price_forecast_not_estimated(Y_df,80,build)
result_mle_1_df <- price_forecast_mle(Y_df,80,build)
result_mle_df   <- price_forecast_mle(Y_df %>% mutate(price=(price+high+low)/3),80,build)

NULL %>% 
  rbind(result_basic_df %>% mutate(error=abs((forecast/price)-1)) %>% mutate(err_mean=mean(error,na.rm=T),err_max=max(error,na.rm=T)) %>% mutate(test="  not estimated")) %>% 
  rbind(result_mle_1_df %>% mutate(error=abs((forecast/price)-1)) %>% mutate(err_mean=mean(error,na.rm=T),err_max=max(error,na.rm=T)) %>% mutate(test=" maximum likelihood")) %>% 
  rbind(result_mle_df   %>% mutate(error=abs((forecast/price)-1)) %>% mutate(err_mean=mean(error,na.rm=T),err_max=max(error,na.rm=T)) %>% mutate(test="mle with typical price")) %>% 
  ggplot(aes(y=test))+ylab("Model")+xlab("Forecast Error Average")+geom_boxplot(aes(x=error))
  #geom_col(aes(x=err_mean))+geom_point(aes(x=err_max),shape=8,color="red")+geom_text(aes(x=err_max,label="Error max."),color="red",hjust=1.3)

result_mle_df %>% ggplot(aes(x=date))+#geom_line(aes(y=dW1))
  geom_ribbon(data=result_basic_df,aes(ymax=(forecast+sd_residuals*qnorm(0.95)),ymin=(forecast-sd_residuals*qnorm(0.95))),alpha=.1,fill="gray")+
  geom_ribbon(aes(ymax=(forecast+sd_residuals*qnorm(0.95)),ymin=(forecast-sd_residuals*qnorm(0.95))),alpha=.1,fill="cyan3")+
  #geom_ribbon(aes(ymax=(forecast+ci*qnorm(0.75)),ymin=(forecast-ci*qnorm(0.75))),alpha=.1)+
  geom_line(data=Y_df,aes(y=price),color="gray")+
  geom_line(aes(y=price   ),color="black")+geom_point(aes(y=price   ),color="black",size=.5)+
  geom_line(aes(y=forecast),color="cyan3")+geom_point(aes(y=forecast),color="cyan3",size=.5)+
  #geom_point(data=result_df %>% mutate(forecast=ifelse(NVR1<100,forecast,NA)),aes(y=forecast),size=3,color="green3",alpha=.5)+
  theme_minimal()






