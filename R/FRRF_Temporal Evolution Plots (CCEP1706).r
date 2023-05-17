source('source.r')
source('source.frrf.r')

get.ship.val = function(time, ship.data, col = 'FL', delta.t = 10) {
    dt = as.numeric(difftime(ship.data$DT, time, units = 'mins'))
    l = which(dt^2 < delta.t^2)
    
    ## Return
    mean(ship.data[l, col], na.rm = TRUE)
}

get.ship.val.closest = function(time, ship.data, col = 'Cycle') {
    dt = as.numeric(difftime(ship.data$DT, time, units = 'mins'))
    l = which.min(dt^2)
    
    ## Return
    ship.data[l, col]
}

input.dir = '../CCE-shipdata/MIMS-TBK/Raw Data/FRRF/'
frrf.files = list.files(input.dir, pattern = '*.csv', full.names = FALSE)
length(frrf.files)

load('../CCE-shipdata/MIMS-TBK/RStates/ship.ncp.2min.rdata')

data = load.frrf(input.dir, frrf.files)

## Models
model = function(A, I, alpha, beta) {    
    A * (1 - exp(-alpha/A * I)) * exp(-beta/A * I)  ## Platt 1980
}

get.platt.fit = function(E, signal, N = 10) {
    
    ## Filter for NA values and for low light (fit is not suited for E = ~0).
    k = which(!is.na(E) & !is.na(signal) & E > 10)
    data = data.frame(signal = signal[k], E = E[k])
    
    m.init = as.numeric(coefficients(lm(signal[1:4] ~ E[1:4]))[2])
    
    ## Claculate general fit
    fit = nls(signal ~ A * (1 - exp(-alpha / A * E)) * exp(-beta / A * E),  
              control = nls.control(maxiter = 100, warnOnly = TRUE),
              data = data,
              algorithm = 'port',
              start = list(A = quantile(data$signal, probs = 0.95), alpha = m.init, beta = m.init/10),
              lower = list(A = median(data$signal, na.rm = TRUE), alpha = m.init/10, beta = 0),
              upper = list(A = max(data$signal, na.rm = TRUE), alpha = 10 * m.init, beta = m.init * 2))
    
    ## Initialize fits (for bootstrap) and df for storing fit parameters
    fits = fit
    
    ## Return fit entry
    list(fit = fit, fit.parameters = fit$m$getPars(), bootstrap = NA, data = data)
}


get.piecewise.fit = function(E, signal, x = seq(1, 4e3, 10)) {
    l = which(!is.na(E) & !is.na(signal))
    pred = data.frame(x = x,
                      y = approx(x = E[l], y = signal[l], xout = x, rule = 2))
    
    fit = list(prediction = pred,
               prediction.f = approx(x = E[l], y = signal[l], xout = x, rule = 2))
    
    ## Return fit entry
    fit
}

## use the boostrap above to generate a look-up table of values.
add.prediction = function(fit, x = seq(1, 4e3, 10)) {
    
    temp = model(A = fit$fit.parameters['A'],
                 I = x,
                 alpha = fit$fit.parameters['alpha'],
                 beta = fit$fit.parameters['beta'])
    
    ## Lookup table
    pred = data.frame(x = x, y = temp)
    fit$prediction = pred
    
    ## Return prediction information
    fit
}

#for (i in 1:5) {
for (i in 1:length(data)) {
    data[[i]]$JVPII = get.platt.fit(E = data[[i]]$A$E, signal = data[[i]]$A$JVPII, N = 10)
    data[[i]]$JVPII = add.prediction(data[[i]]$JVPII)
    data[[i]]$NPQ = get.piecewise.fit(E = data[[i]]$A$E, signal = data[[i]]$A$NSV)
}

l1 = which(ship.data$Cycle == 'Cycle1')
l2 = which(ship.data$Cycle == 'Cycle2')
l3 = which(ship.data$Cycle == 'Cycle3')
l4 = which(ship.data$Cycle == 'Cycle4')

add.frrf.jvpii = function() {
    for (i in 1:length(data)) {
        points(conv.time.unix(data[[i]]$Datetime), data[[i]]$JVPII$fit.parameters['alpha']*1e3, pch = 20)
    }
}

add.frrf.jvpii2 = function() {
    for (i in 1:length(data)) {
        points(conv.time.unix(data[[i]]$Datetime), data[[i]]$JVPII$fit.parameters['beta']*1e3, pch = 20)
    }
}

pdf('../Images/FRRF_Evolution Plots (alpha beta).pdf')

par(mfrow=c(3,1))

plot(ship.data$DT[l1], ship.data$PA[l1] / 250, type = 'l', xlab = '', ylab = 'Alpha (x 10-3)', col = 'grey', main = 'Cycle 1')
add.frrf.jvpii()

plot(ship.data$DT[l2], ship.data$PA[l2] / 250, type = 'l', xlab = '', ylab = 'Alpha (x 10-3)', col = 'grey', main = 'Cycle 2')
add.frrf.jvpii()

plot(ship.data$DT[l3], ship.data$PA[l3] / 250, type = 'l', xlab = '', ylab = 'Alpha (x 10-3)', col = 'grey', main = 'Cycle 3')
add.frrf.jvpii()

plot(ship.data$DT[l4], ship.data$PA[l4] / 250, type = 'l', xlab = '', ylab = 'Alpha (x 10-3)', col = 'grey', main = 'Cycle 4')
add.frrf.jvpii()


plot(ship.data$DT[l1], ship.data$PA[l1] / 4000, type = 'l', xlab = '', ylab = 'Beta  (x 10-3)', col = 'grey', main = 'Cycle 1')
add.frrf.jvpii2()

plot(ship.data$DT[l2], ship.data$PA[l2] / 4000, type = 'l', xlab = '', ylab = 'Beta  (x 10-3)', col = 'grey', main = 'Cycle 2')
add.frrf.jvpii2()

plot(ship.data$DT[l3], ship.data$PA[l3] / 4000, type = 'l', xlab = '', ylab = 'Beta  (x 10-3)', col = 'grey', main = 'Cycle 3')
add.frrf.jvpii2()

plot(ship.data$DT[l4], ship.data$PA[l4] / 4000, type = 'l', xlab = '', ylab = 'Beta  (x 10-3)', col = 'grey', main = 'Cycle 4')
add.frrf.jvpii2()

dev.off()

add.fo.sigma = function() {
    for (i in 1:length(data)) {
        points(conv.time.unix(data[[i]]$Datetime), data[[i]]$A$Fo[1] / data[[i]]$A$Sigma[1], pch = 20)
    }
}

pdf('../Images/FRRF_Evolution Plots (Fo sigma).pdf')

par(mfrow=c(3,1))

plot(ship.data$DT[l1], ship.data$PA[l1] / 600, type = 'l', xlab = '', ylab = 'Fo/Sigma', col = 'grey', main = 'Cycle 1')
add.fo.sigma()

plot(ship.data$DT[l2], ship.data$PA[l2] / 600, type = 'l', xlab = '', ylab = 'Fo/Sigma', col = 'grey', main = 'Cycle 2')
add.fo.sigma()

plot(ship.data$DT[l3], ship.data$PA[l3] / 600, type = 'l', xlab = '', ylab = 'Fo/Sigma', col = 'grey', main = 'Cycle 3')
add.fo.sigma()

plot(ship.data$DT[l4], ship.data$PA[l4] / 600, type = 'l', xlab = '', ylab = 'Fo/Sigma', col = 'grey', main = 'Cycle 4')
add.fo.sigma()

dev.off()


