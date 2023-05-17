library(TheSource)
library(openxlsx)
library(data.table)
options(warn = -1)

#### A utility script that will pull out the datetime of all frrf scans. Makes looking up times easier
get.frrf.times = function(data) {
    times = c()

    for (i in 1:length(data)) {
        times = c(times, data[[i]]$Datetime)
    }

    conv.time.unix(times)
}

#### Extract frrf scans based on time of day
extract.dark = function(data, hours = c(7, 11)) {
    result = list()
    
    for (i in 1:length(data)) {
        h = get.hour(conv.time.unix(data[[i]]$Datetime))
        if (h >= hours[1] & h <= hours[2]) {
            result[[names(data)[i]]] = data[[i]]
        }
    }
    
    result
}

#### Extract scans based on Cycle
extract.cycle = function(data, cycle) {
    result = list()
    
    for (i in 1:length(data)) {
        t = conv.time.unix(data[[i]]$Datetime, tz = 'GMT')
        
        if (t >= cycle[1] & t <= cycle[2]) {
            result[[names(data)[i]]] = data[[i]]
        }
    }
    
    result
}

#### A quick summary of data, should add any and all parameters of interest.
## Useful for diagnostics
get.summary = function(data) {
    summary = data.frame(Time = rep(NA, length(data)), Cycle = NA, f = NA, RCII = NA, Chl = NA,
                         JVPII.250 = NA, JVPII.250.abs = NA, JVPII.250.sigma = NA, ETR.250 = NA,
                         Pm = NA, alpha = NA, beta = NA, Fo = NA, Sigma =NA)
    
    for (i in 1:length(data)) {
        summary$Time[i] = data[[i]]$Datetime
        summary$Cycle[i] = data[[i]]$Params$Cycle
        summary$RCII[i] = data[[i]]$Params$RCII
        summary$f[i] = data[[i]]$Params$f
        summary$Fo[i] = data[[i]]$A$Fo[1]
        summary$Sigma[i] = data[[i]]$A$Sigma[1]
        summary$Chl[i] = data[[i]]$Params$Chl
        summary$JVPII.250[i] = approx(data[[i]]$A$E, data[[i]]$A$JVPII, 250)$y
        summary$JVPII.250.abs[i] = approx(data[[i]]$A$E, data[[i]]$A$JVPII.abs, 250)$y
        summary$JVPII.250.sigma[i] = approx(data[[i]]$A$E, data[[i]]$A$JVPII.sigma, 250)$y
        summary$ETR.250[i] = approx(data[[i]]$A$E, data[[i]]$A$ETR, 250)$y
        summary$Pm[i] = data[[i]]$Fit$p$min$c
        summary$alpha[i] = data[[i]]$Fit$p$min$a
        summary$beta[i] = data[[i]]$Fit$p$min$b
    }
    
    summary
}

## Do this first, we will be using only dark rcii later on.
calc.frrf.rcii = function(scan) { # [mol m-3]
    # [    ]     [nm-2]       [nm2 m-2]    [    ]      [mol-1]
    scan$A$Fo[1] / scan$A$Sigma[1] / 1e-18 * scan$Params$Ka / 6.022e23 ## Oxborough: [] [m-2] [m-1] [mol]
}


calc.frrf.Fo.p = function(scan) {
    Fq.Fm = (scan$A$Fm - scan$A$Fo) / scan$A$Fm
    
    scan$A$Fm[1] * scan$A$Fo[1] / (scan$A$Fm[1] - scan$A$Fo[1]) * Fq.Fm
}


calc.frrf.FqFv = function(scan) {
    Fo.p = calc.frrf.Fo.p(scan)
    
    (scan$A$Fm - scan$A$Fo) / (scan$A$Fm - Fo.p)
}


calc.frrf.jvpii.abs = function(scan) { # [mol e m-3 d-1]
    Fo.p = calc.frrf.Fo.p(scan)
    
    ## [umol photons m-2 s-1] [s d-1] [mol/umol] [kg?] [m-1] = [mol photons m-3 d-1] = [mol electrons m-3 d-1] [mg?]
    scan$A$E * 86400 * 1e-6 * Fo.p * scan$Params$Ka * 1e-6
}


calc.frrf.jvpii.sigma = function(scan) {
    RCII = scan$Params$RCII * 6.022e23 # RCII
    
    # [umol photons s-1] [mol umol-1] [s d-1] [nm2 RCII-1] [m2 nm-2] [unitless] [RCII] [?]
    scan$A$E * 1e-6 * 86400 * scan$A$Sigma * 1e-18 * (1 - scan$A$C) * RCII * 1e-6
}


calc.frrf.etr = function(scan) { ## Schuback 2017: [mol e m-3 d-1]
    Fo = scan$A$Fo[1] / (scan$A$Fv.Fm[1] + scan$A$Fo[1] / scan$A$Fm)
    qP = (scan$A$Fm - scan$A$Fo) / (scan$A$Fm - Fo)
    rcii = scan$Params$RCII * 6.022e23 # [RCII]
    
    ## [umol photons m-2 d-1] [mol umol-1] * [nm2 RCII-1] * qP * [RCII (mol RCII)-1 * m2 nm-2] * [s d-1] * [mol RCII]
    #scan$A$E * 1e-6 * scan$A$Sigma * qP * 6.022e-1 * 86400 * rcii
    # [ mol m-2 s-1 ] [nm2 RCII-1]   [] [m2 nm-2] [s d-1] [RCII] [?]
    scan$A$E * 1e-6 * scan$A$Sigma * qP * 1e-18 * 86400 * rcii * 1e-6
}


calc.frrf.psi.ec = function(scan) { # Schuback
    ## [ emprirical formula] * [mol RCII] * [(mol chl)-1]
    #(486 * scan$A$NSV + 1854) * scan$Params$RCII / scan$Params$Chl
    (486 * scan$A$NSV + 1854) * 0.003
}

## Load Chl data (nighttime only)
chl = read.xlsx('../Documents/CCE/Cruise_CCEP1706/Chlorophyll.xlsx')
chl$Datetime.GMT = conv.time.excel(chl$Datetime.GMT)
chl$Hour = get.hour(chl$Datetime.GMT)
chl = chl[chl$Hour >= 8 & chl$Hour <= 10 & chl$Depth < 10, ]

## Load MLD data
load('../CCE-shipdata/MIMS-TBK/RStates/mld-0.03.rdata')
mld = mld[mld$Cast < 1000,]

## Load event Log
events = read.xlsx('../CCE-shipdata/MIMS-TBK/Raw Data/CCEP1706 Event Log.xlsx')
events$Datetime.UTC = conv.time.excel(events$Datetime.UTC)
events = events[!is.na(events$Datetime.UTC),]

## Load NCP Data
#load('../CCE-shipdata/MIMS-TBK/RStates/ship.ncp.30min.rdata')
load('../Documents/CCE/NCP1706/RStates/ship.ncp.15.rdata')
ship.data = ship.data[ship.data$Cycle != 'Cycle0' & ship.data$Day > 0 & ship.data$Dist.to.Drifter < 8,]

par(mfrow=c(2,2))

chl$ship = NA

for (i in 1:nrow(chl)) {
    l = which(as.numeric(difftime(ship.data$DT, chl$Datetime.GMT[i], units = 'mins'))^2 < 20^2)
    chl$ship[i] = mean(ship.data$FL[l], na.rm = TRUE)
}

ols = lm(chl$Chl ~ chl$ship + 0)
summary(ols)

plot(chl$ship, chl$Chl)
abline(ols)

ship.data$FL = ship.data$FL * coef(ols)[1]
ship.data$FL[ship.data$FL < min(chl$Chl, na.rm = TRUE)] = min(chl$Chl, na.rm = TRUE)

plot(ship.data$DT, ship.data$FL)

nitrate = read.xlsx('../Documents/CCE/Cruise_CCEP1706/Nitrate Uptake.xlsx')

nitrate$Time = conv.time.excel(nitrate$Time)
nitrate$mld = approx(mld$time, mld$mld, nitrate$Time, rule = 2)$y
nitrate$Conser.up.uM[is.na(nitrate$Conser.up.uM)] = -10

nitrate.int = data.frame(Time = unique(nitrate$Time))
nitrate.int$int = NA
nitrate.int$int.ml = NA

for (i in 1:nrow(nitrate.int)) {
    l = which(nitrate$Time == nitrate.int$Time[i])
    
    profile = approx(nitrate$Depth[l], nitrate$Conser.up.uM[l], c(1:100), rule = 2)$y
    
    nitrate.int$int[i] = sum(profile, na.rm = TRUE)
    nitrate.int$int.ml[i] = sum(profile[1:nitrate$mld[l[1]]], na.rm = TRUE)
    
    nitrate.int$MLD[i] = nitrate$mld[l[1]]
}

nitrate$Conser.up.uM[nitrate$Conser.up.uM == -10] = NA

npp = read.xlsx('../Documents/CCE/NCP Manuscripts/P2 - Kranz et al/Datafiles/NPP/NPP Product from Stukel.xlsx', startRow = 2)
npp$Time = conv.time.excel(npp$Time+0.5, tz = 'GMT')
npp[,c(6:8)] = npp[,c(6:8)] * npp$MLD / 12
npp[,c(10:12)] = npp[,c(10:12)] * npp$MLD2 / 12

str(npp)

background = data.frame(Time = seq(as.POSIXct("2017-06-09 07:00:00", tz = 'GMT'), P1706.4()[2], by = '1 hour'),
                        Cycle = NA, NPP.ml = NA, Nitrate.uptake = NA, Nitrate.uptake.ml = NA,
                        Lat = NA, Lon = NA, MLD = NA, T = NA, S = NA, Chl = NA, Chl.ship = NA, PAR = NA, nPSII = NA,
                        NCP.prior.eims = NA, NCP.RT.eims = NA, NCP.NSS.eims = NA,
                        NCP.prior.mims = NA, NCP.RT.mims = NA, NCP.NSS.mims = NA)

background$MLD = approx(mld$time, mld$mld, background$Time, rule = 2)$y
background$Chl = approx(chl$Datetime.GMT, chl$Chl, background$Time, rule = 2)$y
background$Chl.ship = approx(ship.data$DT, ship.data$FL, background$Time, rule = 2)$y

background$Nitrate.uptake = approx(nitrate.int$Time, nitrate.int$int, background$Time, rule = 2, method = 'constant')$y
background$Nitrate.uptake.ml = approx(nitrate.int$Time, nitrate.int$int.ml, background$Time, rule = 2, method = 'constant')$y

background$Nitrate.uptake[background$Nitrate.uptake < 0] = NA
background$Nitrate.uptake.ml[background$Nitrate.uptake.ml < 0] = NA

background$Cycle = ship.data$Cycle[which.closest.time(background$Time, ship.data$DT)]
background$Day = ship.data$Day[which.closest.time(background$Time, ship.data$DT)]

for (i in 1:nrow(background)) {
    
    l.npp = which(npp$Cycle == background$Cycle[i] & npp$Day == background$Day[i])
    background$NPP.ml[i] = npp$NPP[l.npp]
    
    l = which(as.numeric(difftime(ship.data$DT, background$Time[i], units = 'mins'))^2 <= 30^2)
    background$T[i] =   mean(ship.data$TT[l], na.rm = TRUE)
    background$S[i] =   mean(ship.data$SA[l], na.rm = TRUE)
    background$PAR[i] = mean(ship.data$PA[l], na.rm = TRUE)
    background$Lat[i] = mean(ship.data$LA[l], na.rm = TRUE)
    background$Lon[i] = mean(ship.data$LO[l], na.rm = TRUE)
    
    background$NCP.prior.eims[i] = mean(ship.data$NCP.EIMS[l], na.rm = TRUE)
    background$NCP.RT.eims[i] =    mean(ship.data$NCP.EIMS.RT[l], na.rm = TRUE)
    background$NCP.NSS.eims[i] =   mean(ship.data$NCP.EIMS.NSS[l], na.rm = TRUE)
    
    background$NCP.prior.mims[i] = mean(ship.data$NCP.MIMS[l], na.rm = TRUE)
    background$NCP.RT.mims[i] =    mean(ship.data$NCP.MIMS.RT[l], na.rm = TRUE)
    background$NCP.NSS.mims[i] =   mean(ship.data$NCP.MIMS.NSS[l], na.rm = TRUE)
}

par(mfrow = c(2,2))

plot(background$Time, background$NCP.prior.eims, type = 'l')
plot(background$Time, background$NCP.NSS.eims, type= 'l')
plot(background$Time, background$NCP.prior.mims, type = 'l')
plot(background$Time, background$NCP.RT.mims, type = 'l')

input.dir = '../CCE-shipdata/MIMS-TBK/Raw Data/FRRF/'
frrf.files = list.files(input.dir, pattern = '*.csv', full.names = FALSE)
length(frrf.files)

data = load.frrf(input.dir, frrf.files)

print(paste0('Number of raw entries: ', length(data)))
f = approxfun(ship.data$DT, ship.data$Dist.to.Drifter, rule = 2)

for (i in 1:length(data)) {
    time = conv.time.unix(data[[i]]$Datetime)
    
    if (time >= P1706.1()[1] & time <= P1706.1()[2]) {
        data[[i]]$Params$Cycle = 'Cycle1'
    } else if (time >= P1706.2()[1] & time <= P1706.2()[2]) {
        data[[i]]$Params$Cycle = 'Cycle2'
    } else if (time >= P1706.3()[1] & time <= P1706.3()[2]) {
        data[[i]]$Params$Cycle = 'Cycle3'
    } else if (time >= P1706.4()[1] & time <= P1706.4()[2]) {
        data[[i]]$Params$Cycle = 'Cycle4'
    } else {
        data[[i]]$Params$Cycle = NA
    }
    
    data[[i]]$Params$Dist = f(time)
}

## Remove non-cycle measurements
i = 1
while (i <= length(data)) {
    
    check1 = is.na(data[[i]]$Params$Cycle)
    check2 = data[[i]]$Params$Dist > 4000
    
    if (check1 | check2) {
        data[[i]] = NULL
    } else {
         i = i+1
    }
}

print(paste0('Number of filtered entries: ', length(data)))

temp = extract.dark(data, hours = c(8,12))
rcii = data.frame(time = rep(NA, length(temp)), rcii = NA)

for (i in 1:length(temp)) {
    rcii$time[i] = temp[[i]]$Datetime
    rcii$rcii[i] = calc.frrf.rcii(temp[[i]])
}

f.chl = approxfun(ship.data$DT, ship.data$FL, rule = 2)

## Add dark RCII values to all scans!
for (i in 1:length(data)) {
    data[[i]]$Params$RCII = approx(rcii$time, runmed(rcii$rcii, 9, endrule = 'constant'), data[[i]]$Datetime, rule = 2)$y
    
    if (is.na(data[[i]]$Params$RCII)) {
        print(i)
    }
    
    #data[[i]]$Params$Chl = approx(chl$Datetime.GMT, chl$Chl,
    #                              conv.time.unix(data[[i]]$Datetime), rule = 2)$y / 893.51 * 1e3 / 2 ## [mol chl m-3]
    
    data[[i]]$Params$Chl = f.chl(conv.time.unix(data[[i]]$Datetime)) / 893.51 * 1e3 / 2 ## [mol chl m-3]
    
    #data[[i]]$Params$Chl = approx(background$Time,
    #                              runmed(background$Chl, 9, endrule = 'constant'),
    #                              conv.time.unix(data[[i]]$Datetime), rule = 2)$y / 893.51 * 1e3 / 2 ## [mol chl m-3]
    
    ## Set nPSII and use that to calculate the apparent productivity factor
    data[[i]]$Params$nPSII = 0.003
    data[[i]]$Params$f = 0.003 / (data[[i]]$Params$RCII / data[[i]]$Params$Chl)
}


for (i in 1:length(data)) {
    data[[i]]$A$JVPII.abs = calc.frrf.jvpii.abs(data[[i]])
    data[[i]]$A$JVPII.sigma = calc.frrf.jvpii.sigma(data[[i]])
    data[[i]]$A$ETR = calc.frrf.etr(data[[i]])
    data[[i]]$A$FqFv = calc.frrf.FqFv(data[[i]])
}

Webb1974 = function(alpha, Ek, E) {
    alpha * Ek * (1-exp(-E/Ek))
}

Platt1980 = function(alpha, beta, Ps, E) {
    Ps * (1 - exp(-1 * alpha * E / Ps)) * exp(-1 * beta * E / Ps)
}

Jassby1976 = function(alpha, Ek, E) {
    alpha * Ek * tanh(E / Ek)
}

Eilers1988 = function(alpha, Eopt, Pm, E) {
    a = 1 / (alpha * Eopt^2)
    b = 1 / Pm - 2 / (alpha * Eopt)
    c = 1/ alpha
    
    E / (a * E^2 + b * E + c)
}

#### calc.frrf.fit()
## scan: an individual frrf scan object.
## model: the model to be used (i.e. Platt1080 or Webb1974...)
## E = c(1:1e4): the light intensities you want calculated for the "prediction". Just
# leave this unless you have soemthing specific in mind.
## splits = 10: the number of parameter-grid search cuts to make each time. 
## type = 1: the type of JVPII you want to fit. 1 = Instrumental, 2 = Absorbance, and 3 = Sigma-method.
calc.frrf.fit = function(scan, model, E = c(0:1e4), splits = 10, type = 1) {
    
    if (type == 1) { jvpii = scan$A$JVPII}        # Instrument
    if (type == 2) { jvpii = scan$A$JVPII.abs}    # Abs Method
    if (type == 3) { jvpii = scan$A$JVPII.sigma}  # Sigma Method
    
    ## Remove NAs
    l = which(!is.na(jvpii))
    e.measured = scan$A$E[l]
    jvpii = jvpii[l]
    
    ## The general 'fit' object.
    fit = list(model = model, p = NA,
               ideal = data.frame(E = E, pred = NA),
               obs = data.frame(E = e.measured, Obs = jvpii, Pred = NA),
               type = type,
               model.name = deparse(substitute(model)))
        
    ## Next we go through each model possibility and if that is the specified model then
    # use that to fill out the fit object.
    if (deparse(substitute(model)) == 'Webb1974') {
        
        f = function(a, b, E, JVPII) {
            sum(abs(model(alpha = a, Ek = b, E = E) - JVPII)^2)
        }
        
        p = parameter.search(n = 5, cost = f, E = e.measured, JVPII = jvpii,
                bounds = data.frame(min = c(0.001, 100), max = c(0.01, 1000)), splits = splits)
        fit$p = p
        fit$ideal$pred = model(p$min$a[1], p$min$b[1], fit$ideal$E)
        fit$obs$Pred = model(p$min$a[1], p$min$b[1], fit$obs$E)
        
    }
    
    if (deparse(substitute(model)) == 'Platt1980') {
        
        f = function(a, b, c, E, JVPII) {
            sum(abs(model(alpha = a, beta = b, Ps = c, E = E) - JVPII)^2)
        }
        
        p = parameter.search(10, f, E = e.measured, JVPII = jvpii,
                 bounds = data.frame(min = c(0, 0, 0.2), max = c(0.05, 0.05, 15)), splits = splits)
        
        fit$p = p
        fit$ideal$pred = model(p$min[1,1], p$min[1,2], p$min[1,3], fit$ideal$E)
        fit$obs$Pred = model(p$min[1,1], p$min[1,2], p$min[1,3], fit$obs$E)
        
    }
    
    if (deparse(substitute(model)) == 'Jassby1976') {
        f = function(a, b, E, JVPII) {
            sum(abs(model(a, b, E = E) - JVPII)^2)
        }
        
        p = parameter.search(5, f, E = e.measured, JVPII = jvpii,
                 bounds = data.frame(min = c(0.001, 100), max = c(0.01, 1000)), splits = splits)
        
        fit$p = p
        fit$ideal$pred = model(p$min$a[1], p$min$b[1], fit$ideal$E)
        fit$obs$Pred = model(p$min$a[1], p$min$b[1], fit$obs$E)
    }
    
    if (deparse(substitute(model)) == 'Eilers1988') {
        f = function(a, b, c, E, JVPII) {
            sum(abs(model(a, b, c, E = E) - JVPII)^2)
        }
        
        p = parameter.search(5, f, E = e.measured, JVPII = jvpii,
                 bounds = data.frame(min = c(0, 500, 0), max = c(3, 100000, 500)), splits = splits)
        fit$p = p
        fit$ideal$pred = model(p$min$a[1], p$min$b[1], p$min$c[1], fit$ideal$E)
        fit$obs$Pred = model(p$min$a[1], p$min$b[1], p$min$c[1], fit$obs$E)
    }
    
    ## Return the fit object.
    fit
}

## Add a platt fit to each frrf scan.
platt = data.frame(Alpha = rep(NA, length(data)), Beta = NA, Pmax = NA, Cost = NA)

for (i in 1:length(data)) {
    fit.platt = calc.frrf.fit(data[[i]], Platt1980, splits = 10)
    data[[i]]$Fit = fit.platt
    
    platt$Alpha[i] = fit.platt$p$min$a
    platt$Beta[i] = fit.platt$p$min$b
    platt$Pmax[i] = fit.platt$p$min$c
    platt$Cost[i] = fit.platt$p$min$cost
}

par(mfrow=c(2,2))
plot(platt$Alpha, ylab = 'Alpha', xlab = '')
plot(platt$Beta, ylab = 'Beta', xlab = '')
plot(platt$Pmax, ylab = 'Pmax', xlab = '')
plot(platt$Cost, ylab = 'Cost', xlab = '')

par(mfrow = c(2,2))
for (i in which(platt$Cost > 2)) {
    plot(data[[i]]$A$E, data[[i]]$A$JVPII, main = paste('Scan', i), ylab = 'JVPII', xlab = 'E')
    lines(data[[i]]$Fit$ideal$E, data[[i]]$Fit$ideal$pred)
}

summary = get.summary(data)

par(mfrow=c(2,2))
plot(factor(summary$Cycle), summary$Chl, ylab = 'Chl')
plot(factor(summary$Cycle), summary$Fo, ylab = 'Fo')
plot(factor(summary$Cycle), summary$f, ylab = 'f')
plot(factor(summary$Cycle), summary$Sigma, ylab = 'Sigma')

plot(factor(summary$Cycle), summary$JVPII.250, ylab = 'JVPII (instrument)')
plot(factor(summary$Cycle), summary$JVPII.250.abs, ylab = 'JVPII (abs)')
plot(factor(summary$Cycle), summary$JVPII.250.sigma, ylab = 'JVPII (sigma)')
plot(factor(summary$Cycle), summary$ETR, ylab = 'ETR')

frrf.time = get.frrf.times(data)

for (i in 1:nrow(background)) {
    k = which.closest.time(frrf.time, background$Time[i])
    background$nPSII[i] = data[[k]]$Params$f
}



get.FvFm = function(data) {
    res = data.frame(time = rep(NA, length(data)), FvFm = NA)
    
    for (i in 1:length(data)) {
        res$time[i] = data[[i]]$Datetime
        res$FvFm[i] = data[[i]]$A$Fv.Fm[1]
    }
    
    res$time = conv.time.unix(res$time)
    res
}

get.Fo = function(data) {
    res = data.frame(time = rep(NA, length(data)), Fo = NA)
    
    for (i in 1:length(data)) {
        res$time[i] = data[[i]]$Datetime
        res$Fo[i] = data[[i]]$A$Fo[1]
    }
    
    res$time = conv.time.unix(res$time)
    res
}

get.Sigma = function(data) {
    res = data.frame(time = rep(NA, length(data)), Sigma = NA)
    
    for (i in 1:length(data)) {
        res$time[i] = data[[i]]$Datetime
        res$Sigma[i] = data[[i]]$A$Sigma[1]
    }
    
    res$time = conv.time.unix(res$time)
    res
}

get.Chl = function(data) {
    res = data.frame(time = rep(NA, length(data)), Chl = NA, Chl.frrf = NA)
    
    for (i in 1:length(data)) {
        res$time[i] = data[[i]]$Datetime
        res$Chl.frrf[i] = data[[i]]$Params$Chl.multiplier
        
        res$Chl[i] = data[[i]]$Params$Chl
    }
    
    res$time = conv.time.unix(res$time)
    res
}

fvfm = get.FvFm(data)

plot(fvfm$time, fvfm$FvFm, pch = 16, ylim = c(0,0.6), yaxs = 'i', ylab = 'Fv/Fm', xlab = '')

add.shade(P1706.1(), col = '#40004010')
add.shade(P1706.2(), col = '#40004010')
add.shade(P1706.3(), col = '#40004010')
add.shade(P1706.4(), col = '#40004010')

fo.dark = get.Fo(extract.dark(data, hours = c(7,13)))
fo.day = get.Fo(extract.dark(data, hours = c(18,22)))
fo = get.Fo(data)

sigma.dark = get.Sigma(extract.dark(data, hours = c(7,13)))
sigma.day = get.Sigma(extract.dark(data, hours = c(18,22)))
sigma = get.Sigma(data)

chl.dark = get.Chl(extract.dark(data, hours = c(7,13)))
chl.day = get.Chl(extract.dark(data, hours = c(18,22)))
chl.both = get.Chl(data)

par(mfrow=c(2,2))
plot(chl.dark$Chl, fo.dark$Fo/ sigma.dark$Sigma, pch = 20, ylim = c(0, 2.2), xlim = c(0, 6), ylab = 'Fo/Sigma',
     main = 'Night', xaxs='i', yaxs='i', xlab = 'Measured Chl')

plot(chl.day$Chl, fo.day$Fo/ sigma.day$Sigma, pch = 20, ylim = c(0, 2.2), xlim = c(0, 6), ylab = 'Fo/Sigma',
     main = 'Day', xaxs='i', yaxs='i', xlab = 'Measured Chl')

plot(chl.both$Chl, fo$Fo/ sigma$Sigma, pch = 20, ylim = c(0, 2.2), xlim = c(0, 6), ylab = 'Fo/Sigma',
     main = 'All', xaxs='i', yaxs='i', xlab = 'Measured Chl')

par(mfcol = c(3,4), mar = c(4,4,3,1))
for (i in 1:4) {
    if (i == 1) { temp = extract.cycle(data, P1706.1) }
    if (i == 2) { temp = extract.cycle(data, P1706.2) }
    if (i == 3) { temp = extract.cycle(data, P1706.3) }
    if (i == 4) { temp = extract.cycle(data, P1706.4) }
    
    fo.dark = get.Fo(extract.dark(temp, hours = c(6,13)))
    fo.day = get.Fo(extract.dark(temp, hours = c(18,22)))
    fo = get.Fo(temp)

    sigma.dark = get.Sigma(extract.dark(temp, hours = c(6,13)))
    sigma.day = get.Sigma(extract.dark(temp, hours = c(18,22)))
    sigma = get.Sigma(temp)

    chl.dark = get.Chl(extract.dark(temp, hours = c(6,13)))
    chl.day = get.Chl(extract.dark(temp, hours = c(18,22)))
    chl2 = get.Chl(temp)

    
    plot(chl.dark$Chl, fo.dark$Fo/ sigma.dark$Sigma, pch = 20, ylim = c(0, 2), xlim = c(0, 6), ylab = 'Fo/Sigma',
         main = paste('Dark - Cycle', i), xaxs='i', yaxs='i', xlab = 'Measured Chl')

    plot(chl.day$Chl, fo.day$Fo/ sigma.day$Sigma, pch = 20, ylim = c(0, 2), xlim = c(0, 6), ylab = 'Fo/Sigma',
         main = paste('Day - Cycle', i), xaxs='i', yaxs='i', xlab = 'Measured Chl')

    plot(chl2$Chl, fo$Fo/ sigma$Sigma, pch = 20, ylim = c(0, 2), xlim = c(0, 6), ylab = 'Fo/Sigma',
         main = paste('All - Cycle', i), xaxs='i', yaxs='i', xlab = 'Measured Chl')
}

par(mfrow=c(2,3))

plot.boxplot(4, ylim = c(0,2), ylab = 'Fo/Sigma', xlab = 'Cycle', main = 'ALL')

for (i in 1:4) {
    if (i == 1) { temp = extract.cycle(data, P1706.1()) }
    if (i == 2) { temp = extract.cycle(data, P1706.2()) }
    if (i == 3) { temp = extract.cycle(data, P1706.3()) }
    if (i == 4) { temp = extract.cycle(data, P1706.4()) }
    
    fo = get.Fo(temp)
    sigma = get.Sigma(temp)
    add.boxplot.box(i, fo$Fo/sigma$Sigma, col = 'grey')
}

plot.boxplot(4, ylim = c(0,2), ylab = 'Fo/Sigma', xlab = 'Cycle', main = 'Day')

for (i in 1:4) {
    if (i == 1) { temp = extract.cycle(data, P1706.1()) }
    if (i == 2) { temp = extract.cycle(data, P1706.2()) }
    if (i == 3) { temp = extract.cycle(data, P1706.3()) }
    if (i == 4) { temp = extract.cycle(data, P1706.4()) }
    
    fo.day = get.Fo(extract.dark(temp, hours = c(18,22)))
    sigma.day = get.Sigma(extract.dark(temp, hours = c(18,22)))
    add.boxplot.box(i, fo.day$Fo/sigma.day$Sigma, col = 'grey')
}

plot.boxplot(4, ylim = c(0,2), ylab = 'Fo/Sigma', xlab = 'Cycle', main = 'Night')

for (i in 1:4) {
    if (i == 1) { temp = extract.cycle(data, P1706.1()) }
    if (i == 2) { temp = extract.cycle(data, P1706.2()) }
    if (i == 3) { temp = extract.cycle(data, P1706.3()) }
    if (i == 4) { temp = extract.cycle(data, P1706.4()) }
    
    fo.dark = get.Fo(extract.dark(temp, hours = c(7,13)))
    sigma.dark = get.Sigma(extract.dark(temp, hours = c(7,13)))
    add.boxplot.box(i, fo.dark$Fo/sigma.dark$Sigma, col = 'grey')
}


################################################################################################
plot.boxplot(4, ylim = c(0,6), ylab = 'Measured Chl', xlab = 'Cycle', main = 'ALL')

for (i in 1:4) {
    if (i == 1) { temp = extract.cycle(data, P1706.1()) }
    if (i == 2) { temp = extract.cycle(data, P1706.2()) }
    if (i == 3) { temp = extract.cycle(data, P1706.3()) }
    if (i == 4) { temp = extract.cycle(data, P1706.4()) }
    
    chl = get.Chl(temp)
    add.boxplot.box(i, chl$Chl, col = 'grey')
}

plot.boxplot(4, ylim = c(0,6), ylab = 'Measured Chl', xlab = 'Cycle', main = 'Day')

for (i in 1:4) {
    if (i == 1) { temp = extract.cycle(data, P1706.1()) }
    if (i == 2) { temp = extract.cycle(data, P1706.2()) }
    if (i == 3) { temp = extract.cycle(data, P1706.3()) }
    if (i == 4) { temp = extract.cycle(data, P1706.4()) }

    chl.day = get.Chl(extract.dark(temp, hours = c(18,22)))
    add.boxplot.box(i, chl.day$Chl, col = 'grey')
}

plot.boxplot(4, ylim = c(0,6), ylab = 'Measured Chl', xlab = 'Cycle', main = 'Night')

for (i in 1:4) {
    if (i == 1) { temp = extract.cycle(data, P1706.1()) }
    if (i == 2) { temp = extract.cycle(data, P1706.2()) }
    if (i == 3) { temp = extract.cycle(data, P1706.3()) }
    if (i == 4) { temp = extract.cycle(data, P1706.4()) }
    
    chl.dark = get.Chl(extract.dark(temp, hours = c(7,13)))
    add.boxplot.box(i, chl.dark$Chl, col = 'grey')
}

fit.webb = calc.frrf.fit(data[[20]], Webb1974)
fit.platt = calc.frrf.fit(data[[20]], Platt1980)
fit.platt2 = calc.frrf.fit(data[[20]], Platt1980, type = 2)
fit.platt3 = calc.frrf.fit(data[[20]], Platt1980, type = 3)
fit.Jass = calc.frrf.fit(data[[20]], Jassby1976)
fit.Eil = calc.frrf.fit(data[[20]], Eilers1988)

plot.fit = function(fit) {
    plot(fit$obs$E, fit$obs$Obs, ylab = 'JVPII', xlab = 'E')
    lines(fit$ideal$E, fit$ideal$pred, lwd = 2)
    mtext(paste('Fit:', fit$model.name , '- Type:', fit$type), line = 0.25, adj = 1, cex = 0.6)
}

par(mfrow=c(2,2))
plot.fit(fit.webb)
plot.fit(fit.platt)
plot.fit(fit.platt2)
plot.fit(fit.Eil)

#par(mfrow=c(2,2))
plot(data[[20]]$A$E, data[[20]]$A$JVPII, ylab = 'JVPII', xlab = 'E')

x = c(1:1e5)
lines(x, Webb1974(fit.webb$p$min$a[1], fit.webb$p$min$b[1], x), col = 'black', lwd = 2)
lines(x, Platt1980(fit.platt$p$min$a[1], fit.platt$p$min$b[1], fit.platt$p$min$c[1], x), col = 'blue', lwd = 2)
lines(x, Jassby1976(fit.Jass$p$min$a[1], fit.Jass$p$min$b[1], E = x), col = 'red', lwd = 2)
lines(x, Eilers1988(fit.Eil$p$min$a[1], fit.Eil$p$min$b[1], fit.Eil$p$min$c[1], x), col = 'darkgreen', lwd = 2)

legend(500, 2, legend = c('Webb1974', 'Platt1980', 'Jassby1976', 'Eilers1988'),
       text.col = c('black', 'blue', 'red', 'darkgreen'))

load('../CCE-shipdata/MIMS-TBK/RStates/CTD.all.rdata')
ctd$Hour = get.hour(conv.time.excel(ctd$DateTime))

ctd = ctd[ctd$Cast != 71 & ctd$Cast != 70,] ## bad casts
ctd.oxy = ctd
ctd = ctd[ctd$PAR.surface > 10  & !is.na(ctd$Par.percent) & !is.na(ctd$Cast),]
ctd$PAR[ctd$PAR > ctd$PAR.surface] = ctd$PAR.surface[ctd$PAR > ctd$PAR.surface]
ctd$Par.percent[ctd$Par.percent> 1] = 1

for (i in unique(ctd$Cast)) {
    l = which(ctd$Cast == i)
    l.max = which.max(ctd$Depth[l])
    l = l[l.max]
    
    temp = ctd[l,]
    temp$Par.percent = 0
    
    for (k in 1:50) {
        temp$Depth = temp$Depth + 1
        temp$Pressure = temp$Pressure + 1
        ctd = rbind(ctd, temp)
    }
}

get.par.field = function(cycle, xlim) {
    
    l = which(ctd$Cycle == cycle & ctd$PAR.surface > 600)

    l = l[!is.na(ctd$DateTime[l]) & !is.na(ctd$Par.percent[l])]
    xlim = conv.time.excel(xlim, rev = TRUE)
    trans1 = build.section(x = ctd$DateTime[l], y = ctd$Pressure[l], z = ctd$Par.percent[l], xlim = xlim, neighborhood = -1,
                           field.name = 'PAR.percent', ylim = c(0,100), x.scale = 1/24, y.scale = 1/2, x.factor = 8,
                           uncertainty = 0.1)

    trans1$grid$PAR = NA
    for (i in trans1$x) {
        k = which(trans1$grid$x == i)
        l = which.closest.time(conv.time.excel(i), ship.data$DT)
        trans1$grid$PAR[k] = trans1$grid$PAR.percent[k] * ship.data$PA[l]
    }
    
    trans1$grid$x = conv.time.excel(trans1$grid$x, tz = 'GMT')
    trans1$x = conv.time.excel(trans1$x)
    trans1$data$x = conv.time.excel(trans1$data$x)

    trans1
}


trans1 = get.par.field('Cycle1', P1706.1())
trans2 = get.par.field('Cycle2', P1706.2())
trans3 = get.par.field('Cycle3', P1706.3())
trans4 = get.par.field('Cycle4', P1706.4())

plot.section.summary = function(trans, text = 'Cycle 1', ylim = c(100,0)) {
    plot.section(trans, ylim = ylim, mark.points = FALSE, log = TRUE, zlim = c(10, 2000), 
             pal = 'ocean.solar', field = 'PAR', ylab = 'Depth (m)')
    add.contour(trans, field = 'PAR.percent', levels = 0.01, 1e-3, col = 'white', lwd = 2, f = 0.1)
    
    lines(mld$time, mld$mld, lty = 3, col = 'darkgrey', lwd = 2)
    lines(ship.data$DT, -ship.data$PA/150 + 100, col = '#ff3030')
    
    mtext(text, line = 0.25, cex = 0.75, adj = 0)
}

par(plt = c(0.2, 0.8, 0.1, 0.9))
par(mfrow=c(2,2))

plot.section.summary(trans1, text = 'Cycle 1')

par(plt = c(0.1, 0.7, 0.1, 0.9))
plot.section.summary(trans2, text = 'Cycle 2')

par(plt = c(0.2, 0.8, 0.1, 0.9))
plot.section.summary(trans3, text = 'Cycle 3')

par(plt = c(0.1, 0.7, 0.1, 0.9))
plot.section.summary(trans4, text = 'Cycle 4')


par(mfrow=c(1,1))
add.colorbar(10, 2000, pal = 'ocean.solar', log = TRUE, base = 10, units = 'PAR',
             labels = c(10,30,100,300,1000,2000), x.pos = 0.9, width = 0.02)

add.gpp = function(trans) {
    trans$grid$GPP = NA
    trans$grid$psi.ec = NA
    trans$grid$JVPII = NA
    trans$grid$NSV = NA
    trans$grid$Chl = NA
    trans$grid$Dist = NA
    
    frrf.times = get.frrf.times(data)
    
    mixed.layer = approx(mld$time, mld$mld, trans$x)$y
    
    for (j in 1:length(trans$x)) { ## each time point
        t = trans$x[j]
        k = which(trans$grid$x == t) ## which grid points match this time
        i = which.closest.time(conv.time.unix(t), frrf.times)
        
        ## FRRF Stuff
        f = data[[i]]$Params$f ## Adjustment factor based on nPSII shift
        
        ## E:C ratio
        psi = calc.frrf.psi.ec(data[[i]])
        l = which(!is.na(psi))
        psi = approx(x = data[[i]]$A$E[l], y = psi[l], xout = 0:1e4, rule = 2)$y # interpolate
        nsv = approx(x = data[[i]]$A$E[l], y = data[[i]]$A$NSV[l], xout = 0:1e4, rule = 2)$y # interpolate
        
        jvpii = rep(NA, length(k))
        
        for (ii in 1:length(k)) {
            light = which.min((data[[i]]$Fit$ideal$E - trans$grid$PAR[k[ii]])^2)
            jvpii[ii] = data[[i]]$Fit$ideal$pred[light]
            jvpii[ii] = jvpii[ii] 
        }
        
        trans$grid$JVPII[k] = jvpii
        trans$grid$psi.ec[k] = psi[round(trans$grid$PAR[k])+1]
        trans$grid$NSV[k] = nsv[round(trans$grid$PAR[k])+1]
        trans$grid$GPP[k] = jvpii / psi[round(trans$grid$PAR[k])+1] * 1000 * data[[i]]$Params$f # mmol C m-3 d-1
        trans$grid$Chl[k] = data[[i]]$Params$Chl
        trans$grid$Dist[k] = data[[i]]$Params$Dist
        
        #kk = which(trans$grid$y[k] > mixed.layer[j])
        
        #trans$grid$JVPII[k[kk]] = NA
        #trans$grid$psi.ec[k[kk]] = NA
        #trans$grid$NSV[k[kk]] = NA
        #trans$grid$GPP[k[kk]] = NA
        #trans$grid$Chl[k[kk]] = NA
        #trans$grid$Dist[k[kk]] = NA
    }
    trans$grid$GPP = abs(trans$grid$GPP)
    trans
}

options(warn = -1)
trans1 = add.gpp(trans1)
trans2 = add.gpp(trans2)
trans3 = add.gpp(trans3)
trans4 = add.gpp(trans4)

#pdf('../Images/FRRf Diagnostic (cycle4).pdf')

par(mfrow=c(2,2), plt = c(0.1,0.7,0.1,0.9))

temp = trans1

plot.section(temp, field = 'PAR', ylim = c(35,0), pal = 'ocean.solar', mark.points = FALSE, zlim = c(0,3000))
polygon(c(mld$time, rev(mld$time)), c(mld$mld, rev(mld$mld+100)), col = 'white')
add.colorbar(0, 3000, labels = c(0:5)*1000, pal = 'ocean.solar', x.pos = 0.75, units = 'PAR', height = 0.8, cex.units = 0.8)

plot.section(temp, field = 'psi.ec', ylim = c(35,0), pal = 'ocean.thermal', mark.points = FALSE, zlim = c(5,10))
polygon(c(mld$time, rev(mld$time)), c(mld$mld, rev(mld$mld+100)), col = 'white')
add.colorbar(5, 10, labels = c(5:10), pal = 'ocean.thermal', x.pos = 0.75, units = 'e:C', height = 0.8, cex.units = 0.8)

plot.section(temp, field = 'JVPII', ylim = c(35,0), pal = 'ocean.turbid', mark.points = FALSE, zlim = c(0,3))
polygon(c(mld$time, rev(mld$time)), c(mld$mld, rev(mld$mld+100)), col = 'white')
add.colorbar(0, 3, labels = c(0:3), pal = 'ocean.turbid', x.pos = 0.75, units = 'JVPII', height = 0.8, cex.units = 0.8)

plot.section(temp, field = 'NSV', ylim = c(35,0), pal = 'ocean.haline', mark.points = FALSE, zlim = c(0,4))
polygon(c(mld$time, rev(mld$time)), c(mld$mld, rev(mld$mld+100)), col = 'white')
add.colorbar(0, 4, labels = c(0:5), pal = 'ocean.haline', x.pos = 0.75, units = 'NSV', height = 0.8, cex.units = 0.8)

plot.section(temp, field = 'GPP', ylim = c(35,0), pal = 'ocean.haline', mark.points = FALSE, zlim = c(0,40))
polygon(c(mld$time, rev(mld$time)), c(mld$mld, rev(mld$mld+100)), col = 'white')
add.colorbar(0, 40, labels = c(0:10)*5, pal = 'ocean.haline', x.pos = 0.75, units = 'GPP', height = 0.8, cex.units = 0.8)

plot.section(temp, field = 'Dist', ylim = c(35,0), pal = 'ocean.deep', mark.points = FALSE, zlim = c(0,12))
polygon(c(mld$time, rev(mld$time)), c(mld$mld, rev(mld$mld+100)), col = 'white')
add.colorbar(0, 12, labels = c(0:8)*2, pal = 'ocean.deep', x.pos = 0.75, units = 'Dist', height = 0.8, cex.units = 0.8)

plot.section(temp, field = 'Chl', ylim = c(35,0), mark.points = FALSE, zlim = c(0,8), pal = 'ocean.algae')
polygon(c(mld$time, rev(mld$time)), c(mld$mld, rev(mld$mld+100)), col = 'white')
add.colorbar(0, 8, labels = c(0:5)*2, x.pos = 0.75, units = '[Chl]', height = 0.8, cex.units = 0.8)

#dev.off()

par(mfrow=c(2,2), plt = c(0.1,0.7,0.1,0.9))
plot.gpp = function(trans, xlim) {
    plot.section(trans, field = 'GPP', ylim = c(35,0), pal = 'inferno', log = TRUE, zlim = c(2, 1000), mark.points = FALSE,
                 col.high = 'red', col.low = '', xlim = xlim)
    polygon(c(mld$time, rev(mld$time)), c(mld$mld, rev(mld$mld+100)), col = 'white')
}

plot.gpp(trans1, P1706.1())
plot.gpp(trans2, P1706.2())
plot.gpp(trans3, P1706.3())
plot.gpp(trans4, P1706.4())

add.colorbar(10, 1000, log = TRUE, labels = c(1, 10, 100, 1000), x.pos = 0.8, pal = 'inferno',
             units = '[mmol C m-3 d-1]', cex.units = 0.6, height = 0.85)

background$GPP.ml = NA

for (i in 1:nrow(background)) {
    
    if (background$Time[i] >= P1706.1()[1] & background$Time[i] <= P1706.1()[2]) {
        ## Find time
        t = trans1$x[which.closest.time(background$Time[i], trans1$x)]
        
        ## Find grid entries
        k = which(trans1$grid$x == t & trans1$grid$y <= background$MLD[i])
        k2 = which(trans1$grid$x == t)
        
        background$GPP.ml[i] = mean(trans1$grid$GPP[k], na.rm = TRUE) * background$MLD[i]
    }
    
    if (background$Time[i] >= P1706.2()[1] & background$Time[i] <= P1706.2()[2]) {
        ## Find time
        t = trans2$x[which.closest.time(background$Time[i], trans2$x)]
        
        ## Find grid entries
        k = which(trans2$grid$x == t & trans2$grid$y <= background$MLD[i])
        k2 = which(trans2$grid$x == t)
        
        background$GPP.ml[i] = mean(trans2$grid$GPP[k], na.rm = TRUE) * background$MLD[i]
    }
    
    if (background$Time[i] >= P1706.3()[1] & background$Time[i] <= P1706.3()[2]) {
        ## Find time
        t = trans3$x[which.closest.time(background$Time[i], trans3$x)]
        
        ## Find grid entries
        k = which(trans3$grid$x == t & trans3$grid$y <= background$MLD[i])
        k2 = which(trans3$grid$x == t)
        
        background$GPP.ml[i] = mean(trans3$grid$GPP[k], na.rm = TRUE) * background$MLD[i]
    }
    
    if (background$Time[i] >= P1706.4()[1] & background$Time[i] <= P1706.4()[2]) {
        ## Find time
        t = trans4$x[which.closest.time(background$Time[i], trans4$x)]
        
        ## Find grid entries
        k = which(trans4$grid$x == t & trans4$grid$y <= background$MLD[i])
        k2 = which(trans4$grid$x == t)
        
        background$GPP.ml[i] = mean(trans4$grid$GPP[k], na.rm = TRUE) * background$MLD[i]
    }
}

summary(background$GPP.ml)

add.bar = function(x, y, width = 1, col = 'grey', border = TRUE, baseline = 0) {
    for (i in 1:length(x)) {
        rect(x[i] - width/2, baseline, x[i] + width/2, y[i], col = col, border = border)
    }
}

head(npp[npp$Cycle == 'Cycle3',])

par(mfrow=c(2,2))

plot(background$Time, background$GPP.ml, xlab = '', ylab = 'GPP (mmol C m-2 d-1)', yaxs = 'i',
     type = 'l', lwd = 2, ylim = c(0, 5e3), xlim = P1706.1(), xaxs = 'i')
axis(1, at = c(1e4:1e5)*86400, labels = conv.time.unix(c(1e4:1e5)*86400))
mtext('Cycle 1 - ML GPP', line = 0.25, adj = 0, cex = 0.75)
add.night(ship.data$DT, ship.data$PA)
add.bar(npp$Time, npp$NPP, width = 85000, col = '#00000030')

plot(background$Time, background$GPP.ml, xlab = '', ylab = 'GPP (mmol C m-2 d-1)', yaxs = 'i', type = 'l',
     lwd = 2, ylim = c(0, 5e3), xlim = P1706.2(), xaxs = 'i')
axis(1, at = c(1e4:1e5)*86400, labels = conv.time.unix(c(1e4:1e5)*86400))
mtext('Cycle 2 - ML GPP', line = 0.25, adj = 0, cex = 0.75)
add.bar(npp$Time, npp$NPP, width = 85000, col = '#00000030')
add.night(ship.data$DT, ship.data$PA)

plot(background$Time, background$GPP.ml, xlab = '', ylab = 'GPP (mmol C m-2 d-1)', yaxs = 'i',
     type = 'l', lwd = 2, ylim = c(0, 500), xlim = P1706.3(), xaxs = 'i')
axis(1, at = c(1e4:1e5)*86400, labels = conv.time.unix(c(1e4:1e5)*86400))
mtext('Cycle 3 - ML GPP', line = 0.25, adj = 0, cex = 0.75)
add.bar(npp$Time, npp$NPP, width = 86200, col = '#00000030')
add.night(ship.data$DT, ship.data$PA)

plot(background$Time, background$GPP.ml, xlab = '', ylab = 'GPP (mmol C m-2 d-1)', yaxs = 'i',
     type = 'l', lwd = 2, ylim = c(0, 500), xlim = P1706.4(), xaxs = 'i')
axis(1, at = c(1e4:1e5)*86400, labels = conv.time.unix(c(1e4:1e5)*86400))
mtext('Cycle 4 - ML GPP', line = 0.25, adj = 0, cex = 0.75)
add.bar(npp$Time, npp$NPP, width = 86000, col = '#00000030')
add.night(ship.data$DT, ship.data$PA)

l = which(background$Time >= P1706.1()[1] & background$Time <= P1706.1()[2])
write.xlsx(background[l,], file = '../Images/Background - Cycle1 new.xlsx')

l = which(background$Time >= P1706.2()[1] & background$Time <= P1706.2()[2])
write.xlsx(background[l,], file = '../Images/Background - Cycle2 new.xlsx')

l = which(background$Time >= P1706.3()[1] & background$Time <= P1706.3()[2])
write.xlsx(background[l,], file = '../Images/Background - Cycle3 new.xlsx')

l = which(background$Time >= P1706.4()[1] & background$Time <= P1706.4()[2])
write.xlsx(background[l,], file = '../Images/Background - Cycle4 new.xlsx')

frrf = data.frame(Time = get.frrf.times(data),
                  GPP = NA, PAR = NA, MLD = NA, FvFm = NA, sigma = NA,
                  Tau = NA, alpha = NA, beta = NA, Ps = NA, Pm = NA, FqFv = NA)

frrf$GPP = approx(background$Time, background$GPP.ml, frrf$Time, rule = 2)$y
frrf$MLD = approx(background$Time, background$MLD, frrf$Time, rule = 2)$y
frrf$PAR = approx(background$Time, background$PAR, frrf$Time, rule = 2)$y

for (i in 1:length(data)) {
    frrf$Tau[i] = data[[i]]$A$TauES[1]
    frrf$alpha[i] = data[[i]]$Fit$p$min$a
    frrf$beta[i] = data[[i]]$Fit$p$min$b
    frrf$Ps[i] = data[[i]]$Fit$p$min$c
    frrf$Pm[i] = max(data[[i]]$Fit$ideal$pred)
    frrf$FvFm[i] = data[[i]]$A$Fv.Fm[1]
    frrf$sigma[i] = data[[i]]$A$Sigma[1]
    frrf$FqFv[i] = data[[i]]$A$FqFv[1]
}



pdf('../Images/FRRF_Summary of Cycle 34z.pdf')

xlim = P1706.3()
xlim2 = P1706.4()
par(mfrow = c(7,2), plt = c(0.2,0.95, 0, 1))

## A
plot(ship.data$DT, ship.data$PA, xlim = xlim, type = 'l', ylim = c(0,3e3), yaxs='i', ylab = 'PAR', xaxs = 'i')
polygon(x = c(ship.data$DT, rev(ship.data$DT)), y = c(ship.data$PA, rep(0,length(ship.data$PA))), col = '#00000020')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

plot(ship.data$DT, ship.data$PA, xlim = xlim2, type = 'l', ylim = c(0,3e3), yaxs='i', ylab = 'PAR', xaxs = 'i')
polygon(x = c(ship.data$DT, rev(ship.data$DT)), y = c(ship.data$PA, rep(0,length(ship.data$PA))), col = '#00000020')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

## B
plot(frrf$Time, frrf$GPP, xlim = xlim, ylim = c(0, 1e3), yaxs = 'i', ylab = 'GPP (mmol C m-2 d-1)', pch = 16, xaxs = 'i')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

plot(frrf$Time, frrf$GPP, xlim = xlim2, ylim = c(0, 1e3), yaxs = 'i', ylab = 'GPP (mmol C m-2 d-1)', pch = 16, xaxs = 'i')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

## C
plot(frrf$Time, frrf$Pm, xlim = xlim, yaxs = 'i', ylab = 'Pm', pch = 16, ylim = c(0,6), xaxs = 'i')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

plot(frrf$Time, frrf$Pm, xlim = xlim2, yaxs = 'i', ylab = 'Pm', pch = 16, ylim = c(0,6), xaxs = 'i')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

##D
plot(frrf$Time, frrf$alpha, xlim = xlim, yaxs = 'i', ylab = 'Alpha', pch = 16, ylim = c(0, 0.01),
     yaxs = 'i', xaxs = 'i')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

plot(frrf$Time, frrf$alpha, xlim = xlim2, yaxs = 'i', ylab = 'Alpha', pch = 16, ylim = c(0, 0.01),
     yaxs = 'i', xaxs = 'i')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

## E
plot(frrf$Time, frrf$FvFm, xlim = xlim, yaxs = 'i', ylab = 'Fv/Fm', pch = 16, yaxs = 'i',
     ylim = c(0.2, 0.6), xaxs = 'i')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

plot(frrf$Time, frrf$FvFm, xlim = xlim2, yaxs = 'i', ylab = 'Fv/Fm', pch = 16, yaxs = 'i',
     ylim = c(0.2, 0.6), xaxs = 'i')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

## F
plot(frrf$Time, frrf$sigma, xlim = xlim, yaxs = 'i', ylab = 'Sigma (nm2)', pch = 16, ylim = c(3,10),
     yaxs = 'i', xaxs = 'i')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

plot(frrf$Time, frrf$sigma, xlim = xlim2, yaxs = 'i', ylab = 'Sigma (nm2)', pch = 16, ylim = c(3,10),
     yaxs = 'i', xaxs = 'i')
#axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

## G
plot(frrf$Time, 1/frrf$Tau, xlim = xlim, yaxs = 'i', ylab = '1/Tau', pch = 16,
     yaxs = 'i', ylim = c(4e-4, 10e-4), xaxs = 'i')
axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

plot(frrf$Time, 1/frrf$Tau, xlim = xlim2, yaxs = 'i', ylab = '1/Tau', pch = 16,
     yaxs = 'i', ylim = c(4e-4, 10e-4), xaxs = 'i')
axis.POSIXct(1, P1706.1())
add.night(ship.data$DT, ship.data$PA)

dev.off()

write.xlsx(frrf, file = '../Images/FRRF data.xlsx')

temp = data.frame(Time = rep(NA, length(data)), NSV = NA)

plot(NULL,NULL, xlim = c(P1706.1()[1], P1706.4()[2]), ylim = c(0,2.5), ylab = 'NSV @ 500', xlab = 'Time')

add.shade(P1706.1())
add.shade(P1706.2())
add.shade(P1706.3())
add.shade(P1706.4())

for (i in 1:length(data)) {
    points(data[[i]]$Datetime, approx(data[[i]]$A$E, data[[i]]$A$NSV, 500)$y, pch = 20)
    temp$Time[i] = data[[i]]$Datetime
    temp$NSV[i] = approx(data[[i]]$A$E, data[[i]]$A$NSV, 500)$y
}
temp$Time = conv.time.unix(temp$Time)


write.xlsx(temp, '../Images/NSV Summary.xlsx')

load('../Documents/CCE/NCP1706/RStates/ship.ncp.15.rdata')

par(mfrow=c(2,2))
plot(ship.data$DT, ship.data$k, pch = 20, xlim = P1706.1())
add.night(ship.data$DT, ship.data$PA)

plot(ship.data$DT, ship.data$TW, pch = 20, xlim = P1706.1())
add.night(ship.data$DT, ship.data$PA)


plot(ship.data$DT, ship.data$k, pch = 20, xlim = P1706.2())
add.night(ship.data$DT, ship.data$PA)

plot(ship.data$DT, ship.data$TW, pch = 20, xlim = P1706.2())
add.night(ship.data$DT, ship.data$PA)

plot(ship.data$DT, ship.data$k, pch = 20, xlim = P1706.3())
add.night(ship.data$DT, ship.data$PA)

plot(ship.data$DT, ship.data$TW, pch = 20, xlim = P1706.3())
add.night(ship.data$DT, ship.data$PA)

plot(ship.data$DT, ship.data$k, pch = 20, xlim = P1706.4())
add.night(ship.data$DT, ship.data$PA)

plot(ship.data$DT, ship.data$TW, pch = 20, xlim = P1706.4())
add.night(ship.data$DT, ship.data$PA)

plot(ship.data$DT, ship.data$MLD)
plot(ship.data$DT, ship.data$TW)
add.night(ship.data$DT, ship.data$PA)

summary(background$NCP.RT.eims)


