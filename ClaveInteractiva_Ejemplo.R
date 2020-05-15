################################################################################
################################################################################
#
# INTRODUCCIÓN
# El código en este guión ilustra el concepto de una clave interactiva y
# probabilística. Una clave es interactiva (sensu, Dallwitz et al. 2000.
# Principles of interactive keys. http://biodiversity.uno.edu/delta, 3) cuando 
# cualquier combinación de caracteres relevantes puede utilizarse para tratar
# de determinar la especie a la que pertenece un espécimen dado. Una clave es
# probabilística cuando estima la probabilidad de que un espécimen dado
# pertenezca a cada una de las especies consideradas en la clave. La clave
# en este guión se basa en el modelo probabilístico de mezclas normales
# (McLachlan y Peel. 2000. Finite Mixture Models, Willey Series in Probability
# and Statistics, John Wiley & Sons, New York) aplicado al contexto de
# delimitación y descripción de especies (Zapata y Jiménez. 2012. Species
# Delimitation: Inferring Gaps in Morphology across Geography. Syst. Biol. 61:
# 179–194; Cadena et al. 2018. Issues and Perspectives in Species Delimitation
# using Phenotypic Data: Atlantean Evolution in Darwin’s Finches. Syst. Biol.
# 67: 181–194).
#
# ARCHIVOS NECESARIOS PARA UTILIZAR EL CÓDIGO
# Ninguno.
# 
# CONTENIDO
# 1. Preparar paquetes requeridos.
# 2. Simulación de la distribución fenotípica de tres especies
# 3. Representación gráfica de la distribución fenotípica de las tres especies.
# 4. Simulación del fenotipo de un espécimen.
# 5. Representación gráfica del fenotipo del espácimen en el contexto de la
#    distribución fenotípica de las tres especies.
# 6. Simulación del fenotipo medido del espécimen.
# 7. Representación gráfica del fenotipo medido en el espécimen en el contexto
#    de la distribución fenotípica de las tres especies.
# 8. Verosimilitud del fenotipo del espécimen respecto a las distribuciones
#    fenotípicas de cada una de las tres especies.
# 9. Verosimilitud del fenotipo medido en el espécimen respecto a las
#    distribuciones fenotípicas de cada una de las tres especies.
# 10. Comparación las verosimilitudes del fenotipo y del fenotipo medido en el
#     espécimen.
# 11. Probabilidad de muestrear un fenotipo al menos tan extremo como el
#     observado.
# 12. Probabilidad de muestrear un fenotipo medido al menos tan extremo como
#     el observado.
# 13. Comparación las probabilidades de muestrear un fenotipo y un fenotipo
#     medido al menos tan extremos como los observados.
#
################################################################################
################################################################################


################################################################################
# 1. Preparar paquetes requeridos.
################################################################################

#cargue paquetes (asegurese de instalarlos antes de cargarlos)
library(ellipse)
library(mvtnorm)
library(clusterGeneration)


################################################################################
# 2. Simulación de la distribución fenotípica de tres especies.
################################################################################

#defina el número de caracteres fenotípicos en la descripción de las especies
p <- 13

#simule aleatoriamente el vector de medias para cada especie
M.1 <- runif(p, -1, 1)
M.2 <- runif(p, -1, 1)
M.3 <- runif(p, -1, 1)
M <- list(M.1, M.2, M.3)
M

#simule aleatoriamente la matriz de varianza-covarianza para cada especie
VCV.1 <- genPositiveDefMat(covMethod="unifcorrmat", rangeVar=c(1e-10,1), dim=p)$Sigma
VCV.2 <- genPositiveDefMat(covMethod="unifcorrmat", rangeVar=c(1e-10,1), dim=p)$Sigma
VCV.3 <- genPositiveDefMat(covMethod="unifcorrmat", rangeVar=c(1e-10,1), dim=p)$Sigma
VCV <- list(VCV.1, VCV.2, VCV.3)
VCV

################################################################################
# 3. Representación gráfica de la distribución fenotípica de las tres especies.
################################################################################

#selecione los colores para representar cada una de tres especies
color.especies <- c("red", "blue", "green")


################################################################################
# 3.1. Distribución bivariada en dos caracteres particulares.
################################################################################

#escoja los dos caracteres fenotípicos de interés
CF1 <- 13
CF2 <- 7

plot(ellipse(VCV[[1]][c(CF1, CF2), c(CF1, CF2)], centre=M[[1]][c(CF1, CF2)]), type="l", col=color.especies[1],
	xlab=paste("Carácter fenotípico", CF1), ylab=paste("Carácter fenotípico", CF2), xlim=c(-3,3), ylim=c(-3,3))
points(M[[1]][CF1], M[[1]][CF2], col=color.especies[1], pch=19, cex=1.5)
for(i in 2:3){
	points(ellipse(VCV[[i]][c(CF1, CF2), c(CF1, CF2)], centre=M[[i]][c(CF1, CF2)]), type="l", col=color.especies[i])
	points(M[[i]][CF1], M[[i]][CF2], col=color.especies[i], pch=19, cex=1.5)
}

################################################################################
# 3.2. Distribución bivariada en todos los pares de caracteres.
################################################################################

indice.par.caracteres <- combn(1:p,2)
	for(j in 1:ncol(indice.par.caracteres)){
		CF1 <- indice.par.caracteres[1,j]
		CF2 <- indice.par.caracteres[2,j]
		plot(ellipse(VCV[[1]][c(CF1, CF2), c(CF1, CF2)], centre=M[[1]][c(CF1, CF2)]), type="l", col=color.especies[1],
			xlab=paste("Carácter fenotípico", CF1), ylab=paste("Carácter fenotípico", CF2), xlim=c(-3,3), ylim=c(-3,3))
		points(M[[1]][CF1], M[[1]][CF2], col=color.especies[1], pch=19, cex=1.5)
		for(i in 2:3){
			points(ellipse(VCV[[i]][c(CF1, CF2), c(CF1, CF2)], centre=M[[i]][c(CF1, CF2)]), type="l", col=color.especies[i])
			points(M[[i]][CF1], M[[i]][CF2], col=color.especies[i], pch=19, cex=1.5)
		}	
		Sys.sleep(2)
	}


################################################################################
# 4. Simulación del fenotípo de un espécimen.
################################################################################

fenotipo.especimen <- runif(p, -3, 3)
fenotipo.especimen


################################################################################
# 5. Representación gráfica del fenotipo del espécimen en el contexto de la
#    distribución fenotípica de las tres especies.
################################################################################

indice.par.caracteres <- combn(1:p,2)
	for(j in 1:ncol(indice.par.caracteres)){
		CF1 <- indice.par.caracteres[1,j]
		CF2 <- indice.par.caracteres[2,j]
		plot(ellipse(VCV[[1]][c(CF1, CF2), c(CF1, CF2)], centre=M[[1]][c(CF1, CF2)]), type="l", col=color.especies[1],
			xlab=paste("Carácter fenotípico", CF1), ylab=paste("Carácter fenotípico", CF2), xlim=c(-3,3), ylim=c(-3,3))
		points(M[[1]][CF1], M[[1]][CF2], col=color.especies[1], pch=19, cex=1.5)
		for(i in 2:3){
			points(ellipse(VCV[[i]][c(CF1, CF2), c(CF1, CF2)], centre=M[[i]][c(CF1, CF2)]), type="l", col=color.especies[i])
			points(M[[i]][CF1], M[[i]][CF2], col=color.especies[i], pch=19, cex=1.5)
		}	
		points(fenotipo.especimen[CF1], fenotipo.especimen[CF2], pch=3, cex=2)
		Sys.sleep(2)
	}


################################################################################
# 6. Simulación del fenotipo medido del espécimen.
################################################################################

#definición del número de caracteres medidos
numero.caracteres.medidos <- sample(1:p, size=1)
numero.caracteres.medidos
#numero.caracteres.medidos <- 1 #use esta línea si no quiere dejar al azar el número de carácteres medidos

#definición de los caracteres medidos
caracteres.medidos <- sort(sample(1:p, numero.caracteres.medidos))
caracteres.medidos

#definición del fenotipo medido
medidas.especimen <- fenotipo.especimen
medidas.especimen[-caracteres.medidos] <- NA
fenotipo.especimen
medidas.especimen


################################################################################
# 7. Representación gráfica del fenotipo medido en el espécimen en el contexto
#    de la distribución fenotípica de las tres especies.
################################################################################

#si el fenotipo medido sólo incluye un carácter
if(numero.caracteres.medidos<2){
	CF1 <- caracteres.medidos
	max.density <- rep(NA, times=3)
	for(i in 1:3){
		max.density[i] <- max(dnorm(seq(-3,3,0.01), mean = M[[i]][CF1], sd = sqrt(VCV[[i]][CF1,CF1])))
	}
	plot(seq(-3,3,0.01), dnorm(seq(-3,3,0.01), mean = M[[1]][CF1], sd = sqrt(VCV[[1]][CF1,CF1])), type="l", col=color.especies[1],
		ylab="Densidad de probabilidad", xlab=paste("Carácter fenotípico", CF1), ylim=c(0,max(max.density)))
	abline(v= M[[1]][CF1], col= color.especies[1])
	for(i in 2:3){
		points(seq(-3,3,0.01), dnorm(seq(-3,3,0.01), mean = M[[i]][CF1], sd = sqrt(VCV[[i]][CF1,CF1])), type="l", col=color.especies[i])
		abline(v= M[[i]][CF1], col= color.especies[i])
	}
	abline(v=medidas.especimen[CF1], lty=3, lwd=2)
}

#si el fenotipo medido incluye dos o más caracteres
if(numero.caracteres.medidos>1){
	indice.par.caracteres <- combn(caracteres.medidos,2)
	for(j in 1:ncol(indice.par.caracteres)){
		CF1 <- indice.par.caracteres[1,j]
		CF2 <- indice.par.caracteres[2,j]
		plot(ellipse(VCV[[1]][c(CF1, CF2), c(CF1, CF2)], centre=M[[1]][c(CF1, CF2)]), type="l", col=color.especies[1],
			xlab=paste("Carácter fenotípico", CF1), ylab=paste("Carácter fenotípico", CF2), xlim=c(-3,3), ylim=c(-3,3))
		points(M[[1]][CF1], M[[1]][CF2], col=color.especies[1], pch=19, cex=1.5)
		for(i in 2:3){
			points(ellipse(VCV[[i]][c(CF1, CF2), c(CF1, CF2)], centre=M[[i]][c(CF1, CF2)]), type="l", col=color.especies[i])
			points(M[[i]][CF1], M[[i]][CF2], col=color.especies[i], pch=19, cex=1.5)
		}	
		points(medidas.especimen[CF1], medidas.especimen[CF2], pch=3, cex=2)
		Sys.sleep(5)
	}
}


################################################################################
# 8. Verosimilitud del fenotipo del espécimen respecto a las distribuciones
#    fenotípicas de cada una de las tres especies.
################################################################################

#logaritmo de la verosimilitud
logv.fenotipo <- rep(NA, times=3)
logv.fenotipo[1] <- dmvnorm(fenotipo.especimen, mean= M[[1]], sigma = VCV[[1]], log = T) 
logv.fenotipo[2] <- dmvnorm(fenotipo.especimen, mean= M[[2]], sigma = VCV[[2]], log = T)
logv.fenotipo[3] <- dmvnorm(fenotipo.especimen, mean= M[[3]], sigma = VCV[[3]], log = T)
logv.fenotipo

#graficar el logaritmo de la verosimilitud
plot(logv.fenotipo, type="l", pch=19, xaxt="n", xlab="Especies", ylab="Log ( verosimilitud del fenotipo )")
points(logv.fenotipo, type="p", pch=19, cex=2, col=color.especies)
axis(1, at=1:3)

#clasificación del fenotipo medido del espécimen por máxima verosimilitud
which(logv.fenotipo==max(logv.fenotipo))
#incertidumbre de la clasificación por máxima verosimilitud
#(ver página 614, después de ecuación 7 en Fraley, C. and Raftery, A.E., 2002. Model-based clustering,
#discriminant analysis, and density estimation. Journal of the American Statistical Association, 97: 611-631)
1 - max(exp(logv.fenotipo))/sum(exp(logv.fenotipo))


################################################################################
# 9. Verosimilitud del fenotipo medido en el espécimen respecto a las
#    distribuciones fenotípicas de cada una de las tres especies.
################################################################################

#logaritmo de la verosimilitud
logv.fenotipo.medido <- rep(NA, times=3)

#si el fenotipo medido sólo incluye un carácter
if(numero.caracteres.medidos<2){
	logv.fenotipo.medido[1] <- dnorm(medidas.especimen[!is.na(medidas.especimen)], mean= M[[1]][caracteres.medidos], sd = sqrt(VCV[[1]][caracteres.medidos,caracteres.medidos]), log = T) 
	logv.fenotipo.medido[2] <- dnorm(medidas.especimen[!is.na(medidas.especimen)], mean= M[[2]][caracteres.medidos], sd = sqrt(VCV[[2]][caracteres.medidos,caracteres.medidos]), log = T) 
	logv.fenotipo.medido[3] <- dnorm(medidas.especimen[!is.na(medidas.especimen)], mean= M[[3]][caracteres.medidos], sd = sqrt(VCV[[3]][caracteres.medidos,caracteres.medidos]), log = T) 
	plot(logv.fenotipo.medido, type="l", pch=19, xaxt="n", xlab="Especies",
	ylab=paste("Log ( verosimilitud del fenotipo medido ), ", numero.caracteres.medidos, " carácter de ", "p=", p, sep=""))
	points(logv.fenotipo.medido, type="p", pch=19, cex=2, col=color.especies)
	axis(1, at=1:3)
}
logv.fenotipo.medido

#si el fenotipo medido incluye dos o más caracteres
if(numero.caracteres.medidos>1){
	logv.fenotipo.medido[1] <- dmvnorm(medidas.especimen[!is.na(medidas.especimen)], mean= M[[1]][caracteres.medidos], sigma = VCV[[1]][caracteres.medidos,caracteres.medidos], log = T) 
	logv.fenotipo.medido[2] <- dmvnorm(medidas.especimen[!is.na(medidas.especimen)], mean= M[[2]][caracteres.medidos], sigma = VCV[[2]][caracteres.medidos,caracteres.medidos], log = T)
	logv.fenotipo.medido[3] <- dmvnorm(medidas.especimen[!is.na(medidas.especimen)], mean= M[[3]][caracteres.medidos], sigma = VCV[[3]][caracteres.medidos,caracteres.medidos], log = T)
	plot(logv.fenotipo.medido, type="l", pch=19, xaxt="n", xlab="Especies",
		ylab=paste("Log ( verosimilitud del fenotipo medido ), ", numero.caracteres.medidos, " caracteres de ", "p=", p, sep=""))
	points(logv.fenotipo.medido, type="p", pch=19, cex=2, col=color.especies)
	axis(1, at=1:3)
}
logv.fenotipo.medido


################################################################################
# 10. Comparación las verosimilitudes del fenotipo y del fenotipo medido en el
#     espécimen.
################################################################################ 

plot(logv.fenotipo, logv.fenotipo.medido, pch=19, cex=2, col=color.especies,
	xlab=paste("Log ( verosimilitud del fenotipo ), p=", p, sep=""), 
	ylab=paste("Log ( verosimilitud del fenotipo medido ), ", numero.caracteres.medidos, " caracteres", sep=""))
#legend("bottomright", c("Especie 1", "Especie 2", "Especie 3"), pch=19, cex=2, col=color.especies)
#legend("topright", c("Especie 1", "Especie 2", "Especie 3"), pch=19, cex=2, col=color.especies)
#legend("bottomleft", c("Especie 1", "Especie 2", "Especie 3"), pch=19, cex=2, col=color.especies)
#legend("topleft", c("Especie 1", "Especie 2", "Especie 3"), pch=19, cex=2, col=color.especies)


################################################################################
# 11. Probabilidad de muestrear un fenotipo al menos tan extremo como el observado, ver:
#     https://stats.stackexchange.com/questions/97408/relation-of-mahalanobis-distance-to-log-likelihood
#     https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Likelihood_function
################################################################################

P.fenotipo.al.menos.tan.extremo <- rep(NA, times=3)

fenotipo.C.1 <- 0.5*(log(det(VCV[[1]]))+(p*log(2*pi)))
fenotipo.DiMa.1 <- sqrt(-2*(logv.fenotipo[1]+fenotipo.C.1))
P.fenotipo.al.menos.tan.extremo[1] <- pchisq(fenotipo.DiMa.1, p, lower.tail = FALSE, log.p = FALSE)
P.fenotipo.al.menos.tan.extremo[1]

fenotipo.C.2 <- 0.5*(log(det(VCV[[2]]))+(p*log(2*pi)))
fenotipo.DiMa.2 <- sqrt(-2*(logv.fenotipo[2]+fenotipo.C.2))
P.fenotipo.al.menos.tan.extremo[2] <- pchisq(fenotipo.DiMa.2, p, lower.tail = FALSE, log.p = FALSE)
P.fenotipo.al.menos.tan.extremo[2]

fenotipo.C.3 <- 0.5*(log(det(VCV[[3]]))+(p*log(2*pi)))
fenotipo.DiMa.3 <- sqrt(-2*(logv.fenotipo[3]+fenotipo.C.3))
P.fenotipo.al.menos.tan.extremo[3] <- pchisq(fenotipo.DiMa.3, p, lower.tail = FALSE, log.p = FALSE)
P.fenotipo.al.menos.tan.extremo[3]


################################################################################
# 12. Probabilidad de muestrear un fenotipo medido al menos tan extremo como el observado, ver: 
#     https://stats.stackexchange.com/questions/97408/relation-of-mahalanobis-distance-to-log-likelihood
#     https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Likelihood_function
################################################################################

P.fenotipo.medido.al.menos.tan.extremo <- rep(NA, times=3)

#si el fenotipo medido sólo incluye un carácter
if(numero.caracteres.medidos<2){
	logv.fenotipo.medido.C.1 <- 0.5*(log(VCV[[1]][caracteres.medidos,caracteres.medidos])+(numero.caracteres.medidos*log(2*pi)))
	fenotipo.medido.DiMa.1 <- -2*(logv.fenotipo.medido[1]+logv.fenotipo.medido.C.1)
	P.fenotipo.medido.al.menos.tan.extremo[1] <- pchisq(fenotipo.medido.DiMa.1, numero.caracteres.medidos, lower.tail = FALSE, log.p = FALSE)
	P.fenotipo.medido.al.menos.tan.extremo[1]

	logv.fenotipo.medido.C.2 <- 0.5*(log(VCV[[2]][caracteres.medidos,caracteres.medidos])+(numero.caracteres.medidos*log(2*pi)))
	fenotipo.medido.DiMa.2 <- -2*(logv.fenotipo.medido[2]+logv.fenotipo.medido.C.2)
	P.fenotipo.medido.al.menos.tan.extremo[2] <- pchisq(fenotipo.medido.DiMa.2, numero.caracteres.medidos, lower.tail = FALSE, log.p = FALSE)
	P.fenotipo.medido.al.menos.tan.extremo[2]

	logv.fenotipo.medido.C.3 <- 0.5*(log(VCV[[3]][caracteres.medidos,caracteres.medidos])+(numero.caracteres.medidos*log(2*pi)))
	fenotipo.medido.DiMa.3 <- -2*(logv.fenotipo.medido[3]+logv.fenotipo.medido.C.3)
	P.fenotipo.medido.al.menos.tan.extremo[3] <- pchisq(fenotipo.medido.DiMa.3, numero.caracteres.medidos, lower.tail = FALSE, log.p = FALSE)
	P.fenotipo.medido.al.menos.tan.extremo[3]
}

#si el fenotipo medido incluye dos o más caracteres
if(numero.caracteres.medidos>1){
	logv.fenotipo.medido.C.1 <- 0.5*(log(det(VCV[[1]][caracteres.medidos,caracteres.medidos]))+(numero.caracteres.medidos*log(2*pi)))
	fenotipo.medido.DiMa.1 <- sqrt(-2*(logv.fenotipo.medido[1]+logv.fenotipo.medido.C.1))
	P.fenotipo.medido.al.menos.tan.extremo[1] <- pchisq(fenotipo.medido.DiMa.1, p, lower.tail = FALSE, log.p = FALSE)
	P.fenotipo.medido.al.menos.tan.extremo[1]

	logv.fenotipo.medido.C.2 <- 0.5*(log(det(VCV[[2]][caracteres.medidos,caracteres.medidos]))+(numero.caracteres.medidos*log(2*pi)))
	fenotipo.medido.DiMa.2 <- sqrt(-2*(logv.fenotipo.medido[2]+logv.fenotipo.medido.C.2))
	P.fenotipo.medido.al.menos.tan.extremo[2] <- pchisq(fenotipo.medido.DiMa.2, p, lower.tail = FALSE, log.p = FALSE)
	P.fenotipo.medido.al.menos.tan.extremo[2]

	logv.fenotipo.medido.C.3 <- 0.5*(log(det(VCV[[3]][caracteres.medidos,caracteres.medidos]))+(numero.caracteres.medidos*log(2*pi)))
	fenotipo.medido.DiMa.3 <- sqrt(-2*(logv.fenotipo.medido[3]+logv.fenotipo.medido.C.3))
	P.fenotipo.medido.al.menos.tan.extremo[3] <- pchisq(fenotipo.medido.DiMa.3, p, lower.tail = FALSE, log.p = FALSE)
	P.fenotipo.medido.al.menos.tan.extremo[3]
}


################################################################################
# 13. Comparación las probabilidades de muestrear un fenotipo y un fenotipo medido
#     al menos tan extremos como los observados.
################################################################################

plot(P.fenotipo.al.menos.tan.extremo, P.fenotipo.medido.al.menos.tan.extremo, pch=19, cex=2, col=color.especies,
	xlab=paste("P ( fenotipo al menos tan extremo como el observado ), p=", p, sep=""), 
	ylab=paste("P ( fenotipo medido al menos tan extremo como el observado ), ", numero.caracteres.medidos, " carácteres", sep=""))
#legend("bottomright", c("Especie 1", "Especie 2", "Especie 3"), pch=19, cex=2, col=color.especies)
#legend("topright", c("Especie 1", "Especie 2", "Especie 3"), pch=19, cex=2, col=color.especies)
#legend("bottomleft", c("Especie 1", "Especie 2", "Especie 3"), pch=19, cex=2, col=color.especies)
#legend("topleft", c("Especie 1", "Especie 2", "Especie 3"), pch=19, cex=2, col=color.especies)


