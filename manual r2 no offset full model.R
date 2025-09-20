# Carica le librerie necessarie
library(dplyr)
library(readxl)
library(glmmTMB)
library(MuMIn)

# Carica il dataset
file_path <- "Desktop/species/mountainhare3.xlsx"
data <- read_xlsx(file_path)

# Filtra i dati
data <- filter(data, snow_days <= 5)

# Crea le variabili necessarie una per una
data$triangles <- as.factor(data$triangle_id)
data$year_f <- as.factor(data$year)
data$riistakeskus <- as.factor(data$riistakeskus)
data$rhy_id <- as.factor(data$rhy_id)
data$year.c <- data$year - mean(data$year, na.rm = TRUE)
data$latitude.100km <- (data$triangle_y - mean(data$triangle_y, na.rm = TRUE)) / 100000
data$longitude.100km <- (data$triangle_x - mean(data$triangle_x, na.rm = TRUE)) / 100000
data$FLII <- data$indexvalue - mean(data$indexvalue, na.rm = TRUE)
data$forest_cover <- data$forest_pro - mean(data$forest_pro, na.rm = TRUE)
data$canopy <- data$canopy - mean(data$canopy, na.rm = TRUE)
data$ln.snow.days <- log(data$snow_days)
data$ln.covered.di <- log(data$covered_di / 10000)

# Rimuove i valori mancanti
data <- na.omit(data)

# Crea il modello full con glmmTMB
model_nb <- glmmTMB(n_tracks ~ FLII + forest_cover + canopy + latitude.100km + year_f + ln.snow.days + 
                      (1 | riistakeskus/rhy_id/triangles),
                    offset = ln.covered.di,
                    family = nbinom2(link = "log"),
                    data = data)

# Stampa il riassunto del modello
print(summary(model_nb))

# Controlla se ln.covered.di è presente nel data.frame del modello
if ("ln.covered.di" %in% colnames(model_nb$frame)) {
  print("ln.covered.di is present in the model's data frame.")
} else {
  print("ln.covered.di is NOT present in the model's data frame.")
}

# Funzione per calcolare R2 manualmente
glmmTMB.R2 <- function(m, incl_DSV = FALSE, incl_offset = FALSE, offset_var = NULL, exclude = NULL) {
  # Data frame per la predizione
  new.data.m <- m$frame
  
  # Aggiungi l'offset al data frame se necessario
  if (!is.null(offset_var) & incl_offset == FALSE) {
    new.data.m[[offset_var]] <- data[[offset_var]]
  }
  
  # Definisci i valori possibili da escludere dalla predizione degli effetti fissi
  if (!is.null(exclude)) {
    new.data.m[, match(exclude, colnames(new.data.m))] <- 0
  }
  
  # Valori previsti sulla scala logaritmica
  ln.yhat <- predict(m, type = "link", re.form = ~0, newdata = new.data.m)
  n <- length(ln.yhat)
  
  # Rimuovi l'offset dagli effetti fissi se necessario
  if (!is.null(offset_var) & incl_offset == FALSE) {
    ln.yhat <- ln.yhat - new.data.m[[offset_var]]
  }
  
  # Calcola le componenti di varianza
  s2_f <- var(ln.yhat) * (n-1) / n
  s2_a <- unlist(VarCorr(m)$cond)
  s2_e <- trigamma(sigma(m))
  s2_d <- log(1 / exp(mean(ln.yhat)) + 1)
  
  # Calcola R2, con o senza la varianza specifica della distribuzione
  if (incl_DSV == TRUE) {
    s2_tot <- s2_f + sum(s2_a) + s2_e + s2_d
    R2_marg <- s2_f / s2_tot
    R2_cond <- (s2_f + sum(s2_a)) / s2_tot
  } else {
    s2_tot <- s2_f + sum(s2_a) + s2_e
    R2_marg <- s2_f / s2_tot
    R2_cond <- (s2_f + sum(s2_a)) / s2_tot
  }
  
  # Restituisce il risultato
  R2 <- c(R2_marg = unname(R2_marg), R2_cond = unname(R2_cond))
  return(R2)
}

# Calcola R2 manualmente
r2_custom <- glmmTMB.R2(model_nb, incl_offset = FALSE, offset_var = "ln.covered.di", exclude = "ln.snow.days")
print(r2_custom)

# Calcola R2 usando MuMIn
r2_mumin <- r.squaredGLMM(model_nb)
cat("R2 MuMIn:\n")
print(r2_mumin)
