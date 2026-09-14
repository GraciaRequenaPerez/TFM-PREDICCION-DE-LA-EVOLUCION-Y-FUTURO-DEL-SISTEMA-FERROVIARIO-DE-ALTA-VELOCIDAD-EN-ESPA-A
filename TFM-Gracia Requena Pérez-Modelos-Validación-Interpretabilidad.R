library(readxl)
library(forecast)
library(lmtest)

# Carga de datos
Viajeros <- read_excel("C:/Users/graci/Desktop/TFM/Datos/DATOS AVE y EXO.xlsx")

# 1. Serie temporal principal (Variable Objetivo)
y_ts <- ts(Viajeros$ViajerosAVyLD, start = c(1996, 1), frequency = 12)

# 2. Matriz de Variables Exógenas
x_mat <- as.matrix(Viajeros[, 3:14])
x_ts  <- ts(x_mat, start = c(1996, 1), frequency = 12)


##############################################################################
##             1. PARTICIÓN DE DATOS (UNIFICADA Y ESTÁNDAR)                ##
##############################################################################

# --- PARTICIÓN A: Modelos Tradicionales / Directos (HW, SARIMA) ---
y_train <- window(y_ts, start = c(1996, 1), end = c(2023, 12))
y_test  <- window(y_ts, start = c(2024, 1), end = c(2025, 12))

# --- PARTICIÓN B: Modelos con Validación/Ajuste (Prophet, XGBoost, ML) ---
# Subdivisión con conjunto de Validación
y_train_sub <- window(y_ts, start = c(1996, 1), end = c(2021, 12))
y_val       <- window(y_ts, start = c(2022, 1), end = c(2023, 12))

# Unión de Train + Val para el reentrenamiento final antes de predecir Test
y_train_val <- window(y_ts, start = c(1996, 1), end = c(2023, 12)) 

# --- Exógenas Alineadas con Partición B ---
x_train     <- window(x_ts, start = c(1996, 1), end = c(2021, 12))
x_val       <- window(x_ts, start = c(2022, 1), end = c(2023, 12))
x_train_val <- window(x_ts, start = c(1996, 1), end = c(2023, 12))
x_test      <- window(x_ts, start = c(2024, 1), end = c(2025, 12))


##############################################################################
##                    2. FUNCIÓN DE CÁLCULO DE MÉTRICAS                      ##
##############################################################################

calc_metrics <- function(real, pred) {
  real <- as.numeric(real)
  pred <- as.numeric(pred)
  
  mae  <- mean(abs(real - pred))
  rmse <- sqrt(mean((real - pred)^2))
  mape <- mean(abs((real - pred) / real)) * 100
  
  data.frame(
    MAE  = round(mae, 2),
    RMSE = round(rmse, 2),
    MAPE = round(mape, 2)
  )
}


##############################################################################
##                    MODELO 1: HOLT-WINTERS (ADITIVO)                      ##
##############################################################################

# ==============================================================================
# MODELO 1: HOLT-WINTERS (ADITIVO)
# ==============================================================================

fit_hw <- HoltWinters(y_train, seasonal = "additive")

# Extracción de parámetros
hw_params <- data.frame(
  Parametro = c("Alpha (Nivel)", "Beta (Tendencia)", "Gamma (Estacionalidad)"),
  Valor     = c(fit_hw$alpha, fit_hw$beta, fit_hw$gamma)
)
hw_params

# Diagnóstico de residuos
checkresiduals(fit_hw)

# Predicción en Test y evaluación
pred_hw_test <- forecast(fit_hw, h = 24)
metrics_hw_test <- calc_metrics(real = y_test, pred = pred_hw_test$mean)
metrics_hw_test

# ==============================================================================
# MODELO 2: SARIMA
# ==============================================================================
library(tseries)
library(forecast)

# 1. Diagnóstico de Estacionariedad y Diferenciación
adf.test(y_train)
ndiffs(y_train)   
nsdiffs(y_train)
var(y_train, na.rm = TRUE)

# Primera diferenciación regular (d=1)
y_diff <- diff(y_train, lag = 1)
adf.test(y_diff)
var(y_diff, na.rm = TRUE)

par(mfrow = c(1, 2))
Acf(y_diff, main = "ACF (d=1)")
Pacf(y_diff, main = "PACF (d=1)")

# Diferenciación regular y estacional (d=1, D=1)
y_train_diff <- diff(diff(y_train, lag = 1), lag = 12)
adf.test(y_train_diff)
var(y_train_diff)

Acf(y_train_diff, main = "ACF (d=1, D=1)")
Pacf(y_train_diff, main = "PACF (d=1, D=1)")
par(mfrow = c(1, 1)) # Restablecer ventana gráfica


# 2. Evaluación de candidatos SARIMA
sarima_models <- list(
  "SARIMA(0,1,1)(0,1,1)[12]" = Arima(y_train, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12)),
  "SARIMA(1,1,1)(0,1,1)[12]" = Arima(y_train, order = c(1, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12)),
  "SARIMA(1,1,1)(1,1,1)[12]" = Arima(y_train, order = c(1, 1, 1), seasonal = list(order = c(1, 1, 1), period = 12)),
  "SARIMA(1,1,0)(0,1,1)[12]" = Arima(y_train, order = c(1, 1, 0), seasonal = list(order = c(0, 1, 1), period = 12)),
  "SARIMA(0,1,1)(0,1,2)[12]" = Arima(y_train, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 2), period = 12))
)

# Comparativa de AICc automatizada
aicc_tabla <- data.frame(
  Modelo = names(sarima_models),
  AICc   = sapply(sarima_models, function(m) m$aicc)
)
aicc_tabla


# 3. Selección y validación del modelo final
fit_sarima <- sarima_models[["SARIMA(0,1,1)(0,1,2)[12]"]]

coeftest(fit_sarima)
checkresiduals(fit_sarima)


# 4. Predicción en Test y evaluación
pred_sarima_test <- forecast(fit_sarima, h = 24)
metrics_sarima_test <- calc_metrics(real = y_test, pred = pred_sarima_test$mean)
metrics_sarima_test

# ==============================================================================
# MODELO 3: SARIMAX
# ==============================================================================
library(forecast)
library(lmtest)
library(car)

# 1. Preparación de datasets específicos (smax_*) alineados (1996-2023)
smax_train_df <- data.frame(Y = as.numeric(y_train), x_train_val)
smax_train_x  <- x_train_val
smax_test_x   <- x_test


# 2. Paso previo: MRL y selección de exógenas por VIF
mrl_inicial <- lm(Y ~ ., data = smax_train_df)
summary(mrl_inicial)
vif(mrl_inicial)

# Función de filtrado por VIF silenciosa
filtrar_por_vif <- function(datos, var_y = "Y", umbral = 10) {
  df <- datos
  
  repeat {
    mod <- lm(as.formula(paste(var_y, "~ .")), data = df)
    
    coefs <- summary(mod)$coefficients
    if (any(is.na(coefs[, 1]))) {
      var_na <- names(which(is.na(coefs[, 1])))[1]
      df <- df[, !names(df) %in% var_na]
      next
    }
    
    vifs <- tryCatch({
      v <- vif(mod)
      if (is.matrix(v)) v[, 3]^2 else v
    }, error = function(e) NULL)
    
    if (is.null(vifs) || length(vifs) <= 1) break
    
    vif_max <- max(vifs, na.rm = TRUE)
    if (vif_max <= umbral) break
    
    var_eliminar <- names(which.max(vifs))
    df <- df[, !names(df) %in% var_eliminar]
  }
  
  return(df)
}

smax_train_df_limpio <- filtrar_por_vif(smax_train_df, umbral = 10)
names(smax_train_df_limpio)[-1]

# MRL depurado y selección Backward
mrl_limpio   <- lm(Y ~ ., data = smax_train_df_limpio)
mrl_depurado <- step(mrl_limpio, direction = "backward", trace = FALSE)

names(coef(mrl_depurado))[-1]
vars_finales <- c("FestivosN", "COVID")

# Filtrado definitivo de exógenas
smax_x_train_sel <- smax_train_x[, vars_finales]
smax_x_test_sel  <- smax_test_x[, vars_finales]

mrl_limpio_final <- lm(Y ~ ., data = smax_train_df_limpio[, c("Y", vars_finales)])
summary(mrl_limpio_final)


# 3. Extracción de residuos y diagnóstico de diferenciación (d, D)
residuos_mrl <- ts(residuals(mrl_limpio_final), start = start(y_train), frequency = 12)
#### Estimación del parámetro lambda de Box-Cox para saber si es aditiva o multiplicativa ###
lambdaR <- BoxCox.lambda(residuos_mrl)
cat("Valor de Lambda estimado:", round(lambdaR, 3), "\n")

d_sugerida <- ndiffs(residuos_mrl)
D_sugerida <- nsdiffs(residuos_mrl)

# Inspección de varianzas según nivel de diferenciación
v_nivel   <- var(residuos_mrl, na.rm = TRUE)
v_d1      <- var(diff(residuos_mrl, lag = 1), na.rm = TRUE)
v_d1_D1   <- var(diff(diff(residuos_mrl, lag = 1), lag = 12), na.rm = TRUE)
c(Varianza_Nivel = v_nivel, Varianza_d1 = v_d1, Varianza_d1_D1 = v_d1_D1)

x_diff_eval <- diff(diff(as.matrix(smax_x_train_sel), lag = 1), lag = 12)
sd_exogenas <- apply(x_diff_eval, 2, sd, na.rm = TRUE)
print(sd_exogenas)

# Análisis gráfico de ACF/PACF
par(mfrow = c(1, 2))
Acf(residuos_mrl, main = "ACF Residuos Nivel")
Pacf(residuos_mrl, main = "PACF Residuos Nivel")

residuos_diff <- diff(diff(residuos_mrl, lag = 1), lag = 12)
Acf(residuos_diff, main = "ACF Residuos (d=1, D=1)")
Pacf(residuos_diff, main = "PACF Residuos (d=1, D=1)")
par(mfrow = c(1, 1))


# 4. Ajuste de modelos SARIMAX candidatos (y_train y smax_x_train_sel alineados 1996-2023)
m1 <- Arima(y_train, order = c(0,1,1), seasonal = list(order = c(0,1,1), period = 12), xreg = smax_x_train_sel)
m2 <- Arima(y_train, order = c(2,1,1), seasonal = list(order = c(0,1,1), period = 12), xreg = smax_x_train_sel)
m3 <- Arima(y_train, order = c(2,1,1), seasonal = list(order = c(0,1,2), period = 12), xreg = smax_x_train_sel)
m4 <- Arima(y_train, order = c(2,1,2), seasonal = list(order = c(0,1,2), period = 12), xreg = smax_x_train_sel)

# Tabla comparativa de AICc
aicc_tabla_sarimax <- data.frame(
  Modelo = c("M1: (0,1,1)(0,1,1)", "M2: (2,1,1)(0,1,1)", "M3: (2,1,1)(0,1,2)", "M4: (2,1,2)(0,1,2)"),
  AICc   = c(m1$aicc, m2$aicc, m3$aicc, m4$aicc)
)
aicc_tabla_sarimax


# 5. Modelo final, validación y predicción
fit_sarimax <- m4  # Seleccionamos directamente el objeto m4

summary(fit_sarimax)
coeftest(fit_sarimax)
checkresiduals(fit_sarimax)

# Predicción en Test y evaluación
pred_sarimax_test <- forecast(fit_sarimax, xreg = smax_x_test_sel)
metrics_sarimax_test <- calc_metrics(real = y_test, pred = pred_sarimax_test$mean)
metrics_sarimax_test

# ==============================================================================
# MODELO 4: PROPHET (UNIVARIANTE Y MULTIVARIANTE OPTIMIZADO)
# ==============================================================================

library(prophet)
library(dplyr)

# ----------------------------------------------------------------------------
# 3.0 FUNCIÓN AUXILIAR DE CONVERSIÓN A DATAFRAME PROPHET
# ----------------------------------------------------------------------------
crear_df_prophet <- function(y_ts, xreg_mat = NULL) {
  time_vals <- time(y_ts)
  years     <- floor(time_vals)
  months    <- round((time_vals - years) * 12) + 1
  fechas    <- as.Date(sprintf("%d-%02d-01", years, months))
  
  df <- data.frame(
    ds = fechas,
    y  = as.numeric(y_ts)
  )
  
  if (!is.null(xreg_mat)) {
    xreg_df <- as.data.frame(xreg_mat)
    df      <- cbind(df, xreg_df)
  }
  
  return(df)
}

# ----------------------------------------------------------------------------
# 3.1 PROPHET BASE (UNIVARIANTE) - CONFIGURACIÓN POR DEFECTO
# ----------------------------------------------------------------------------
pro_df_train_val_no_exo <- crear_df_prophet(y_train_val)
pro_df_test_no_exo      <- crear_df_prophet(y_test)

fit_prophet_base <- prophet(
  yearly.seasonality = TRUE,
  weekly.seasonality = FALSE,
  daily.seasonality  = FALSE,
  seasonality.mode   = "additive"
)

suppressMessages(suppressWarnings({
  fit_prophet_base <- fit.prophet(fit_prophet_base, pro_df_train_val_no_exo)
}))

pred_prophet_base_test <- predict(fit_prophet_base, pro_df_test_no_exo)

metrics_prophet_base_test <- calc_metrics(
  real = pro_df_test_no_exo$y, 
  pred = pred_prophet_base_test$yhat
)

print(metrics_prophet_base_test)


# ----------------------------------------------------------------------------
# 3.2 PROPHET CON EXÓGENAS Y TUNING DE HIPERPARÁMETROS
# ----------------------------------------------------------------------------
pro_df_train_exo     <- crear_df_prophet(y_train_sub, x_train)
pro_df_val_exo       <- crear_df_prophet(y_val, x_val)
pro_df_train_val_exo <- crear_df_prophet(y_train_val, x_train_val)
pro_df_test_exo      <- crear_df_prophet(y_test, x_test)

nombres_exogenas <- colnames(x_train)

# --- Malla de Búsqueda (Grid Search) ---
param_grid <- expand.grid(
  changepoint_prior_scale = c(0.3, 0.5, 0.7, 0.9, 0.95),
  seasonality_prior_scale = c(0.005, 0.01, 0.05),
  stringsAsFactors        = FALSE
)

results_grid <- data.frame()

for (i in 1:nrow(param_grid)) {
  p <- param_grid[i, ]
  
  m_tune <- prophet(
    changepoint.prior.scale = p$changepoint_prior_scale,
    seasonality.prior.scale = p$seasonality_prior_scale,
    seasonality.mode        = "additive",
    yearly.seasonality      = TRUE,
    weekly.seasonality      = FALSE,
    daily.seasonality       = FALSE
  )
  
  for (col in nombres_exogenas) {
    m_tune <- add_regressor(m_tune, name = col)
  }
  
  suppressMessages(suppressWarnings({
    m_tune <- fit.prophet(m_tune, pro_df_train_exo)
  }))
  
  pred_val    <- predict(m_tune, pro_df_val_exo)
  metrics_val <- calc_metrics(real = pro_df_val_exo$y, pred = pred_val$yhat)
  
  res          <- cbind(p, metrics_val)
  results_grid <- rbind(results_grid, res)
}

# Selección del mejor hiperparámetro (menor RMSE)
results_grid <- results_grid[order(results_grid$RMSE), ]
best_params  <- results_grid[1, ]

print(best_params[, c("changepoint_prior_scale", "seasonality_prior_scale", "MAE", "RMSE", "MAPE")])

# --- Reentrenamiento Modelo Optimizado en Train + Val ---
fit_prophet_exo_tuned <- prophet(
  changepoint.prior.scale = best_params$changepoint_prior_scale,
  seasonality.prior.scale = best_params$seasonality_prior_scale,
  seasonality.mode        = "additive",
  yearly.seasonality      = TRUE,
  weekly.seasonality      = FALSE,
  daily.seasonality       = FALSE
)

for (col in nombres_exogenas) {
  fit_prophet_exo_tuned <- add_regressor(fit_prophet_exo_tuned, name = col)
}

suppressMessages(suppressWarnings({
  fit_prophet_exo_tuned <- fit.prophet(fit_prophet_exo_tuned, pro_df_train_val_exo)
}))

pred_prophet_exo_tuned <- predict(fit_prophet_exo_tuned, pro_df_test_exo)

metrics_prophet_exo_tuned <- calc_metrics(
  real = pro_df_test_exo$y, 
  pred = pred_prophet_exo_tuned$yhat
)

print(metrics_prophet_exo_tuned)

# ----------------------------------------------------------------------------
# 3.3 COMPARATIVA RÁPIDA PROPHET
# ----------------------------------------------------------------------------
comparativa_prophet <- rbind(
  Base_Univariante = metrics_prophet_base_test,
  Tuned_Exogenas   = metrics_prophet_exo_tuned
)

print(comparativa_prophet)

##############################################################################
##                          MODELO 5: XGBOOST                               ##
##############################################################################

library(xgboost)
library(dplyr)

# ----------------------------------------------------------------------------
# 5.0 FUNCIÓN AUXILIAR DE PREPARACIÓN DE DATOS (LAGGED + EXÓGENAS + CALENDARIO)
# ----------------------------------------------------------------------------
preparar_datos_xgboost <- function(y_ts, xreg_mat = NULL, num_lags = 12) {
  time_vals <- time(y_ts)
  years     <- floor(time_vals)
  months    <- round((time_vals - years) * 12) + 1
  
  df <- data.frame(
    Fecha = as.Date(sprintf("%d-%02d-01", years, months)),
    Mes   = as.numeric(months),
    y     = as.numeric(y_ts)
  )
  
  # Generar Lags de la variable objetivo (y_t-1 a y_t-12)
  for (l in 1:num_lags) {
    df[[paste0("lag_", l)]] <- dplyr::lag(df$y, n = l)
  }
  
  # Añadir variables exógenas con nombres limpios
  if (!is.null(xreg_mat)) {
    xreg_df <- as.data.frame(xreg_mat)
    colnames(xreg_df) <- make.names(colnames(xreg_df))
    df <- cbind(df, xreg_df)
  }
  
  # Eliminar filas iniciales generadas por los lags
  df <- na.omit(df)
  return(df)
}

# ----------------------------------------------------------------------------
# 5.1 CONSTRUCCIÓN Y DIVISIÓN DEL DATASET XGBOOST
# ----------------------------------------------------------------------------

# Crear dataset consolidado
df_xgb_completo <- preparar_datos_xgboost(y_ts, x_mat, num_lags = 12)

# Filtrar subconjuntos manteniendo las ventanas temporales exactas
df_train_xgb     <- df_xgb_completo %>% filter(Fecha <= as.Date("2021-12-01"))
df_val_xgb       <- df_xgb_completo %>% filter(Fecha >= as.Date("2022-01-01") & Fecha <= as.Date("2023-12-01"))
df_test_xgb      <- df_xgb_completo %>% filter(Fecha >= as.Date("2024-01-01") & Fecha <= as.Date("2025-12-01"))
df_train_val_xgb <- df_xgb_completo %>% filter(Fecha <= as.Date("2023-12-01"))

# Definir predictores (features) excluyendo clave temporal y objetivo
features_xgb <- setdiff(colnames(df_train_xgb), c("Fecha", "y"))

# Extracción de matrices/vectores para XGBoost
xgb_X_train     <- as.matrix(df_train_xgb[, features_xgb])
xgb_y_train     <- df_train_xgb$y

xgb_X_val       <- as.matrix(df_val_xgb[, features_xgb])
xgb_y_val       <- df_val_xgb$y

xgb_X_test      <- as.matrix(df_test_xgb[, features_xgb])
xgb_y_test      <- df_test_xgb$y

xgb_X_train_val <- as.matrix(df_train_val_xgb[, features_xgb])
xgb_y_train_val <- df_train_val_xgb$y

# Construcción de objetos DMatrix
xgb_dtrain     <- xgb.DMatrix(data = xgb_X_train, label = xgb_y_train)
xgb_dval       <- xgb.DMatrix(data = xgb_X_val, label = xgb_y_val)
xgb_dtrain_val <- xgb.DMatrix(data = xgb_X_train_val, label = xgb_y_train_val)
xgb_dtest      <- xgb.DMatrix(data = xgb_X_test, label = xgb_y_test)

# ----------------------------------------------------------------------------
# 5.2 MODELO BASE XGBOOST
# ----------------------------------------------------------------------------
params_base <- list(
  booster          = "gbtree",
  objective        = "reg:squarederror",
  eta              = 0.3,
  max_depth        = 6,
  subsample        = 0.8,
  colsample_bytree = 0.8)

# 1. Evaluación rápida en Validación (2022-2023)
set.seed(123)
fit_xgb_val <- xgb.train(
  params  = params_base,
  data    = xgb_dtrain,
  nrounds = 250,
  verbose = 0)

pred_xgb_val    <- predict(fit_xgb_val, xgb_dval)
metrics_xgb_val <- calc_metrics(real = xgb_y_val, pred = pred_xgb_val)

print(metrics_xgb_val)

# Importancia de Variables (Modelo Base)
importance_matrix_base <- xgb.importance(feature_names = features_xgb, model = fit_xgb_val)
xgb.plot.importance(importance_matrix_base[1:min(10, nrow(importance_matrix_base)), ], 
                    main = "Top 10 Importancia de Variables (XGBoost Base)")

# 2. Evaluación del modelo Base en TEST (Reentrenando en Train + Val)
set.seed(123)
fit_xgb_base_test <- xgb.train(
  params  = params_base,
  data    = xgb_dtrain_val,
  nrounds = 250,
  verbose = 0)

pred_xgb_base_test    <- predict(fit_xgb_base_test, xgb_dtest)
metrics_xgb_base_test <- calc_metrics(real = xgb_y_test, pred = pred_xgb_base_test)

print(metrics_xgb_base_test)

# ----------------------------------------------------------------------------
# 5.3 OPTIMIZACIÓN DE HIPERPARÁMETROS (GRID SEARCH)
# ----------------------------------------------------------------------------
grid_xgb <- expand.grid(
  eta              = c(0.3, 0.4, 0.5),
  max_depth        = c(3, 4, 5),
  subsample        = c(0.6, 0.7, 0.8),
  colsample_bytree = c(0.6, 0.7, 0.8),
  nrounds          = c(250, 300, 350),
  stringsAsFactors = FALSE)

results_xgb_grid <- data.frame()

for (i in 1:nrow(grid_xgb)) {
  p <- grid_xgb[i, ]
  
  params_tmp <- list(
    booster          = "gbtree",
    objective        = "reg:squarederror",
    eta              = p$eta,
    max_depth        = p$max_depth,
    subsample        = p$subsample,
    colsample_bytree = p$colsample_bytree)
  
  set.seed(123)
  fit_tmp <- xgb.train(
    params  = params_tmp,
    data    = xgb_dtrain,
    nrounds = p$nrounds,
    verbose = 0 )
  
  pred_tmp    <- predict(fit_tmp, xgb_dval)
  metrics_tmp <- calc_metrics(real = xgb_y_val, pred = pred_tmp)
  
  res              <- cbind(p, metrics_tmp)
  results_xgb_grid <- rbind(results_xgb_grid, res)
}

# Selección del mejor set de hiperparámetros por RMSE
results_xgb_grid <- results_xgb_grid[order(results_xgb_grid$RMSE), ]
top5_xgb_params  <- results_xgb_grid[1:5, ]

print(top5_xgb_params[, c("eta", "max_depth", "subsample", "colsample_bytree", "nrounds", "MAE", "RMSE", "MAPE")])

best_xgb_params <- results_xgb_grid[1, ]

print(best_xgb_params[, c("eta", "max_depth", "subsample", "colsample_bytree", "nrounds", "MAE", "RMSE", "MAPE")])

# ----------------------------------------------------------------------------
# 5.4 REENTRENAMIENTO Y EVALUACIÓN DEFINITIVA (TEST)
# ----------------------------------------------------------------------------
params_optimos <- list(
  booster          = "gbtree",
  objective        = "reg:squarederror",
  eta              = best_xgb_params$eta,
  max_depth        = best_xgb_params$max_depth,
  subsample        = best_xgb_params$subsample,
  colsample_bytree = best_xgb_params$colsample_bytree)

set.seed(123)
fit_xgb_tuned <- xgb.train(
  params  = params_optimos,
  data    = xgb_dtrain_val,
  nrounds = best_xgb_params$nrounds,
  verbose = 0)

# Predicción final en TEST
pred_xgb_tuned    <- predict(fit_xgb_tuned, xgb_dtest)
metrics_xgb_tuned <- calc_metrics(real = xgb_y_test, pred = pred_xgb_tuned)

print(metrics_xgb_tuned)

# Importancia de Variables (Modelo Tuned Final)
importance_matrix_tuned <- xgb.importance(feature_names = features_xgb, model = fit_xgb_tuned)
xgb.plot.importance(importance_matrix_tuned[1:min(10, nrow(importance_matrix_tuned)), ], 
                    main = "Top 10 Importancia de Variables (XGBoost Tuned)")

# ----------------------------------------------------------------------------
# 5.5 COMPARATIVA RÁPIDA XGBOOST
# ----------------------------------------------------------------------------
comparativa_xgb <- rbind(
  XGBoost_Base_Test  = metrics_xgb_base_test,
  XGBoost_Tuned_Test = metrics_xgb_tuned)
print(comparativa_xgb)

##############################################################################
##                            MODELO 6: MLP                                 ##
##############################################################################

library(caret)
library(nnet)
library(dplyr)

# ==============================================================================
# 6.1 MODELO BASELINE (ENTRADA ESTÁNDAR: SOLO EXÓGENAS)
# ==============================================================================
mlp_preProc_exog <- preProcess(x_train_val, method = c("center", "scale"))
mlp_x_tr_exog_sc <- predict(mlp_preProc_exog, x_train_val)
mlp_x_te_exog_sc <- predict(mlp_preProc_exog, x_test)

mlp_df_exog_tr   <- data.frame(mlp_x_tr_exog_sc, ViajerosAVyLD = as.numeric(y_train_val))

set.seed(123)
mlp_fit_base <- nnet(
  ViajerosAVyLD ~ .,
  data   = mlp_df_exog_tr,
  size   = 3,
  decay  = 0.1,
  linout = TRUE,
  trace  = FALSE,
  maxit  = 1500)

mlp_pred_base_test <- predict(mlp_fit_base, newdata = as.data.frame(mlp_x_te_exog_sc))
metrics_mlp_base   <- calc_metrics(real = as.numeric(y_test), pred = as.numeric(mlp_pred_base_test))

print(metrics_mlp_base)

# ==============================================================================
# 6.2 CONSTRUCCIÓN DE CONJUNTOS CON RETARDOS
# ==============================================================================
mlp_num_lags <- 12

crear_lags_matriz <- function(vec, num_lags = 12) {
  n <- length(vec)
  mat_lags <- matrix(NA, nrow = n, ncol = num_lags)
  for (i in 1:num_lags) {
    mat_lags[(i + 1):n, i] <- vec[1:(n - i)]
  }
  colnames(mat_lags) <- paste0("lag_", 1:num_lags)
  return(mat_lags)
}

mlp_y_completo <- c(as.numeric(y_train_val), as.numeric(y_test))
mlp_x_completo <- rbind(x_train_val, x_test)

mlp_mat_lags   <- crear_lags_matriz(mlp_y_completo, num_lags = mlp_num_lags)
mlp_X_completo <- cbind(mlp_mat_lags, mlp_x_completo)

mlp_idx_validos <- (mlp_num_lags + 1):length(mlp_y_completo)
mlp_X_clean     <- mlp_X_completo[mlp_idx_validos, , drop = FALSE]
mlp_y_clean     <- mlp_y_completo[mlp_idx_validos]

mlp_n_tr_sub_clean <- length(y_train_sub) - mlp_num_lags
mlp_n_tr_val_clean <- length(y_train_val) - mlp_num_lags

# Particiones con prefijo propio
mlp_x_train_sub <- mlp_X_clean[1:mlp_n_tr_sub_clean, , drop = FALSE]
mlp_y_train_sub <- mlp_y_clean[1:mlp_n_tr_sub_clean]

mlp_x_val       <- mlp_X_clean[(mlp_n_tr_sub_clean + 1):mlp_n_tr_val_clean, , drop = FALSE]
mlp_y_val       <- mlp_y_clean[(mlp_n_tr_sub_clean + 1):mlp_n_tr_val_clean]

mlp_x_train_val <- mlp_X_clean[1:mlp_n_tr_val_clean, , drop = FALSE]
mlp_y_train_val <- mlp_y_clean[1:mlp_n_tr_val_clean]

mlp_x_test      <- mlp_X_clean[(mlp_n_tr_val_clean + 1):nrow(mlp_X_clean), , drop = FALSE]
mlp_y_test      <- mlp_y_clean[(mlp_n_tr_val_clean + 1):length(mlp_y_clean)]

# Escalado Z-Score para modelo con retardos
mlp_preProc_sub  <- preProcess(mlp_x_train_sub, method = c("center", "scale"))
mlp_X_tr_scaled  <- predict(mlp_preProc_sub, mlp_x_train_sub)
mlp_X_val_scaled <- predict(mlp_preProc_sub, mlp_x_val)
mlp_df_train     <- data.frame(mlp_X_tr_scaled, ViajerosAVyLD = as.numeric(mlp_y_train_sub))

mlp_preProc_tr_val  <- preProcess(mlp_x_train_val, method = c("center", "scale"))
mlp_X_tr_val_scaled <- predict(mlp_preProc_tr_val, mlp_x_train_val)
mlp_X_test_scaled   <- predict(mlp_preProc_tr_val, mlp_x_test)
mlp_df_tr_val       <- data.frame(mlp_X_tr_val_scaled, ViajerosAVyLD = as.numeric(mlp_y_train_val))

# ==============================================================================
# 6.3 BÚSQUEDA EN MALLA (OPTIMIZADO POR RMSE) Y EVALUACIÓN EN TEST
# ==============================================================================
mlp_grid <- expand.grid(
  size  = c(2, 4, 6, 8),
  decay = c(0.001, 0.01, 0.1))

mlp_resultados_malla <- data.frame()

for (i in 1:nrow(mlp_grid)) {
  s <- mlp_grid$size[i]
  d <- mlp_grid$decay[i]
  
  set.seed(123)
  fit_temp <- nnet(
    ViajerosAVyLD ~ .,
    data   = mlp_df_train,
    size   = s,
    decay  = d,
    linout = TRUE,
    trace  = FALSE,
    maxit  = 1500)
  
  pred_val <- predict(fit_temp, newdata = as.data.frame(mlp_X_val_scaled))
  m_val    <- calc_metrics(real = as.numeric(mlp_y_val), pred = as.numeric(pred_val))
  
  mlp_resultados_malla <- rbind(
    mlp_resultados_malla, 
    data.frame(size = s, decay = d, MAE = m_val$MAE, RMSE = m_val$RMSE, MAPE = m_val$MAPE))
}

# Selección por menor RMSE
mlp_best_params <- mlp_resultados_malla %>% arrange(RMSE) %>% slice(1)

print(mlp_best_params)

# --- Reentrenamiento Modelo Optimizado en Train + Val ---
set.seed(123)
mlp_model_opt <- nnet(
  ViajerosAVyLD ~ .,
  data   = mlp_df_tr_val,
  size   = mlp_best_params$size,
  decay  = mlp_best_params$decay,
  linout = TRUE,
  trace  = FALSE,
  maxit  = 1500)

mlp_pred_opt_test <- predict(mlp_model_opt, newdata = as.data.frame(mlp_X_test_scaled))
metrics_mlp_opt   <- calc_metrics(real = as.numeric(mlp_y_test), pred = as.numeric(mlp_pred_opt_test))

print(metrics_mlp_opt)

# ==============================================================================
# 6.4 COMPARATIVA RÁPIDA MLP
# ==============================================================================
comparativa_mlp <- rbind(
  MLP_Baseline_SoloExog = metrics_mlp_base,
  MLP_Opt_Lags_Exog     = metrics_mlp_opt)

print(comparativa_mlp)

##############################################################################
##                            MODELO 7: NNAR                                ##
##############################################################################

library(forecast)
library(dplyr)

# ----------------------------------------------------------------------------
# 7.0 PREPARACIÓN DE DATOS (SINCRONIZACIÓN DE ÍNDICES)
# ----------------------------------------------------------------------------
nnar_y_train     <- window(y_ts, start = c(1996, 1), end = c(2021, 12))
nnar_y_val       <- window(y_ts, start = c(2022, 1), end = c(2023, 12))
nnar_y_train_val <- window(y_ts, start = c(1996, 1), end = c(2023, 12))
nnar_y_test      <- window(y_ts, start = c(2024, 1), end = c(2025, 12))

# Matrices exógenas específicas para NNAR
nnar_xreg_train     <- as.matrix(window(x_ts, start = c(1996, 1), end = c(2021, 12)))
nnar_xreg_val       <- as.matrix(window(x_ts, start = c(2022, 1), end = c(2023, 12)))
nnar_xreg_train_val <- as.matrix(window(x_ts, start = c(1996, 1), end = c(2023, 12)))
nnar_xreg_test      <- as.matrix(window(x_ts, start = c(2024, 1), end = c(2025, 12)))

# Vectores numéricos limpios para la evaluación de métricas
nnar_y_val_vec  <- as.numeric(nnar_y_val)
nnar_y_test_vec <- as.numeric(nnar_y_test)

nnar_h_val  <- length(nnar_y_val_vec)
nnar_h_test <- length(nnar_y_test_vec)


# ----------------------------------------------------------------------------
# 7.1 MODELO BASELINE (NNAR UNIVARIANTE SIN OPTIMIZAR)
# ----------------------------------------------------------------------------
set.seed(123)
fit_nnar_base <- nnetar(
  nnar_y_train_val, 
  p       = 12, 
  size    = 10, 
  repeats = 20)

pred_nnar_base <- forecast(fit_nnar_base, h = nnar_h_test)

metrics_nnar_base <- calc_metrics(real = nnar_y_test_vec, pred = as.numeric(pred_nnar_base$mean))

print(metrics_nnar_base)


# ----------------------------------------------------------------------------
# 7.2 BÚSQUEDA EN MALLA SOBRE VALIDACIÓN (2022-2023) - NNAR MULTIVARIANTE
# ----------------------------------------------------------------------------
grid_nnar <- expand.grid(
  size    = c(4, 6, 8, 10, 12),
  repeats = c(10, 20, 30) )

resultados_malla_nnar <- data.frame()

for (i in 1:nrow(grid_nnar)) {
  s <- grid_nnar$size[i]
  r <- grid_nnar$repeats[i]
  
  set.seed(123)
  fit_temp <- nnetar(
    nnar_y_train,
    p       = 12,
    size    = s,
    repeats = r,
    xreg    = nnar_xreg_train )
  
  # Predicción en periodo de Validación con exógenas de validación
  pred_val_obj <- forecast(fit_temp, h = nnar_h_val, xreg = nnar_xreg_val)
  pred_val     <- as.numeric(pred_val_obj$mean)
  
  # Evaluación con la función centralizada de métricas
  m_val <- calc_metrics(real = nnar_y_val_vec, pred = pred_val)
  
  res_row <- data.frame(
    size    = s,
    repeats = r,
    MAE     = m_val$MAE,
    RMSE    = m_val$RMSE,
    MAPE    = m_val$MAPE )
  
  resultados_malla_nnar <- rbind(resultados_malla_nnar, res_row)
}

# Selección del mejor conjunto de hiperparámetros por RMSE
best_params_nnar <- resultados_malla_nnar %>% arrange(RMSE) %>% slice(1)
print(best_params_nnar)


# ----------------------------------------------------------------------------
# 7.3 REENTRENAMIENTO OPTIMIZADO Y EVALUACIÓN EN TEST (2024-2025)
# ----------------------------------------------------------------------------
set.seed(123)
fit_nnar_opt <- nnetar(
  nnar_y_train_val,
  p       = 12,
  size    = best_params_nnar$size,
  repeats = best_params_nnar$repeats,
  xreg    = nnar_xreg_train_val)

pred_opt_obj  <- forecast(fit_nnar_opt, h = nnar_h_test, xreg = nnar_xreg_test)
pred_opt_test <- as.numeric(pred_opt_obj$mean)

metrics_nnar_opt <- calc_metrics(real = nnar_y_test_vec, pred = pred_opt_test)
print(metrics_nnar_opt)

# ----------------------------------------------------------------------------
# 7.4 COMPARATIVA RÁPIDA NNAR
# ----------------------------------------------------------------------------
comparativa_nnar <- rbind(
  NNAR_Univariante_Base = metrics_nnar_base,
  NNAR_Multivariante_Opt = metrics_nnar_opt)
print(comparativa_nnar)

##############################################################################
##                 MODELO 8: MODELO HÍBRIDO CEEMDAN + XGBOOST               ##
##############################################################################

library(Rlibeemd)
library(forecast)
library(xgboost)
library(dplyr)

# ----------------------------------------------------------------------------
# 8.0 DESCOMPOSICIÓN CEEMDAN SOBRE TRAIN + VAL (1996-2023)
# ----------------------------------------------------------------------------

h_test  <- length(y_test)
n_train <- length(y_train_val)

set.seed(123)
matriz_ceemdan <- ceemdan(
  y_train_val, 
  ensemble_size  = 250, 
  noise_strength = 0.2, 
  S_number       = 4)

num_componentes <- ncol(matriz_ceemdan)
num_componentes

# Convertir columnas en lista de series temporales
imfs_list <- lapply(1:num_componentes, function(i) {
  ts(matriz_ceemdan[, i], start = start(y_train_val), frequency = frequency(y_train_val))
})

names(imfs_list) <- c(paste0("IMF_", 1:(num_componentes - 1)), "Residuo")

# ----------------------------------------------------------------------------
# 8.1 ENTRENAMIENTO Y PREDICCIÓN POR COMPONENTE CON XGBOOST
# ----------------------------------------------------------------------------

# Creamos la lista para almacenar los modelos sin alterar tu lógica
modelos_ceemdan_xgb <- list()
X_train_imf_last    <- NULL

pred_imfs_xgb <- sapply(names(imfs_list), function(imf_name) {
  imf_ts <- imfs_list[[imf_name]]
  
  fit_imf_ar   <- auto.arima(imf_ts)
  pred_imf_ext <- forecast(fit_imf_ar, h = h_test)$mean
  
  imf_completa <- c(as.numeric(imf_ts), as.numeric(pred_imf_ext))
  
  mat_lags_completa   <- crear_lags_matriz(imf_completa, num_lags = 12)
  x_exogenas_completa <- rbind(x_train_val, x_test)
  
  X_completo <- cbind(
    Mes = as.numeric(cycle(ts(1:length(imf_completa), start = start(y_train_val), frequency = 12))),
    mat_lags_completa,
    x_exogenas_completa)
  
  y_completo_imf <- c(as.numeric(imf_ts), rep(NA, h_test))
  idx_validos    <- 13:length(imf_completa)
  
  X_clean <- X_completo[idx_validos, , drop = FALSE]
  y_clean <- y_completo_imf[idx_validos]
  
  n_train_val_clean <- length(imf_ts) - 12
  
  X_train_imf <- X_clean[1:n_train_val_clean, , drop = FALSE]
  y_train_imf <- y_clean[1:n_train_val_clean]
  X_test_imf  <- X_clean[(n_train_val_clean + 1):nrow(X_clean), , drop = FALSE]
  
  dtrain_imf <- xgb.DMatrix(data = apply(X_train_imf, 2, as.numeric), label = as.numeric(y_train_imf))
  dtest_imf  <- xgb.DMatrix(data = apply(X_test_imf, 2, as.numeric))
  
  set.seed(123)
  fit_xgb_imf <- xgb.train(
    params  = params_base, 
    data    = dtrain_imf, 
    nrounds = 150, 
    verbose = 0)
  
  modelos_ceemdan_xgb[[imf_name]] <<- fit_xgb_imf
  X_train_imf_last                <<- X_train_imf
  
  return(predict(fit_xgb_imf, dtest_imf))
})

# Reconstrucción: suma de las predicciones de XGBoost sobre las IMFs
pred_ceemdan_xgb <- rowSums(pred_imfs_xgb)

# ----------------------------------------------------------------------------
# 8.2 EVALUACIÓN Y RESULTADOS CEEMDAN + XGBOOST
# ----------------------------------------------------------------------------
metrics_ceemdan_xgb <- calc_metrics(
  real = as.numeric(y_test), 
  pred = pred_ceemdan_xgb)

print(metrics_ceemdan_xgb)

##############################################################################
##         PLOTS INDIVIDUALES INDEPENDIENTES EN TEST (1996-2025)            ##
##############################################################################

library(ggplot2)
library(dplyr)

# ----------------------------------------------------------------------------
# 1. PREPARACIÓN DE FECHAS Y DATOS REALES (360 MESES)
# ----------------------------------------------------------------------------
fechas_completas <- seq.Date(from = as.Date("1996-01-01"), to = as.Date("2025-12-01"), by = "month")
fechas_test      <- seq.Date(from = as.Date("2024-01-01"), to = as.Date("2025-12-01"), by = "month")

y_real_tv   <- as.numeric(y_train_val)
y_real_test <- as.numeric(y_test)

df_base_plot <- data.frame(
  Fecha = fechas_completas,
  Real  = c(y_real_tv, y_real_test),
  Tipo  = c(rep("Histórico (1996-2023)", length(y_real_tv)), 
            rep("Real Test (2024-2025)", length(y_real_test)))
)

# ----------------------------------------------------------------------------
# 2. FUNCIÓN AUXILIAR PARA EXTRAER VECTORES NUMÉRICOS DE CUALQUIER TIPO DE OBJETO
# ----------------------------------------------------------------------------
extraer_vector_num <- function(obj) {
  if (is.null(obj)) return(rep(NA, length(fechas_test)))
  
  if (is.numeric(obj)) {
    return(as.numeric(obj))
  } else if ("mean" %in% names(obj)) {
    return(as.numeric(obj$mean))       # Para objetos de la librería forecast
  } else if ("pred" %in% names(obj)) {
    return(as.numeric(obj$pred))       # Para predict.Arima / predict.arima
  } else if ("yhat" %in% names(obj)) {
    return(as.numeric(obj$yhat))       # Para data frames de Prophet
  } else {
    return(as.numeric(obj))
  }
}

# ----------------------------------------------------------------------------
# 3.  PREDICCIONES DE LOS 8 MODELOS
# ----------------------------------------------------------------------------
predicciones_modelos <- list(
  "Modelo 1: Holt-Winters"       = extraer_vector_num(pred_hw_test),
  "Modelo 2: SARIMA"             = extraer_vector_num(pred_sarima_test),
  "Modelo 3: SARIMAX"            = extraer_vector_num(pred_sarimax_test),
  "Modelo 4: Prophet"            = extraer_vector_num(pred_prophet_exo_tuned),
  "Modelo 5: XGBoost"            = extraer_vector_num(pred_xgb_tuned),
  "Modelo 6: MLP Optimizado"     = extraer_vector_num(mlp_pred_opt_test),
  "Modelo 7: NNAR Multivariante" = extraer_vector_num(pred_opt_obj),
  "Modelo 8: CEEMDAN + XGBoost"  = extraer_vector_num(pred_ceemdan_xgb)
)

colores_tfm <- c(
  "Histórico (1996-2023)" = "#832D51",
  "Real Test (2024-2025)" = "#447A5F",
  "Predicción Modelo"     = "#A6C048"
)

# ----------------------------------------------------------------------------
# 4. FUNCIÓN PARA GENERAR UN GRÁFICO INDIVIDUAL
# ----------------------------------------------------------------------------
generar_plot_unico <- function(nombre_modelo, vector_pred) {
  
  df_pred <- data.frame(
    Fecha      = fechas_test,
    Prediccion = vector_pred
  )
  
  ggplot() +
    geom_line(data = filter(df_base_plot, Tipo == "Histórico (1996-2023)"), 
              aes(x = Fecha, y = Real, color = "Histórico (1996-2023)"), 
              linewidth = 0.8) +
    geom_line(data = filter(df_base_plot, Tipo == "Real Test (2024-2025)"), 
              aes(x = Fecha, y = Real, color = "Real Test (2024-2025)"), 
              linewidth = 0.8) +
    geom_line(data = df_pred, 
              aes(x = Fecha, y = Prediccion, color = "Predicción Modelo"), 
              linewidth = 0.8) +
    scale_color_manual(values = colores_tfm) +
    labs(
      title = nombre_modelo,
      x = "Año", 
      y = "Viajeros"
    ) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    theme_minimal(base_size = 11) + 
    theme(
      plot.title        = element_text(face = "bold", size = 13, color = "Black", hjust = 0),
      legend.position   = "bottom",
      legend.title      = element_blank(),
      legend.text       = element_text(size = 10),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(color = "black", size = 9),
      axis.title        = element_text(face = "bold", size = 10)
    )
}

# ----------------------------------------------------------------------------
# 5. GENERACIÓN Y GUARDADO DE FICHEROS PNG INDIVIDUALES (300 DPI)
# ----------------------------------------------------------------------------
for (i in seq_along(predicciones_modelos)) {
  
  nombre   <- names(predicciones_modelos)[i]
  pred_vec <- predicciones_modelos[[nombre]]
  
  p <- generar_plot_unico(nombre, pred_vec)
  
  print(p)
  
  nombre_archivo <- paste0("Figura_Modelo_", i, ".png")
  ggsave(
    filename = nombre_archivo,
    plot     = p,
    width    = 16,
    height   = 7,
    dpi      = 300
  )
}


#=======================================================================
#-------------VALIDACIÓN CRUZADA EXPANSIVA COMPLETA (8 MODELOS)---------
#=======================================================================

library(forecast)
library(dplyr)
library(ggplot2)
library(prophet)
library(xgboost)
library(nnet)

evaluar_cv_modelo_fijo <- function(y_ts, n_train, paso_h = 12, predict_fun, nombre_modelo = "Modelo") {
  n_total <- length(y_ts)
  num_iteraciones <- n_total - paso_h - n_train + 1
  
  mat_errores <- matrix(NA, nrow = num_iteraciones, ncol = paso_h)
  mat_reales  <- matrix(NA, nrow = num_iteraciones, ncol = paso_h)
  
  idx <- 1
  for (i in n_train:(n_total - paso_h)) {
    datos_hasta_hoy <- window(y_ts, end = time(y_ts)[i])
    preds_h         <- predict_fun(datos_hasta_hoy, h = paso_h, pos_i = i)
    reales_h        <- window(y_ts, start = time(y_ts)[i + 1], end = time(y_ts)[i + paso_h])
    
    mat_errores[idx, ] <- as.numeric(reales_h) - as.numeric(preds_h)
    mat_reales[idx, ]  <- as.numeric(reales_h)
    
    idx <- idx + 1
  }
  
  # A. Métricas Promedio Globales (promediando todos los pasos h=1 a h=12)
  mae_prom  <- mean(abs(mat_errores), na.rm = TRUE)
  rmse_prom <- sqrt(mean(mat_errores^2, na.rm = TRUE))
  mape_prom <- mean(abs(mat_errores / mat_reales), na.rm = TRUE) * 100
  
  res_resumen <- data.frame(
    Modelo    = nombre_modelo,
    MAE_prom  = round(mae_prom, 2),
    RMSE_prom = round(rmse_prom, 2),
    MAPE_prom = round(mape_prom, 2)
  )
  
  # B. Curva completa de MAPE (h=1 a h=12) -> Se mantiene intacta para el gráfico
  mape_por_h <- colMeans(abs(mat_errores / mat_reales), na.rm = TRUE) * 100
  
  df_curva   <- data.frame(
    Modelo = nombre_modelo,
    h      = 1:paso_h, 
    MAPE   = round(mape_por_h, 2)
  )
  
  return(list(resumen = res_resumen, curva = df_curva))
}

# ==============================================================================
# 2. HOLT-WINTERS 
# ==============================================================================
pred_hw_fijo <- function(datos_hasta_hoy, h, pos_i) {
  pred_obj <- suppressWarnings(
    tryCatch({
      fit_actual <- HoltWinters(datos_hasta_hoy, optim.control = list(maxit = 1000))
      predict(fit_actual, n.ahead = h)
    }, error = function(e) {
      fit_actual <- HoltWinters(datos_hasta_hoy, gamma = FALSE) 
      predict(fit_actual, n.ahead = h)
    })
  )
  
  return(as.numeric(pred_obj))
}

cv_hw <- evaluar_cv_modelo_fijo(
  y_ts = y_ts, n_train = length(y_train_val), paso_h = 12, 
  predict_fun = pred_hw_fijo, nombre_modelo = "Holt-Winters"
)

# ==============================================================================
# 3. SARIMA 
# ==============================================================================
pred_sarima_fijo <- function(datos_hasta_hoy, h, pos_i) {
  # Refit con el objeto fit_sarima previamente ajustado
  fit_actual <- Arima(datos_hasta_hoy, model = fit_sarima)
  pred_obj   <- forecast(fit_actual, h = h)
  return(as.numeric(pred_obj$mean))
}

cv_sarima <- evaluar_cv_modelo_fijo(
  y_ts = y_ts, n_train = length(y_train_val), paso_h = 12, 
  predict_fun = pred_sarima_fijo, nombre_modelo = "SARIMA"
)

# ==============================================================================
# 4. SARIMAX 
# ==============================================================================
# 1. Creamos la matriz filtrada con las 2 exógenas finales (cambia los nombres por las tuyas)
x_mat_sarimax <- x_mat[, c("FestivosN", "COVID"), drop = FALSE]

# 2. Adaptador de SARIMAX corregido usando la matriz de 2 variables
pred_sarimax_fijo <- function(datos_hasta_hoy, h, pos_i) {
  # Extraemos del subconjunto de 2 variables
  xreg_train_i <- x_mat_sarimax[1:pos_i, , drop = FALSE]
  xreg_test_i  <- x_mat_sarimax[(pos_i + 1):(pos_i + h), , drop = FALSE]
  
  fit_actual <- Arima(
    datos_hasta_hoy, 
    model = fit_sarimax, 
    xreg  = xreg_train_i
  )
  
  pred_obj <- forecast(fit_actual, h = h, xreg = xreg_test_i)
  return(as.numeric(pred_obj$mean))
}

cv_sarimax <- evaluar_cv_modelo_fijo(
  y_ts = y_ts, n_train = length(y_train_val), paso_h = 12, 
  predict_fun = pred_sarimax_fijo, nombre_modelo = "SARIMAX"
)

# ==============================================================================
# 5. PROPHET
# ==============================================================================
pred_prophet_fijo <- function(datos_hasta_hoy, h, pos_i) {
  # 1. Obtener el año y mes inicial de la ventana actual
  st <- start(datos_hasta_hoy)
  fecha_inicio <- as.Date(paste(st[1], st[2], "01", sep = "-"))
  
  # 2. Generar el vector de fechas mensual correcto para Prophet
  fechas <- seq(from = fecha_inicio, length.out = length(datos_hasta_hoy), by = "month")
  
  df_p <- data.frame(
    ds = fechas,
    y  = as.numeric(datos_hasta_hoy)
  )
  
  # 3. Ajustar Prophet de forma silenciosa
  m <- prophet(df_p, yearly.seasonality = TRUE,
               weekly.seasonality = FALSE, 
               daily.seasonality = FALSE,
               changepoint.prior.scale = best_params$changepoint_prior_scale,
               seasonality.prior.scale = best_params$seasonality_prior_scale,
               seasonality.mode        = "additive",
               verbose = FALSE)
  futuro <- make_future_dataframe(m, periods = h, freq = "month")
  forecast_p <- predict(m, futuro)
  
  return(tail(forecast_p$yhat, h))
}

cv_prophet <- evaluar_cv_modelo_fijo(
  y_ts = y_ts, n_train = length(y_train_val), paso_h = 12, 
  predict_fun = pred_prophet_fijo, nombre_modelo = "Prophet Tuned"
)

# ==============================================================================
# 6. XGBOOST TUNED
# ==============================================================================
pred_xgb_fijo <- function(datos_hasta_hoy, h, pos_i) {
  i <- length(datos_hasta_hoy)
  
  y_ext_ts <- ts(
    c(as.numeric(datos_hasta_hoy), rep(NA, h)), 
    start     = start(datos_hasta_hoy), 
    frequency = frequency(datos_hasta_hoy)
  )
  
  df_xgb_i <- preparar_datos_xgboost(
    y_ts     = y_ext_ts, 
    xreg_mat = x_mat[1:(i + h), , drop = FALSE], 
    num_lags = 12
  )
  
  features_cols <- setdiff(colnames(df_xgb_i), c("Fecha", "y"))
  
  df_train_i <- df_xgb_i[1:(i - 12), ]
  df_test_i  <- tail(df_xgb_i, h)
  
  mat_train <- matrix(as.numeric(as.matrix(df_train_i[, features_cols])), nrow = nrow(df_train_i))
  mat_test  <- matrix(as.numeric(as.matrix(df_test_i[, features_cols])),  nrow = nrow(df_test_i))
  
  dtrain <- xgb.DMatrix(data = mat_train, label = as.numeric(df_train_i$y))
  dtest  <- xgb.DMatrix(data = mat_test)
  
  num_rounds <- if (exists("params_optimos") && !is.null(params_optimos$nrounds)) {
    params_optimos$nrounds
  } else if (exists("nrounds_opt")) {
    nrounds_opt
  } else if (exists("best_nrounds")) {
    best_nrounds
  } else {
    100
  }
  

  set.seed(123)
  fit_xgb <- xgb.train(
    params  = params_optimos, 
    data    = dtrain, 
    nrounds = best_xgb_params$nrounds,
    verbose = 0
  )
  
  return(predict(fit_xgb, dtest))
}

cv_xgb <- evaluar_cv_modelo_fijo(
  y_ts = y_ts, n_train = length(y_train_val), paso_h = 12, 
  predict_fun = pred_xgb_fijo, nombre_modelo = "XGBoost"
)

# ==============================================================================
# 7. NNAR MULTIVARIANTE
# ==============================================================================
pred_nnar_exo_fijo <- function(datos_hasta_hoy, h, pos_i) {
  i <- length(datos_hasta_hoy)
  
  xreg_train_i <- x_mat[1:i, , drop = FALSE]
  xreg_test_i  <- x_mat[(i + 1):(i + h), , drop = FALSE]
  
  fit_actual <- nnetar(
    datos_hasta_hoy, 
    p       = 12, 
    size    = 4, 
    repeats = 30, 
    xreg    = xreg_train_i
  )
  
  pred_obj <- forecast(fit_actual, h = h, xreg = xreg_test_i)
  return(as.numeric(pred_obj$mean))
}

cv_nnar <- evaluar_cv_modelo_fijo(
  y_ts = y_ts, n_train = length(y_train_val), paso_h = 12, 
  predict_fun = pred_nnar_exo_fijo, nombre_modelo = "NNAR + Exógenas"
)

# ==============================================================================
# 8. MLP
# ==============================================================================
pred_mlp_fijo <- function(datos_hasta_hoy, h, pos_i) {
  i <- length(datos_hasta_hoy)
  
  xreg_train_i <- x_mat[1:i, , drop = FALSE]
  xreg_test_i  <- x_mat[(i + 1):(i + h), , drop = FALSE]
  
  # Red feedforward MLP sin componente autorregresivo estacional explícito (P=0)
  fit_actual <- nnetar(
    datos_hasta_hoy, 
    p       = 12, 
    P       = 0, 
    size    = 8, 
    decay   = 0.01,
    xreg    = xreg_train_i
  )
  
  pred_obj <- forecast(fit_actual, h = h, xreg = xreg_test_i)
  return(as.numeric(pred_obj$mean))
}

cv_mlp <- evaluar_cv_modelo_fijo(
  y_ts = y_ts, n_train = length(y_train_val), paso_h = 12, 
  predict_fun = pred_mlp_fijo, nombre_modelo = "MPL"
)

# ==============================================================================
# 9. CEEMDAN + XGBOOST 
# ==============================================================================
pred_ceemdan_xgb_fijo <- function(datos_hasta_hoy, h, pos_i) {
  
  # 1. Descomposición CEEMDAN devuelve una matriz donde las columnas son las IMFs + Residuo
  decomp_matrix <- ceemdan(as.numeric(datos_hasta_hoy), num_siftings = 50)
  
  preds_imf_sum <- numeric(h)
  
  # 2. Predecir cada componente (IMFs) de forma independiente con XGBoost
  for (j in 1:ncol(decomp_matrix)) {
    imf_j <- ts(
      decomp_matrix[, j], 
      start     = start(datos_hasta_hoy), 
      frequency = frequency(datos_hasta_hoy)
    )
    
    # Predecimos la componente j con el adaptador XGBoost
    pred_j <- pred_xgb_fijo(imf_j, h = h, pos_i = pos_i)
    
    # Sumamos la predicción al total agregado
    preds_imf_sum <- preds_imf_sum + pred_j
  }
  
  return(preds_imf_sum)
}

# Ejecutar validación cruzada para CEEMDAN + XGBoost
cv_ceemdan <- evaluar_cv_modelo_fijo(
  y_ts        = y_ts, 
  n_train     = length(y_train_val), 
  paso_h      = 12, 
  predict_fun = pred_ceemdan_xgb_fijo, 
  nombre_modelo = "CEEMDAN + XGBoost"
)
# ==============================================================================
# 10. TABLA RESUMEN GLOBAL (8 MODELOS)
# ==============================================================================
tabla_cv_global_completa <- bind_rows(
  cv_hw$resumen,
  cv_sarima$resumen,
  cv_sarimax$resumen,
  cv_prophet$resumen,
  cv_xgb$resumen,
  cv_nnar$resumen,
  cv_mlp$resumen,
  cv_ceemdan$resumen
) 

cat("\n=== RESUMEN GLOBAL VALIDACIÓN CRUZADA (MÉTRICAS PROMEDIO) ===\n")
print(tabla_cv_global_completa, row.names = FALSE)

# ==============================================================================
# 11. GRÁFICA COMPARATIVA DE ERRORES (MAPE)
# ==============================================================================
df_curvas_global_completa <- bind_rows(
  cv_hw$curva,
  cv_sarima$curva,
  cv_sarimax$curva,
  cv_prophet$curva,
  cv_xgb$curva,
  cv_nnar$curva,
  cv_mlp$curva,
  cv_ceemdan$curva
)

ggplot(df_curvas_global_completa, aes(x = h, y = MAPE, color = Modelo, group = Modelo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:12) +
  theme_minimal() +
  labs(
    title = "Evolución del MAPE según el horizonte de predicción (h = 1 a 12)",
    x = "Horizonte de Predicción (h)",
    y = "MAPE (%)",
    color = "Modelo"
  ) +

  theme(
 
    plot.title = element_text(size = 24, , face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 12))
    
# ==============================================================================
# ANÁLISIS DE EXPLICABILIDAD SHAP CON FASTSHAP (5 MODELOS)
# ==============================================================================

library(kernelshap)
library(shapviz)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. SHAP PARA PROPHET TUNED
# ------------------------------------------------------------------------------
X_eval_prophet   <- tail(pro_df_train_val_exo, 60)
df_preds_prophet <- predict(fit_prophet_exo_tuned, X_eval_prophet)

# 1. Matriz de contribuciones aditivas (S) y de valores reales (X)
shap_matrix_prophet    <- as.matrix(df_preds_prophet[, nombres_exogenas])
feature_matrix_prophet <- as.matrix(X_eval_prophet[, nombres_exogenas])

# 2. Construcción del objeto shapviz indicando 'object = shap_matrix_prophet'
sv_prophet <- shapviz(
  object   = shap_matrix_prophet, 
  X        = feature_matrix_prophet,
  baseline = mean(df_preds_prophet$trend)
)

# 3. Gráfica Beeswarm
p_shap_prophet <- sv_importance(sv_prophet, kind = "beeswarm") +
  scale_color_gradient(low = "#A6C048", high = "#832D51") +
  theme_minimal() +
  labs(title = "Explicabilidad SHAP - Prophet Tuned") +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p_shap_prophet)
ggsave("SHAP_Prophet.png", plot = p_shap_prophet, width = 8, height = 4, dpi = 300)

# ------------------------------------------------------------------------------
# 2. SHAP PARA XGBOOST TUNED
# ------------------------------------------------------------------------------
sv_xgb <- shapviz(fit_xgb_tuned, X_pred = xgb_X_train_val)

p_shap_xgb <- sv_importance(sv_xgb, kind = "beeswarm") +
  scale_color_gradient(low = "#A6C048", high = "#832D51") +
  theme_minimal() +
  labs(
    title = "Explicabilidad SHAP - XGBoost Tuned",
    x = "Valor SHAP (Impacto en la predicción)",
    y = "Característica"
  ) +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p_shap_xgb)
ggsave("SHAP_XGBoost.png", plot = p_shap_xgb, width = 8, height = 4, dpi = 300)

# ------------------------------------------------------------------------------
# 3. SHAP PARA MLP (nnet)
# ------------------------------------------------------------------------------
pred_wrapper_mlp <- function(object, newdata, ...) {
  as.numeric(predict(object, newdata = as.data.frame(newdata)))
}

X_shap_mlp <- mlp_df_tr_val %>% dplyr::select(-ViajerosAVyLD)

set.seed(123)
bg_mlp <- dplyr::sample_n(X_shap_mlp, size = 20)
X_eval_mlp <- tail(X_shap_mlp, 60) # Evaluamos sobre los últimos 5 años

ks_mlp <- kernelshap(
  object   = mlp_model_opt,
  X        = X_eval_mlp,
  pred_fun = pred_wrapper_mlp,
  bg_X     = bg_mlp
)

sv_mlp <- shapviz(ks_mlp)
p_shap_mlp <- sv_importance(sv_mlp, kind = "beeswarm") +
  scale_color_gradient(low = "#A6C048", high = "#832D51") +
  theme_minimal() +
  labs(
    title = "Explicabilidad SHAP - MLP (nnet)",
    x = "Valor SHAP (Impacto en la predicción)",
    y = "Característica"
  ) +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p_shap_mlp)
ggsave("SHAP_MPL.png", plot = p_shap_mlp, width = 8, height = 4, dpi = 300)
# ------------------------------------------------------------------------------
# 4. SHAP PARA NNAR OPTIMIZADO
# ------------------------------------------------------------------------------
pred_wrapper_nnar <- function(object, newdata, ...) {
  xreg_m <- as.matrix(newdata)
  as.numeric(forecast(object, h = 1, xreg = xreg_m)$mean)
}

X_shap_nnar <- as.data.frame(nnar_xreg_train_val)

set.seed(123)
bg_nnar <- dplyr::sample_n(X_shap_nnar, size = 20)
X_eval_nnar <- tail(X_shap_nnar, 60)

ks_nnar <- kernelshap(
  object   = fit_nnar_opt,
  X        = X_eval_nnar,
  pred_fun = pred_wrapper_nnar,
  bg_X     = bg_nnar
)

sv_nnar <- shapviz(ks_nnar)
p_shap_nnar <- sv_importance(sv_nnar, kind = "beeswarm") +
  scale_color_gradient(low = "#A6C048", high = "#832D51") +
  theme_minimal() +
  labs(
    title = "Explicabilidad SHAP - NNAR",
    x = "Valor SHAP (Impacto en la predicción)",
    y = "Característica"
  ) +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p_shap_nnar)
ggsave("SHAP_NNAR.png", plot = p_shap_nnar, width = 8, height = 4, dpi = 300)

# ------------------------------------------------------------------------------
# 5. SHAP PARA CEEMDAN + XGBOOST
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RECONSTRUCCIÓN DE LA MATRIZ DE ENTRADA DE CEEMDAN
# ------------------------------------------------------------------------------
# 1. Reconstruir lags sobre la serie histórica de IMFs
imf_1_ts <- imfs_list[[1]]
n_train_val_raw <- length(imf_1_ts)

# Crear matriz de lags para la primera IMF como referencia de dimensiones
mat_lags_ceemdan <- crear_lags_matriz(as.numeric(imf_1_ts), num_lags = 12)

# Unir con la variable Mes y las exógenas de entrenamiento
X_completo_ceemdan <- cbind(
  Mes = as.numeric(cycle(ts(1:n_train_val_raw, start = start(y_train_val), frequency = 12))),
  mat_lags_ceemdan,
  x_train_val
)

# Filtrar las observaciones válidas 
idx_validos <- 13:n_train_val_raw
X_shap_ceemdan <- as.data.frame(X_completo_ceemdan[idx_validos, , drop = FALSE])

# ------------------------------------------------------------------------------
# SHAP CEEMDAN + XGBoost
# ------------------------------------------------------------------------------

pred_wrapper_ceemdan_xgb <- function(object_list, newdata, ...) {
  dmat <- xgb.DMatrix(data = apply(as.matrix(newdata), 2, as.numeric))
  preds_totales <- numeric(nrow(newdata))
  
  for (m in object_list) {
    preds_totales <- preds_totales + predict(m, dmat)
  }
  return(preds_totales)
}

X_shap_ceemdan <- as.data.frame(X_train_imf_last)

set.seed(123)
bg_ceemdan     <- dplyr::sample_n(X_shap_ceemdan, size = 20)
X_eval_ceemdan <- tail(X_shap_ceemdan, 60)

ks_ceemdan <- kernelshap(
  object   = modelos_ceemdan_xgb,
  X        = X_eval_ceemdan,
  pred_fun = pred_wrapper_ceemdan_xgb,
  bg_X     = bg_ceemdan
)

sv_ceemdan <- shapviz(ks_ceemdan)


# Gráfica con tus colores personalizados (#447A5F y #832D51)
p_shap_ceemdan <- sv_importance(sv_ceemdan, kind = "beeswarm") +
  scale_color_gradient(low = "#A6C048", high = "#832D51") +
  theme_minimal() +
  labs(
    title = "Explicabilidad SHAP - CEEMDAN + XGBoost",
    x = "Valor SHAP (Impacto en la predicción)",
    y = "Característica",
    color = "Valor de la\nvariable"
  ) +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p_shap_ceemdan)
ggsave("SHAP_CEEMDAN_XGBoost.png", plot = p_shap_ceemdan, width = 8, height = 4, dpi = 300)


# Guardar los objetos shapviz calculados para no repetir el bucle pesado
saveRDS(sv_prophet, "sv_prophet.rds")
saveRDS(sv_xgb, "sv_xgb.rds")
saveRDS(sv_mlp, "sv_mlp.rds")
saveRDS(sv_nnar, "sv_nnar.rds")
saveRDS(sv_ceemdan, "sv_ceemdan.rds")

sv_prophet <- readRDS("sv_prophet.rds")
sv_xgb     <- readRDS("sv_xgb.rds")
sv_mlp     <- readRDS("sv_mlp.rds")
sv_nnar    <- readRDS("sv_nnar.rds")
sv_ceemdan <- readRDS("sv_ceemdan.rds")
