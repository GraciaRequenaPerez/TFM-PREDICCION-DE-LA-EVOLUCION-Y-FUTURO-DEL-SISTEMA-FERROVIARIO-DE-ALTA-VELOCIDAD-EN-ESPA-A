library(readxl)
library(forecast)
library(lmtest)
library(ggplot2)
PREDFIN <- read_excel("C:/Users/graci/Desktop/TFM/Datos/PREDFIN.xlsx")

PREDFIN_train <- subset(PREDFIN, !is.na(ViajerosAVyLD))
PREDFIN_pred  <- subset(PREDFIN, is.na(ViajerosAVyLD))

# 3. Serie temporal histórica para el modelo
ts_viajeros_FIN <- ts(PREDFIN_train$ViajerosAVyLD, start = c(1996, 1), frequency = 12)
plot(ts_viajeros_FIN, col = "#832D51")


# 4. Matriz exógena de entrenamiento
xreg_train_FIN <- cbind(
  FestivosN  = PREDFIN_train$FestivosN,
  COVID      = PREDFIN_train$COVID,
  Accidentes = PREDFIN_train$Accidentes
)

# 5. Reestimación del modelo SARIMAX sobre el train completo
fit_shock <- Arima(
  ts_viajeros_FIN,
  order = c(2, 1, 2),
  seasonal = list(order = c(0, 1, 2), period = 12),
  xreg = xreg_train_FIN,
  method = "ML"
)

AIC(fit_shock)

# 6. Test de significación de los coeficientes
coeftest(fit_shock)

# 1. Matriz exógena futura (julio 2026 - julio 2031)
xreg_pred <- cbind(
  FestivosN  = PREDFIN_pred$FestivosN,
  COVID      = PREDFIN_pred$COVID,      # Fijado a 0
  Accidentes = PREDFIN_pred$Accidentes # Según el escenario configurado en el Excel
)

# 2. Generar la predicción fuera de muestra
forecast_5y <- forecast(fit_shock, xreg = xreg_pred, h = nrow(PREDFIN_pred))

# 3. Representación gráfica de la proyección
autoplot(forecast_5y) +
  labs(
    title = "Proyección de número de viajeros en Alta Velocidad y Larga Distancia a 5 Años: Modelo SARIMAX",
    subtitle = "Horizonte 2026-2031 ",
    x = "Año",
    y = "Viajeros (en miles)"
  ) +
  theme_minimal()

library(ggplot2)
library(forecast)
library(scales)

# Versión monocromática robusta en tono VINO (#832D51)
autoplot(forecast_5y, fcol = "#62223E") +
  
  # Serie histórica en el mismo color vino
  geom_line(color = "#62223E", linewidth = 0.6) +
  
  # Formato numérico en español (puntos para miles: 2.500, 5.000...)
  scale_y_continuous(labels = label_number(big.mark = ".", decimal.mark = ",")) +
  
  # Títulos y etiquetas
  labs(
    title = "Proyección de número de viajeros en alta Velocidad y larga distancia a 5 Años",
    subtitle = "Modelo SARIMAX (2026-2031) | Intervalos de confianza del 80% y 95%",
    x = "Año",
    y = "Viajeros (en miles)",
  ) +
  
  # Tema minimalista con textos en NEGRO y tamaños editables inline
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 13, color = "black", margin = margin(b = 6)),
    plot.subtitle = element_text(size = 10, color = "#333333", margin = margin(b = 14)),
    
    axis.title.x  = element_text(face = "bold", size = 11, color = "black", margin = margin(t = 8)),
    axis.title.y  = element_text(face = "bold", size = 11, color = "black", margin = margin(r = 8)),
    axis.text     = element_text(size = 9, color = "black"),
    
    panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    plot.margin   = margin(15, 15, 15, 15)
  )

library(dplyr)
library(lubridate)
library(knitr)

# 1. Crear el vector de fechas para el horizonte proyectado (agosto 2026 - julio 2031)
fechas_futuras <- seq(from = as.Date("2026-08-01"), by = "month", length.out = length(forecast_5y$mean))

# 2. Construir el dataframe mensual base en miles de viajeros
df_mensual <- data.frame(
  Fecha          = fechas_futuras,
  Año            = year(fechas_futuras),
  Viajeros_Miles = as.numeric(forecast_5y$mean)
)

# 3. Agrupar por año, convertir a MILLONES de viajeros y aplicar escenarios (+/- 3%)
tabla_escenarios_millones <- df_mensual %>%
  group_by(Año) %>%
  summarise(
    Meses         = n(),
    # Conversión directa de miles a millones de viajeros (dividiendo entre 1.000)
    Base_Millones = sum(Viajeros_Miles) / 1000
  ) %>%
  mutate(
    # i representa el incremento anual acumulado a partir del primer año proyectado
    i = row_number() - 1,
    
    # Escenario Pesimista: -3% anual acumulado sobre la base
    Pesimista_Millones = round(Base_Millones * ((1 - 0.03)^i), 2),
    
    # Escenario Base: Redondeo a 2 decimales
    Base_Millones      = round(Base_Millones, 2),
    
    # Escenario Optimista: +3% anual acumulado sobre la base
    Optimista_Millones = round(Base_Millones * ((1 + 0.03)^i), 2)
  ) %>%
  select(Año, Meses, Pesimista_Millones, Base_Millones, Optimista_Millones)

# 4. Mostrar la tabla
kable(tabla_escenarios_millones, 
      col.names = c("Año", "Meses", "Pesimista (-3%/año)", "Escenario Base (SARIMAX)", "Optimista (+3%/año)"),
      caption = "Proyección Anual de Demanda por Escenarios (2026-2031) [Millones de Viajeros]")


