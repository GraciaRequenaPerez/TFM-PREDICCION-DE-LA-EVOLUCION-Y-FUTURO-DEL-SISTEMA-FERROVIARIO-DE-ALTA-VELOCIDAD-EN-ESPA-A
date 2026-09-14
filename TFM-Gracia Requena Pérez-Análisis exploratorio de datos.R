
library(readxl)
Viajeros <- read_excel("C:/Users/graci/Desktop/TFM/Datos/DATOS AVE y EXO.xlsx")
View(Viajeros)
names(Viajeros)

viajeros_ts <- ts(Viajeros$ViajerosAVyLD,
                  start = c(1996,1),
                  frequency = 12)

###############################################################
############### ANALISIS EXPLORATORIO DE DATOS ###############
library(forecast)
library(ggplot2)

###Plot la serie temporal
col_principal <- "#832D51"
autoplot(viajeros_ts, color = col_principal, linewidth = 0.8) +
  labs(
    title = "Evolución Mensual de Viajeros",
    subtitle = "Serie histórica 1996 - 2025 AV + LD",
    x = "Año",
    y = "Número de Viajeros"
  ) +
  scale_y_continuous(labels = scales::label_number(big.mark = ".")) +
  theme_minimal()

# Métricas descriptivas sobre la serie completa (viajeros_ts)
estadisticos_aed <- data.frame(
  Estadístico = c("Media", "Desviación Estándar", "Mínimo", "Mediana", "Máximo"),
  Valor       = c(
    round(mean(viajeros_ts, na.rm = TRUE), 2),
    round(sd(viajeros_ts, na.rm = TRUE), 2),
    round(min(viajeros_ts, na.rm = TRUE), 2),
    round(median(viajeros_ts, na.rm = TRUE), 2),
    round(max(viajeros_ts, na.rm = TRUE), 2)
  )
)

# Visualizar tabla
print(estadisticos_aed)

###Descomposición de la serie temporal
decomp <- stl(viajeros_ts, s.window = "periodic")
plot(decomp)



### ANALISIS DE TENDENCIA
summary(decomp)
trend <- decomp$time.series[, "trend"]
autoplot(trend, color = "#832D51", linewidth = 0.8) +
  labs(
    title = "Tendencia de la Serie de Viajeros AV + LD",
    x = "Año",
    y = "Valor"
  ) +
  theme_minimal()

##Analisis de puntos de quiebre
library(strucchange)

tiempo <- 1:length(viajeros_ts)

bp <- breakpoints( viajeros_ts~ tiempo, h = 0.1)

plot(bp)
summary(bp)

##Según el BIC mínimo hay 3 puntos de quiebre
plot(viajeros_ts, 
  col = "#832D51", 
  lwd = 1.2,
  xlab = "Año",
  ylab = "Viajeros AV + LD")

lines(fitted(bp, breaks = 3), col = "#447A5F", lwd = 1)

#1 Quiebre COVID
##2 2013 - Linea Alicante, conexión con Francia e introducción del sistema de precios dinamicos que redujo los precios de los billetes
###3 2008 Madrid Barcelona

###### ANALISIS ESTACIONALIDAD######
seasonal <- decomp$time.series[, "seasonal"]
plot(seasonal)
autoplot(seasonal, color = "#832D51", linewidth = 0.5) +
  labs(
    title = "Estacionalidad de la Serie de Viajeros AV + LD",
    x = "Año",
    y = "Valor"
  ) +
  theme_minimal()

library(ggplot2)
library(patchwork) # Para juntar los dos gráficos en 1 fila

# ==============================================================================
# 1. PREPARACIÓN DE DATOS
# ==============================================================================
meses_labels <- c("Ene", "Feb", "Mar", "Abr", "May", "Jun", 
                  "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")

df_estacional <- data.frame(
  Valor   = as.numeric(viajeros_ts),
  Mes_Num = cycle(viajeros_ts)
)
df_estacional$Mes <- factor(df_estacional$Mes_Num, levels = 1:12, labels = meses_labels)

# Resumen para el perfil estacional medio (media y desviación estándar)
df_perfil_medio <- aggregate(Valor ~ Mes, data = df_estacional, FUN = function(x) {
  c(media = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
})
df_perfil_medio <- data.frame(
  Mes   = df_perfil_medio$Mes,
  Media = df_perfil_medio$Valor[, "media"],
  SD    = df_perfil_medio$Valor[, "sd"]
)

# ==============================================================================
# 2. GRÁFICO IZQUIERDA: BOXPLOT MENSUAL
# ==============================================================================
p1 <- ggplot(df_estacional, aes(x = Mes, y = Valor)) +
  geom_boxplot(fill = "#447A5F", color = "#335B47", alpha = 0.7, outlier.size = 1.5) +
  labs(
    title = "Distribución mensual (boxplot)",
    x = NULL,
    y = "Viajeros"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank()
  )

# ==============================================================================
# 3. GRÁFICO DERECHA: PERFIL ESTACIONAL EN LÍNEA + BARRAS DE ERROR
# ==============================================================================
y_limits <- c(0, 3500)
y_breaks <- seq(0, 3500, by = 500)

p2 <- ggplot(df_perfil_medio, aes(x = Mes, y = Media, group = 1)) +
  # Barras de error (desviación estándar)
  geom_errorbar(aes(ymin = Media - SD, ymax = Media + SD), width = 0.2, color = "#333333") +
  # Línea conectando las medias de cada mes
  geom_line(color = "#832D51", size = 1) +
  # Puntos en las medias
  geom_point(color = "#832D51", size = 2.5) +
  scale_y_continuous(limits = y_limits, breaks = y_breaks) +
  labs(
    title = "Perfil estacional medio",
    x = NULL,
    y = "Viajeros medios"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank()
  )

# ==============================================================================
# 4. UNIR AMBOS PANELES
# ==============================================================================
p2 + p1

#### Estimación del parámetro lambda de Box-Cox para saber si es aditiva o multiplicativa ###
lambda <- BoxCox.lambda(viajeros_ts)
cat("Valor de Lambda estimado:", round(lambda, 3), "\n")
## Lambda = 1,243 Serie aditiva, no hace falta realizar ninguna transformación logaritmica.


#####Estacionariedad####
library(tseries)
library(forecast)

# 1. Test Augmented Dickey-Fuller (ADF) Si p-value > 0.05  No estacionaria
adf_res <- adf.test(viajeros_ts, alternative = "stationary")
print(adf_res)
#Para sorpresa de nadie, p > 0.05, así que son NO ESTACIONARIA

# 2. Test KPSS Si p-value <  0.05  No estacionaria
kpss_res <- kpss.test(viajeros_ts, null = "Trend")
print(kpss_res)
#Para sorpresa de nadie, p < 0.05, así que son NO ESTACIONARIA

########ACF y PACF###############

library(forecast)
library(gridExtra)

g1 <- ggAcf(viajeros_ts, lag.max = 36) + 
  labs(title = "ACF - Viajeros AV+LD")

g2 <- ggPacf(viajeros_ts, lag.max = 36) + 
  labs(title = "PACF - Viajeros AV+LD")

grid.arrange(g1, g2, ncol = 2)

library(ggplot2)

############OTROS DATOS·############################
# ==============================================================================
# 1. CÁLCULO DE LA VARIACIÓN INTERANUAL (%) 
# ==============================================================================

# Serie original en vector numérico
y <- as.numeric(viajeros_ts)
n <- length(y)

# T_t comparado con T_{t-12} (mismo mes del año anterior)
y_actual <- y[13:n]
y_hace_un_ano <- y[1:(n - 12)]

# Calculamos el porcentaje
var_val <- ((y_actual - y_hace_un_ano) / y_hace_un_ano) * 100

# Extraemos las fechas correspondientes a partir del mes 13 (año 2)
fechas_var <- as.Date(time(viajeros_ts))[13:n]

# Crear data frame limpio
df_var <- data.frame(
  Fecha = fechas_var,
  Var   = var_val,
  Signo = ifelse(var_val >= 0, "Positivo", "Negativo")
)

#=======================================================================
# 2.VARIACIÓN INTERANUAL DE VIAJEROS (%) CON MÚLTIPLES HITOS DESTACADOS
#=======================================================================

# Hito 1: Apertura Madrid-Barcelona
h1_inicio <- as.Date("2008-03-01")
h1_fin    <- as.Date("2009-03-01")

# Hito 2: Introducción precios dinámicos
h2_inicio <- as.Date("2013-03-01")
h2_fin    <- as.Date("2014-03-01")

# Hito 3: Impacto COVID-19
h3_inicio <- as.Date("2020-03-01")
h3_fin    <- as.Date("2022-03-01")

ggplot(df_var, aes(x = Fecha, y = Var, fill = Signo)) +
  
  # --- HITO 1: MADRID-BARCELONA ---
  annotate("rect", xmin = h1_inicio, xmax = h1_fin, ymin = -100, ymax = 250, 
           alpha = 0.2, fill = "#A6C048") +
  geom_vline(xintercept = as.numeric(h1_inicio), linetype = "dashed", color = "#A6C048", linewidth = 0.7) +
  
  annotate("text", x = h1_inicio - 100, y = 85, label = "LAV Madrid-BCN", 
           size = 6, fontface = "bold", color = "#335B47", hjust = 1) +
  # --- HITO 2: PRECIOS DINAMICOS ---
  annotate("rect", xmin = h2_inicio, xmax = h2_fin, ymin = -100, ymax = 250, 
           alpha = 0.2, fill = "#A6C048") +
  geom_vline(xintercept = as.numeric(h2_inicio), linetype = "dashed", color = "#A6C048", linewidth = 0.7) +
  
  annotate("text", x = h2_inicio - 100, y = 85, label = "Precios dinámicos", 
           size = 6, fontface = "bold", color = "#335B47", hjust = 1) +
  
  # --- HITO 3: IMPACTO COVID-19 ---
  annotate("rect", xmin = h3_inicio, xmax = h3_fin, ymin = -100, ymax = 250, 
           alpha = 0.25, fill = "#A6C048") +
  geom_vline(xintercept = as.numeric(h3_inicio), linetype = "dashed", color = "#A6C048", linewidth = 0.7) +
  annotate("text", x = h3_inicio - 100, y = -85, label = "Shock COVID-19", 
           size = 6, fontface = "bold", color = "#335B47", hjust = 1) +
  
  # --- BARRAS DE VARIACIÓN INTERANUAL ---
  geom_bar(stat = "identity", width = 25) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Positivo" = "#447A5F", "Negativo" = "#832D51")) +
  
  # --- EJES Y LIMITES ---
  coord_cartesian(ylim = c(-100, 200)) + 
  scale_y_continuous(breaks = seq(-100, 200, by = 25)) +
  
  # --- TÍTULOS Y ESTILO ---
  labs(
    title = "Variación interanual de viajeros (%)",
    subtitle = "Eventos destacados: Inauguración LAV Madrid-Barcelona (2008), Introducción precios dinámicos (2013) y Shock COVID-19 (2020)",
    x = NULL,
    y = "% var. interanual"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 24, face = "bold"),
    plot.subtitle = element_text(size = 18, color = "gray30", margin = margin(b = 8)),
    axis.title   = element_text(size = 18),
    axis.text   = element_text(size = 16),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )


library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(patchwork)

df_exog_anual <- Viajeros %>%
  mutate(
    # Intenta convertir a fecha primero; si falla o si ya es texto/fecha, extrae los 4 dígitos del año
    Anyo = case_when(
      inherits(Mes, "Date") ~ as.numeric(format(Mes, "%Y")),
      inherits(Mes, "POSIXt") ~ as.numeric(format(Mes, "%Y")),
      TRUE ~ as.numeric(gsub("^.*?(\\d{4}).*", "\\1", as.character(Mes)))
    )
  ) %>% 
  group_by(Anyo) %>%
  summarise(
    Km_Linea   = max(KM, na.rm = TRUE),
    Materiales = max(Materiales, na.rm = TRUE)
  ) %>%
  filter(!is.na(Anyo)) %>% # Elimina posibles filas vacías al final del Excel
  ungroup()

# ==============================================================================
# 2. GRÁFICO 1 (IZQUIERDA): KILÓMETROS DE LÍNEA (KM)
# ==============================================================================
p_km <- ggplot(df_exog_anual, aes(x = Anyo, y = Km_Linea)) +
  geom_col(fill = "#832D51", width = 0.7, alpha = 0.85) +
  geom_line(color = "#447A5F", linewidth = 0.9, group = 1) +
  geom_point(color = "#447A5F", size = 2) +
  scale_x_continuous(breaks = seq(min(df_exog_anual$Anyo), max(df_exog_anual$Anyo), by = 5)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "A. Evolución de la red (KM)",
    subtitle = "Kilometros a cierre de ejercicio",
    x = NULL,
    y = "Kilómetros"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "black", margin = margin(b = 4)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
  )

# ==============================================================================
# 3. GRÁFICO 2 (DERECHA): MATERIALES
# ==============================================================================
# 1. Creación del desglose completo desde 1996
df_mat_desglosado <- df_exog_anual %>%
  mutate(
    # Flota de Ouigo según año de entrada
    Ouigo = case_when(
      Anyo == 2021 ~ 4,
      Anyo >= 2022 & Anyo <= 2025 ~ 20,
      TRUE ~ 0
    ),
    # Flota de Iryo según año de entrada
    Iryo = case_when(
      Anyo == 2022 ~ 16,
      Anyo >= 2023 & Anyo <= 2025 ~ 18,
      TRUE ~ 0
    ),
    # Renfe asume la totalidad de la serie desde 1996 hasta 2020, 
    # y la diferencia del total a partir de 2021
    Renfe = Materiales - Ouigo - Iryo
  ) %>%
  # Pivotamos a formato largo para el geom_col apilado
  pivot_longer(
    cols = c(Renfe, Ouigo, Iryo),
    names_to = "Operador",
    values_to = "Unidades"
  ) %>%
  # Orden visual de las barras: Renfe abajo, luego Ouigo e Iryo encima
  mutate(Operador = factor(Operador, levels = c("Iryo", "Ouigo", "Renfe")))

# 2. Gráfico histórico completo (1996 - 2025)
p_mat <- ggplot() +
  # Barras apiladas (para 1996-2020 serán 100% Renfe; en 2021-2025 se dividen)
  geom_col(
    data = df_mat_desglosado,
    aes(x = Anyo, y = Unidades, fill = Operador),
    width = 0.7, 
    alpha = 0.9
  ) +
  # Línea y puntos del total global sobre toda la serie histórica
  geom_line(
    data = df_exog_anual, 
    aes(x = Anyo, y = Materiales), 
    color = "#447A5F", 
    linewidth = 0.9, 
    group = 1
  ) +
  geom_point(
    data = df_exog_anual, 
    aes(x = Anyo, y = Materiales), 
    color = "#447A5F", 
    size = 2
  ) +
  # Paleta de colores para los operadores
  scale_fill_manual(
    values = c(
      "Renfe" = "#832D51", 
      "Ouigo" = "#0096C9",
      "Iryo"  = "red"  
    )
  ) +
  # Marcas del eje X cada 4 o 5 años para cubrir desde 1996 a 2025
  scale_x_continuous(
    breaks = seq(1996, max(df_exog_anual$Anyo), by = 4)
  ) +
  scale_y_continuous(labels = scales::label_number(big.mark = ".", decimal.mark = ",")) +
  labs(
    title = "B. Evolución de la flota (ud)",
    subtitle = "Material rodante de los operadores",
    x = NULL,
    y = "Nº de materiales",
    fill = "Operador"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, margin = margin(b = 4)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.text = element_text(size = 9)
  )


# ==============================================================================
# 4. MOSTRAR PANELES EN PARALELO
# ==============================================================================
p_km + p_mat


###Plot la serie desempleo
parados_ts <- ts(Viajeros$ParadosESP,
                  start = c(1996,1),
                  frequency = 12)
col_principal <- "#832D51"
autoplot(parados_ts, color= col_principal, linewidth = 1) +
  labs(
    title = "Evolución Mensual de Desempleo",
    subtitle = "Serie histórica 1996 - 2025",
    x = "Año",
    y = "Número de Desempleados (en miles)"
  ) +
  scale_y_continuous(labels = scales::label_number(big.mark = ".")) +
  theme_minimal()+
theme(
  plot.title    = element_text(size = 18, face = "bold"),
  plot.subtitle = element_text(size = 12, margin = margin(b = 8)),)



###VARIABLES DE CALENDARIO
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Asignar el factor de Mes a partir de la estructura temporal de viajeros_ts
Viajeros_df <- Viajeros %>%
  mutate(
    Mes_num = cycle(viajeros_ts),
    Mes = factor(Mes_num, levels = 1:12, 
                 labels = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", 
                            "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"))
  )

# 2. Agregación mensual de ViajerosAVyLD y las 3 variables de festivos
df_perfil_completo <- Viajeros_df %>%
  group_by(Mes) %>%
  summarise(
    Media = mean(ViajerosAVyLD, na.rm = TRUE),
    SD = sd(ViajerosAVyLD, na.rm = TRUE),
    FestivosN = mean(FestivosN, na.rm = TRUE),
    FestivosM = mean(FestivosM, na.rm = TRUE),
    FestivosCCVA = mean(FestivosCCVA, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Pivotar festivos para generar las barras apiladas
df_festivos_medio <- df_perfil_completo %>%
  pivot_longer(
    cols = c(FestivosN, FestivosM, FestivosCCVA),
    names_to = "Tipo_Festivo",
    values_to = "Media_Dias"
  ) %>%
  mutate(
    Tipo_Festivo = factor(
      Tipo_Festivo,
      levels = c("FestivosN", "FestivosM", "FestivosCCVA"),
      labels = c("Nacionales", "Madrid", "CAT/CCVA/AND")
    )
  )

# 4. Ajustes visuales y factor de conversión para el eje secundario
y_limits <- c(0, 3500)
y_breaks <- seq(0, 3500, by = 500)
factor_escala <- 3500 / 5  # Mapea un máximo de ~5 festivos al techo de 3500

# 5. Construcción del gráfico p2 integrado
p2 <- ggplot() +
  # Fondo: Barras apiladas de la media de festivos
  geom_col(
    data = df_festivos_medio,
    aes(x = Mes, y = Media_Dias * factor_escala, fill = Tipo_Festivo),
    alpha = 0.45,
    width = 0.55
  ) +
  # Barras de error (Desviación estándar de ViajerosAVyLD)
  geom_errorbar(
    data = df_perfil_completo,
    aes(x = Mes, ymin = Media - SD, ymax = Media + SD),
    width = 0.2, color = "#333333"
  ) +
  # Línea conectando las medias de viajero
  geom_line(
    data = df_perfil_completo,
    aes(x = Mes, y = Media, group = 1),
    color = "#832D51", size = 1
  ) +
  # Puntos en las medias de viajero
  geom_point(
    data = df_perfil_completo,
    aes(x = Mes, y = Media),
    color = "#832D51", size = 2.5
  ) +
  # Ejes primario (Viajeros) y secundario (Festivos)
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    sec.axis = sec_axis(~ . / factor_escala, name = "Días festivos medios")
  ) +
  # Paleta de colores ajustada a la serie borgoña (#832D51)
  scale_fill_manual(
    values = c(
      "Nacionales" = "#447A5F",
      "Madrid" = "#EA6993",
      "CAT/CCVA/AND" = "#ACC556"
    )
  ) +
  labs(
    title = "B. Perfil estacional y distribución anual de festivos",
    x = NULL,
    y = "Viajeros medios",
    fill = "Tipo de Festivo"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size =8 )
  )

library(patchwork)

# 1. Agregación mensual de Viajeros y la variable 'vacaciones'
df_perfil_vacaciones <- Viajeros_df %>%
  group_by(Mes) %>%
  summarise(
    Media = mean(ViajerosAVyLD, na.rm = TRUE),
    SD = sd(ViajerosAVyLD, na.rm = TRUE),
    Vacaciones_Media = mean(Vacaciones, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Factor de conversión (mapea el valor 1.0 de vacaciones al techo de 3500)
factor_escala_vac <- 3500 / 1

# 3. Construcción del gráfico p1
p1 <- ggplot(df_perfil_vacaciones) +
  # Fondo: Barras de intensidad vacacional (0, 0.5, 1)
  geom_col(
    aes(x = Mes, y = Vacaciones_Media * factor_escala_vac),
    fill = "#447A5F",
    alpha = 0.35,
    width = 0.55
  ) +
  # Barras de error (Desviación estándar)
  geom_errorbar( aes(x = Mes, ymin = Media - SD, ymax = Media + SD), width = 0.2, color = "#333333") +
  geom_line( aes(x = Mes, y = Media, group = 1), color = "#832D51", size = 1 ) +
  geom_point(aes(x = Mes, y = Media), color = "#832D51", size = 2.5) +
  scale_y_continuous(
    limits = c(0, 3500),
    breaks = seq(0, 3500, by = 500),
  ) +
  labs(
    title = "A. Perfil estacional y distribución de vacaciones", 
    x = NULL,
    y = "Viajeros medios"  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank(),
  )

# 4. Combinar con p2 (asumiendo que p2 ya está creado en tu sesión)
grafico_compuesto <- p1 + p2
print(grafico_compuesto)


library(ggcorrplot)
library(dplyr)

# Seleccionar solo las variables numéricas exógenas y la variable objetivo
matriz_vars <- Viajeros_df %>%
  select(
    ViajerosAVyLD,
    KM,
    Materiales,
    Compañías,
    ParadosESP,
    Vacaciones,
    FestivosN,
    FestivosM,
    FestivosCCVA,
    COVID,
    PrDin,
    Crisis2008,
    Accidentes
    
  ) %>%
  cor(use = "complete.obs", method = "pearson")

# Generar el mapa de calor de correlaciones
p_corr <- ggcorrplot(
  matriz_vars,
  hc.order = FALSE,
  type = "lower",
  lab = TRUE,
  lab_size = 3.5,
  method = "square",
  colors = c("#447A5F", "#ECF1D8", "#832D51"),
  title = "Matriz de Correlaciones de Pearson (Exógenas vs Viajeros)",
  ggtheme = theme_minimal()+
  theme(
    plot.title    = element_text(size = 18, face = "bold"))
)

print(p_corr)

