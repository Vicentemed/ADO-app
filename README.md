# ADO-app · Asistente de análisis guiado

Aplicación web en **R / Shiny** que acompaña al libro **_Análisis de Datos en Odontología. Un enfoque pragmático con R_**. Reúne todo el código estadístico usado en el libro en un asistente paso a paso: el usuario elige el tipo de análisis, carga sus datos (o usa los ejemplos del libro), selecciona las variables, y la app verifica los supuestos, decide la prueba adecuada siguiendo la regla del libro, ejecuta el análisis, genera los gráficos y produce un reporte en Word listo para editar.

## Contenido del repositorio

```
ADO-app/
├── app.R                  # Interfaz y servidor Shiny (flujo de 9 pasos)
├── R/
│   ├── 00_config.R        # Configuración global: alfa, paleta, formateadores, tamaños de efecto
│   ├── 01_registry.R      # Registro genérico de análisis (ado_register) y regla de decisión
│   ├── 02_data.R          # Lectura de datos (.txt/.tsv/.csv/.xlsx), saneado de nombres, roles de variables
│   ├── 03_assumptions.R   # Verificación de supuestos (normalidad, homogeneidad de varianzas, etc.)
│   ├── 04_descriptives.R  # Estadística descriptiva (table1)
│   ├── 05_plots.R         # Gráficos con ggplot2 y parámetros cosméticos
│   ├── 06_report.R        # Construcción del reporte Word (rmarkdown)
│   ├── 07_ui_helpers.R    # Componentes reutilizables de la interfaz
│   ├── 08_preproc.R       # Módulo de preprocesamiento de datos
│   └── analyses/          # Un archivo por familia de análisis
│       ├── medias.R
│       ├── correlacion.R
│       ├── categoricos.R
│       ├── concordancia.R
│       ├── regresion.R
│       └── supervivencia.R
├── data/                  # Conjuntos de datos de ejemplo del libro (tabulados, .txt)
├── report/
│   ├── report_template.Rmd    # Plantilla del reporte Word
│   └── Ejemplo_*.docx         # Reportes de ejemplo generados por la app
├── tests/                 # Pruebas y generadores de reportes de ejemplo (ver abajo)
└── www/style.css          # Estilos de la interfaz
```

## Análisis disponibles

Cada análisis está vinculado al capítulo del libro donde se explica con un conjunto de datos de ejemplo precargado. Cuando la prueba paramétrica tiene una alternativa no paramétrica, la app la propone automáticamente si se viola algún supuesto.

| Familia | Análisis | Alternativa no paramétrica |
|---|---|---|
| **Estadística descriptiva** | Tabla 1 con la lógica de `table1(~ respuestas \| explicativas)`: una o más variables de respuesta en las filas y una o más variables explicativas o de grupo en las columnas (dos se anidan, como en `~ VH \| Group*Cure.time`) | - |
| **Comparación de medias** | Prueba t independiente | Wilcoxon suma de rangos (Mann-Whitney) |
| | Prueba t apareada | Wilcoxon de rangos con signo |
| | ANOVA de una vía | Kruskal-Wallis |
| | ANOVA de medidas repetidas (ezANOVA con corrección GG/HF, o modelo multinivel `lme()` con razón de verosimilitud) | ANOVA de Friedman |
| | Curvas de crecimiento por grupo (`drc`, razón de verosimilitud) | - |
| | ANOVA factorial (dos factores) | - |
| | ANCOVA (grupo + covariable) | - |
| **Correlación** | Correlación de Pearson | Correlación de Spearman |
| | Tau de Kendall | - |
| | Correlación parcial (controlando Z) | - |
| | Predictive Power Score (PPS) | - |
| **Datos categóricos** | Chi-cuadrado de Pearson | Prueba exacta de Fisher |
| | Prueba de McNemar (apareada) | - |
| **Concordancia** | Kappa de Cohen (ponderado si ordinal; Fleiss con > 2 evaluadores) | - |
| | Coeficiente de correlación intraclase (ICC) | - |
| | Bland-Altman | Coeficiente de concordancia de Lin (CCC) |
| **Regresión** | Regresión lineal múltiple | - |
| | Regresión logística | - |
| **Supervivencia** | Kaplan-Meier y prueba log-rank | - |
| | Regresión de Cox | - |

Los post-hoc y contrastes se eligen en el paso 6, según el análisis: Tukey, Dunnett (con grupo control a elección), Bonferroni, Holm, Benjamini-Hochberg, Games-Howell, contrastes ortogonales o personalizados y análisis de tendencia en el ANOVA de una vía; Siegel y Castellan o Wilcoxon por pares en Kruskal-Wallis; `friedmanmc` en Friedman; residuos estandarizados ajustados en tablas mayores de 2×2.

Las reglas de decisión siguen el libro: la normalidad se evalúa con Shapiro-Wilk y, con n > 50, con el sesgo y la curtosis estandarizados (|valor| > 1); chi-cuadrado frente a Fisher se decide solo por las frecuencias esperadas (alguna < 5 → Fisher); en medidas repetidas, si la prueba de Mauchly es significativa se aplica la corrección de Greenhouse-Geisser o Huynh-Feldt (o el modelo multinivel); Friedman queda para respuestas ordinales o marcadamente no normales.

## Flujo de uso

La barra lateral permite elegir la **familia** y el **análisis**, y la fuente de datos. El panel central guía por nueve pasos; los botones **Anterior / Siguiente** avanzan entre ellos.

1. **Datos.** Elija *Ejemplo del libro* (se carga el conjunto de datos del capítulo) o *Subir mi archivo* (`.txt` o `.tsv` separados por tabulador, `.csv`, `.xlsx`). Los nombres de columna con espacios o símbolos se corrigen automáticamente y se avisa del cambio.
2. **Preprocesar.** Operaciones opcionales que se aplican en orden y quedan registradas: eliminar filas o variables con valores faltantes, remover variables, filtrar filas (por categorías o por expresión), convertir entre formato ancho y largo, y cambiar el tipo de una variable (factor o numérica). Se puede deshacer la última operación o restablecer los datos originales.
3. **Variables.** Asigne cada variable al rol que pide el análisis (respuesta, grupo, mediciones, sujeto, etc.). La app valida la selección e indica qué tipo requiere cada rol. Opcionalmente, elija una **variable de estratificación**: todo el análisis se repite por separado dentro de cada nivel y los resultados se muestran en pestañas.
4. **Descriptiva.** Tabla descriptiva (table1) de las variables del análisis. En el panel de la izquierda se pueden elegir otras variables de respuesta (filas) y explicativas o de grupo (columnas), siguiendo `table1(~ respuestas | explicativas)`; la tabla elegida es la que entra al reporte.
5. **Supuestos.** Verificación de supuestos y decisión según la regla del libro, con el gráfico Q-Q (`car::qqPlot`) de la respuesta, de las diferencias apareadas o de los residuos, según el análisis. Puede dejar la ruta *automática* o forzar la prueba paramétrica o la alternativa.
6. **Análisis.** Ajuste las opciones de la prueba y pulse **Ejecutar análisis**. Se muestran las tablas de resultados.
7. **Post-hoc.** Contrastes múltiples cuando el análisis los requiere.
8. **Gráficos.** Elija el gráfico y ajuste colores, títulos y otros parámetros cosméticos.
9. **Reporte.** Escriba el título y el **contexto de los datos** (objetivo, población, diseño, significado de las variables), marque las secciones a incluir y pulse **Generar reporte Word**. La app también muestra un texto de *Cómo reportar* prellenado con los resultados, donde los corchetes `[ ]` señalan lo que debe completar el usuario.

## Requisitos e instalación

Se necesita **R** (4.x o superior) y **Pandoc** para generar el reporte Word. Pandoc viene incluido con RStudio; si la app se ejecuta fuera de RStudio, `R/00_config.R` intenta localizarlo en las rutas habituales de RStudio o Quarto en Windows.

Instale los paquetes necesarios:

```r
install.packages(c(
  # Interfaz y reporte
  "shiny", "bslib", "DT", "ggplot2", "htmltools", "scales",
  "rmarkdown", "knitr", "flextable", "table1", "readxl",
  # Análisis
  "car", "multcomp", "ez", "effects", "pgirmess", "psych", "pastecs", "pwr",
  "irr", "irrCAC", "epiR", "ggm", "ppsr", "pROC", "ResourceSelection",
  "nlme", "drc", "gam", "leaps", "caret", "lsr", "rms", "qcc", "dcurves",
  "survival", "survminer"
))

# ggm necesita el paquete graph de Bioconductor
install.packages("BiocManager"); BiocManager::install("graph")
```



## Ejecución

Desde la carpeta del proyecto:

```r
shiny::runApp()
```

O desde la terminal:

```bash
Rscript -e "shiny::runApp('.', launch.browser = TRUE)"
```

## Pruebas y reportes de ejemplo

Los scripts de `tests/` se ejecutan desde la carpeta del proyecto con `Rscript`:

| Script | Qué hace |
|---|---|
| `tests/smoke_test.R` | Ejecuta **cada análisis registrado** sobre su ejemplo del libro: descriptiva, supuestos, análisis, post-hoc, gráficos y texto de reporte. |
| `tests/verificacion_libro/run_all.R` | Corre todos los análisis sobre los ejemplos del libro y compara cada cifra con la impresa en el libro; escribe `resultados.csv`. Requiere el texto del libro en esa carpeta (no se incluye en el repositorio). |
| `tests/test_preproc.R` | Prueba las utilidades de preprocesamiento sin Shiny. |
| `tests/make_example_report.R` | Genera un reporte Word de ejemplo (ANOVA de una vía) en `report/`. |
| `tests/make_preproc_report.R` | Genera un reporte Word con historial de preprocesamiento. |
| `tests/test_stratified.R` | Prueba del análisis estratificado. Requiere un archivo `.xlsx` local que no está en el repositorio. |
| `tests/test_nonsyntactic_names.R` | Prueba con nombres de columna no válidos en R. Requiere el mismo archivo local. |

```bash
Rscript tests/smoke_test.R
```

## Cómo agregar un análisis

La interfaz es genérica: cada análisis se **declara** con `ado_register()` en un archivo de `R/analyses/`. La declaración indica la familia, el capítulo, los roles de variables que pide, el ejemplo del libro, y las funciones de descriptiva, supuestos, ejecución, post-hoc, gráficos y texto de reporte. El esquema completo está documentado al inicio de `R/01_registry.R`. Al agregar una entrada, la app la muestra automáticamente en el selector y el `smoke_test.R` la ejecuta.

## Autor

Dr. Vicente Esparza-Villalpando.

Universidad Autónoma de Aguascalientes

Centro de Ciencias de la Salud

Departamento de Estomatología
