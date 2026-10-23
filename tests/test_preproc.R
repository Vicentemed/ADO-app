# tests/test_preproc.R — prueba de las utilidades de preprocesamiento (sin Shiny)
# Uso: Rscript tests/test_preproc.R   (desde CDO-app)
for (f in c("00_config.R", "02_data.R")) source(file.path("R", f), encoding = "UTF-8")
# cargar solo las funciones puras del módulo (sin UI)
src <- readLines(file.path("R", "08_preproc.R"), encoding = "UTF-8")
end <- grep("^# ---- UI", src)[1] - 1
eval(parse(text = src[1:end]))
ok <- function(cond, msg) cat(if (isTRUE(cond)) "OK   " else "FALLA", msg, "\n")

# 1) NA summary + eliminación
d <- read_book_data("data/datBF.txt"); d$Basal[c(2, 5)] <- NA
ns <- na_summary(d)
ok(ns[ns$Variable == "Basal", "n NA"] == 2, "na_summary cuenta NA")
ok(nrow(d[complete.cases(d), ]) == nrow(d) - 2, "eliminar filas con NA (cualquier variable)")
ok(nrow(d[complete.cases(d["Group"]), ]) == nrow(d), "eliminar filas con NA solo en variable sin NA")

# 2) subset por expresión y por niveles
d2 <- read_book_data("data/DXKTw.txt")
keep <- eval(parse(text = 'Grupo == 1 & Ed > 20')[[1]], d2)
ok(all(d2[keep, "Grupo"] == 1) && all(d2[keep, "Ed"] > 20), "subset por expresión R")
ok(nrow(d2[as.character(d2$Sexo) %in% "F", ]) + nrow(d2[as.character(d2$Sexo) %in% "M", ]) == nrow(d2), "subset por categorías")

# 3) wide → long (DXKTw: X0..X2 por Px) y vuelta
long <- pp_wide_to_long(d2, c("X0", "X1", "X2"), ids = c("Px", "Grupo"))
ok(nrow(long) == 3 * nrow(d2) && all(c("Px", "Grupo", "Times", "Response") %in% names(long)), "wide→long dimensiones y columnas")
ok(is.factor(long$Times) && identical(levels(long$Times), c("X0", "X1", "X2")), "wide→long conserva orden de tiempos")
ok(all(long$Response[long$Times == "X1"] == d2$X1), "wide→long valores correctos")
wide <- pp_long_to_wide(long, id = "Px", key = "Times", value = "Response")
ok(nrow(wide) == nrow(d2) && all(c("Response_X0", "Response_X1", "Response_X2") %in% names(wide)), "long→wide dimensiones y nombres <valor>_<nivel>")
ok(all(wide$Response_X2[match(d2$Px, wide$Px)] == d2$X2), "long→wide round-trip exacto")

# 4) wide → long sin id crea Px
long2 <- pp_wide_to_long(read_book_data("data/Carbuse.txt")[c("Ra", "ID")], c("Ra", "ID"))
ok("Px" %in% names(long2), "wide→long sin id crea Px")
