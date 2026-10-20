# AACclustWind
######## En R ##########
# Ver si el paquete existe en la biblioteca
"AACclustWind" %in% rownames(installed.packages())

###
remove.packages("AACclustWind")


setwd("C:\\.......\\AACclustWind")

#devtools::document()

# 2. Volver a hacer el check (debería dar 0 errors, 0 warnings)
#devtools::check()

# 3. Instalar el paquete formalmente (CRUCIAL para que funcione el paralelo)
devtools::install()

# Listar todas las funciones exportadas
ls("package:AACclustWind")
###### Para  detectar el numero de los cores.  ######
library(parallel)

# Detectar todos los cores lógicos
n_cores <- detectCores()
print(n_cores)

# Detectar solo cores físicos (excluye hyperthreading)
n_cores_fisicos <- detectCores(logical = FALSE)
print(n_cores_fisicos)


# 4. Probar la función que te fallaba
library(AACclustWind)
data(Felinos)
resultado <- aac_parallel(data = Felinos, k_range = 2:5, n_cores = 4)
print(resultado)


