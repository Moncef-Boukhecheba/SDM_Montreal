# Decommenter cette ligne pour installer le package BlockCV (introuvable sur Anaconda)
# install.packages('blockCV')

# Importer les bibliotheques necessaires 

# Tidyverse : une bibliotheque populaire de manipulation de donnees
# Elle contient plusieurs autres packages telles que dplyr, readr, ggplot...etc.
library('tidyverse')

# Sp : Fournit des classes et des methodes permettant de representer et manipuler les donnees geospatiales
library('sp')

# Raster : Api haut niveau permettant la lecture / la manipulation des donnees raster et vecteurs
library('raster')

# RandomForest : Permet la creation d'un modele Random Forest
library('randomForest')

# RGDAL : Implementation Haut niveau de la bibliotheque publique GDAL
library('rgdal')

# Caret : Bibliotheque d'evalation des modeles de Machine Learning
library('caret')

# Conflicted : Affiche les conflits entre les packages si il y en a (lors de l'appel d'une fonction par exemple)
library('conflicted')

# Rgeos : Est utilisee ici pour creer des polygones utiles Ã  la generation de donnees de pseudo-absence
library('rgeos')

# BlockCV : Permet la separation des donnees de test/entrainement en utilsant la Block Cross Validation
library('blockCV')

# Dismo : Bibliotheque tres populaire permettant d'implementer des fonctions et des modeles associes aux SDMs
library('dismo')

# ROCR : Bibliotheque de generation des courbes ROC
library('ROCR')

# ========== PRETRAITEMENT DES DONNEES D'OCCURRENCE =============

# Lecture du fichier CSV
occurrence_data_raw <- read_csv('./data/arbres-publics.csv')

# Creation de la classe 
occurrence_data_raw$Present <- 1

# Selection des attributs pertinents
occurrence_data <- occurrence_data_raw %>% dplyr::select('Essence_latin', 'Latitude', 'Longitude', 'Present')
occurrence_data$Essence_latin <- as.factor(occurrence_data$Essence_latin)

# ===== 1. Valeurs manquantes =======
before_count <- sum(count(occurrence_data)$n)
occurrence_data <- occurrence_data %>% dplyr::filter((!is.na(Longitude) | Longitude == 0) & 
                                                     !is.na(Latitude) & !is.na(Essence_latin))
after_count <- sum(count(occurrence_data)$n)

# Afficher le nombre de lignes supprimees
sprintf("Nombre de lignes avec valeurs manquantes enlevees : %d", before_count - after_count)

# ===== 2. Enregistrements doublons =======
before_count <- sum(count(occurrence_data)$n)
occurrence_data <- occurrence_data %>% dplyr::distinct()
after_count <- sum(count(occurrence_data)$n)

# Afficher le nombre de lignes supprimees
sprintf("Nombre de lignes doublons enlevees : %d", before_count - after_count)

# ===== 3. Especes pas assez representees ======
threshold <- 10
before_count <- length(levels(occurrence_data$Essence_latin))
occurrence_data = occurrence_data %>%
      group_by(Essence_latin) %>%
      dplyr::filter(n() > threshold)
num_species <- length(levels(droplevels(occurrence_data$Essence_latin)))

# Afficher le nombre d'especes enlevees
sprintf("Nombre d'especes enlevees : %d", before_count - num_species)


# ======= CHARGEMENT DES VARIABLES GLOBALES ========
# Charger la carte de Montreal (Limites administratives)
montreal_map <- shapefile('data/montreal/limadmin/LIMADMIN.shp')

# Fichiers des données climatiques 
worldclim_path <- "data/worldclim_cropped"
climate_filenames <- list.files(worldclim_path, pattern="*.tif", full.names=TRUE)

# Charger une carte des données climatiques pour référence (génération des pseudo-absences)
climate_raster <- raster(file.path(worldclim_path, "wc2.1_30s_bio_1.tif"))

# Definir une valeur de 'seed' pour que les resultats soient les memes à  chaque tentative
set.seed(12)

# Fonction qui regroupe tout le processus d'entrainement (sera appelée pour chaque espèce)
train_sdm <- function (species_data, species_name, species_count){
    # Obtenir le "carre" entourant la ville de Montreal (minlong, maxlong, minlat, maxlat)
    montreal_boundaries <- extent(montreal_map)
    
    # Transformer les donnees d'occurrence en SpatialPointsDataFrame
    species_data <- SpatialPointsDataFrame(data=species_data[3], coords=species_data[2:1])
    
    # Filtrer les donnees incoherentes
    species_data <- crop(species_data, montreal_boundaries)
        
    # ============ GENERATION DE PSEUDO ABSENCES ==========
    # Generer les points d'absence
    r_layer <- raster(ext=montreal_boundaries, res=0.00833333333333342)
    raster_mask = raster::rasterize(montreal_map, r_layer, res=res(r_layer))
    absence_data = SpatialPoints(randomPoints(mask=raster_mask, n=length(species_data), p=species_data, tryf=100, 
                                              lonlatCorrection=TRUE))
    
    # Recentrer les points d'absence
    crs(absence_data) <- crs(montreal_map)
    absence_data_indexes <- over(absence_data, geometry(montreal_map))
    absence_data <- absence_data[!is.na(absence_data_indexes), ]
    
    # Transformer les donnees d'absence en SpatialPointsDataFrame (avec l'attribut Present à 0)
    absence_data <- SpatialPointsDataFrame(data = data.frame(Present = rep(0, length(absence_data))), coords=absence_data)
    
    # Combiner les donnees d'absence et de presence
    crs(absence_data) <- crs(species_data)
    species_data <- rbind(species_data, absence_data)
    
    # =========== DONNEES CLIMATIQUES ============        
    # Dataset complet 
    species_data_final <- species_data
    
    # La fonction qui retourne la moyenne des élements adjacents à chaque cellule
    fill.na <- function(x, i=5) {return(mean(x, na.rm=TRUE))} 
    
    # Pour chaque fichier...
    for(filename in climate_filenames){
        # Charger le fichier en memoire
        climate_raster <- raster(filename)
        
        # Enlever les valeurs manquantes
        climate_raster <- focal(climate_raster, w = matrix(1,5,5), fun = fill.na, NAonly=TRUE)
        
        # Stocker le nom de la variable (ex : bio_1)
        variable_name <- strsplit(filename, '_')[[1]]
        variable_name <- strsplit(variable_name[length(variable_name)], '.tif')
        variable_name <- paste('bio_', variable_name, sep="")
        
        # Extraire la variable climatique pour chaque point d'occurrence
        climate_variable <- raster::extract(climate_raster, species_data_final, col=variable_name, method="simple")
        
        # Ajouter la nouvelle feature dans le dataset
        species_data_final@data[variable_name] <- climate_variable
    }
    
    # ---- Suppression des valeurs manquantes / attributs constants -----

    # Enlever les valeurs manquantes
    length_before_na_omit <- count(species_data_final@data)$n
    species_data_final = species_data_final[apply(species_data_final@data, FUN=function(x){!any(is.na(x))}, MARGIN=1), ]
    length_after_na_omit <- count(species_data_final@data)$n
    print(paste('Nombre de valeurs manquantes enlevees :', length_before_na_omit - length_after_na_omit))
    
    # Identifier et enlever les variables qui ont 0 de variance
    zero_variance_variables <- which(apply(species_data_final@data, 2, var) == 0)
    
    num_zero_var_variables <- length(zero_variance_variables)
    print(paste('Nombre de variables avec zero variance enlevees :', num_zero_var_variables))
    if (num_zero_var_variables > 0){
        print(paste('Attribut supprimes', colnames(species_data_final@data[as.numeric(zero_variance_variables)])))
        species_data_final@data <- species_data_final@data[ -as.numeric(zero_variance_variables)]
    }
    
    # ---- Suppression des attributs redondants -----
    
    # Créer la matrice de correlation
    cor_matrix <- cor(species_data_final@data)
    
    # Eliminer le triangle superieur de la matrice de correlation (car elle est symetrique)
    cor_matrix[upper.tri(cor_matrix)] <- 0
    
    # Fixer les valeurs dans la diagonale car les variables avec elles-mÃªmes ont un coefficient de 1
    diag(cor_matrix) <- 0
      
    # Enlever les attributs corrélés                                                                 
    species_data_final@data <- species_data_final@data[, !apply(cor_matrix, 2, function(x) any(abs(x) > 0.9))]
    
    # ---- Suppression des enregistrements redondants -----

    # Compter le nombre d'enregistrements uniques dans le dataset
    length_sp_data <- count(species_data_final@data)
    print(paste("Nombre d'enregistrements uniques en incluant les classes :", count(species_data_final@data %>% distinct())$n, 
                "/", length_sp_data))
    
    # Compter le nombre d'enregistrements sans la valeur cible 
    print(paste("Nombre d'enregistrements uniques sans les classes :",
          count(species_data_final@data[2:length(species_data_final@data)] %>% distinct())$n, 
               "/", length_sp_data))
    
    # Compter le nombre d'enregistrements par classes
    print("==== Nombre d'enregistrements uniques par classe ====")
    print(paste("Present = 1 :", count(species_data_final@data[species_data_final@data$Present == 1, ] %>% distinct())$n))
    print(paste("Present = 0 :", count(species_data_final@data[species_data_final@data$Present == 0, ] %>% distinct())$n))
    
    # Enlever les données dupliquées dans le dataset
    species_data_final_trimmed <- species_data_final[!duplicated(species_data_final@data), ]
    
    # ========== ENTRAINEMENT DU MODELE =============
    
    # Definir le CRS du dataset (pour eviter une erreur de BlockCV)
    crs(climate_raster) <- NA
    crs(species_data_final_trimmed) <- NA
    
    # Strategie de block spatiaux 
    folds <- spatialBlock(speciesData = species_data_final_trimmed, # Donnees d'occurrence
                       species = "Present", # La colonne cible (la classe)
                       theRange = 5000, # Taille des blocks (en m)
                       k = 5, # Nombre de folds
                       # rasterLayer = climate_raster, # a raster for background (optional)
                       selection = "random", # Selection aleatoires de blocs 
                       iteration = 100, # Nombre d'iterations avant de trouver les blocs optimaux (bonne distribution des donnees)
                       biomod2Format = FALSE)
    
    # Metriques 
    # --- Entrainement ----
    train_accuracy <- numeric(length(folds$folds))
    train_kappa <- numeric(length(folds$folds))
    train_auc <- numeric(length(folds$folds))
    train_f1_measure <- numeric(length(folds$folds))
    train_no_information_rate <- numeric(length(folds$folds))
    
    # --- Test ----
    accuracy <- numeric(length(folds$folds))
    kappa <- numeric(length(folds$folds))
    auc <- numeric(length(folds$folds))
    f1_measure <- numeric(length(folds$folds))
    no_information_rate <- numeric(length(folds$folds))
    
    # Modeles et predictions
    rf <- vector("list", (length(folds$folds)))
    # predictions <- vector("list", (length(folds$folds)))
    
    # Pour chaque fold...
    for (i in 1:length(folds$folds)){
        fold <- folds$folds[[i]]
        
        # Extraire les donnees d'entrainement
        xy_train <- species_data_final_trimmed[fold[1][[1]], ]@data
        y_train <- xy_train$Present
        x_train <- xy_train[2:length(colnames(xy_train))]
        
        # Extraire les donnees de test
        xy_test <- species_data_final_trimmed[fold[2][[1]], ]@data
        
        # Entrainer le modele Random Forest (parametres par default)
        # sdm_rf <- randomForest(x_train, as.factor(y_train))
        sdm_rf <- tuneRF(x_train, as.factor(y_train), trace=FALSE, plot=FALSE,
                         stepFactor=2, improve=1e-6, ntree=501, doBest=TRUE) 
        rf[[i]] <- sdm_rf
        
        # ===== METRIQUES ======
        # Nous utilisons la fonction confusionMatrix de Caret 
        # pour obtenir les metriques suivantes : accuracy, kappa et no_information_rate
        
        # ------ ENTRAINEMENT ------
        prediction <- predict(sdm_rf, xy_train)
        
        # Accuracy / Kappa / NIR / F-Measure
        cm <- confusionMatrix(as.factor(xy_train$Present), as.factor(prediction))
        overall <- cm[3][[1]]
        train_accuracy[i] <- overall[1]
        train_kappa[i] <- overall[2]
        train_no_information_rate[i] <- overall[5]
        
        byClass <- cm[4][[1]]
        train_f1_measure[i] <- byClass[7]
        
        # AUC
        rf_p_train <- predict(sdm_rf, type="prob", newdata=xy_train)[,2]
        rf_pr_train <- prediction(rf_p_train, as.factor(xy_train$Present))
        train_auc[i] <- performance(rf_pr_train, measure = "auc")@y.values[[1]] 
        
        # ----- TEST ------
        prediction <- predict(sdm_rf, xy_test)
        
        cm <- confusionMatrix(as.factor(xy_test$Present), as.factor(prediction))
        overall <- cm[3][[1]]
        accuracy[i] <- overall[1]
        kappa[i] <- overall[2]
        no_information_rate[i] <- overall[5]
        
        byClass <- cm[4][[1]]
        f1_measure[i] <- byClass[7]
        
        # Calculer l'AUC
        rf_p_train <- predict(sdm_rf, type="prob", newdata=xy_test)[,2]
        rf_pr_train <- prediction(rf_p_train, as.factor(xy_test$Present))
        auc[i] <- performance(rf_pr_train, measure = "auc")@y.values[[1]]
    }
    
    # Data frame contenant les statistiques de nos modeles
    model_stats = data.frame(
        Species_Name=species_name, 
        Species_Count=species_count,
        Final_Count=count(species_data_final_trimmed@data)$n,
        Train_Accuracy=mean(train_accuracy), 
        Train_Kappa=mean(train_kappa), 
        Train_NIR=mean(train_no_information_rate),
        Train_Auc=mean(train_auc), 
        Train_F1_Measure=mean(train_f1_measure),
        Accuracy=mean(accuracy), 
        Kappa=mean(kappa), 
        No_Information_Rate=mean(no_information_rate),
        Auc=mean(auc),
        F1_Measure=mean(f1_measure)
    )
                                      
    return(model_stats)
}
                                                                
# Pour chaque espèce...
species <- occurrence_data %>% count(Essence_latin, sort=TRUE)
models_statistics <- NA
for (i in 1:num_species){
    # Extraire les informations et les données de l'espèce
    species_name <- species$Essence_latin[i]
    species_count <- species$n[i]
    species_data <- (occurrence_data %>% dplyr::filter(Essence_latin==species_name) %>% ungroup())[2:4] 
    
    # Entrainer les modèles et sauvegarder le résultat
    trained_model_stats <- train_sdm(species_data, species_name, species_count)
    if (is.na(models_statistics)){
        models_statistics <- trained_model_stats
    } else {
        models_statistics <- rbind(models_statistics, trained_model_stats)
    }
}

# Sauvegarder le résultat des tests dans un fichier csv
write_csv(models_statistics, "species_models_statistics.csv")