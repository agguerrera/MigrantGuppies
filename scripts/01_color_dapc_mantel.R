# Load required packages
library(psych)
library(stats)
library(adegenet)
library(factoextra)
library(vegan)
library(dplyr)
library(tidyr)
library(purrr)
library(here)

# Necessary input files:
# rgb_uv_pool5DT1PR.rds
# xy5DT1PR.csv

#Read in color data and xy coordinates
rgb_uv_pool5DT1PR <- readRDS(here::here("data", "rgb_uv_pool5DT1PR.rds"))
xy5DT1PR <- read.csv(here("data", "xy5DT1PR.csv"))

############################################
# DAPC for individual pools within streams #
############################################
myCols <- c("red4","red3","red1","darkorange3","darkorange","orange","goldenrod2","yellow3","yellow","navy","blue","blue3","cornflowerblue","steelblue","cyan4","cyan3","cyan","lightblue","purple4","orchid4","plum4","plum3","orchid","rosybrown3","rosybrown2","pink4","pink2")

rgb_uv_pool5DT1PR$StrPool <- paste0(substr(rgb_uv_pool5DT1PR$Stream, 1, 1), rgb_uv_pool5DT1PR$PoolLetter)
rgb_uv_pool5DT1PR <- rgb_uv_pool5DT1PR %>% relocate(StrPool, .after = 1)

dapc_StrPool = dapc(rgb_uv_pool5DT1PR[,c(10:28833)], rgb_uv_pool5DT1PR$StrPool, var.contrib = TRUE, scale = FALSE, n.pca = 100, n.da = 24, var.loadings = TRUE) #n.da is the number of discriminant functions to use, which defaults to 1 less than number of groups
scatter(dapc_StrPool, clab=0,scree.da=TRUE, main = "DAPC by Pool", col=myCols, bg="white",leg=FALSE, cex=1, lwd=1)

# Calculating variance explained by each LD
eigenvalues <- dapc_StrPool$eig
total_variance <- sum(eigenvalues)
var_explained <- eigenvalues / total_variance

dapc_StrPool$ind.coord
axis(1, at = c(-8,-4,4,8) ,labels = TRUE, tick = TRUE, pos = 0, mgp=c(3,0.75,0))
axis(2, at = c(-8,-4,4,8),labels = TRUE, tick = TRUE, pos = 0, las=2, mgp=c(3,0.75,0))

myInset <- function(){
  temp <- dapc_StrPool$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc_StrPool$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated variance (%)",
       cex=1, pch=20, type="h", lwd=2)
}
add.scatter(
  myInset(),
  posi = "topright",
  inset = c(0.02, 0.02),
  ratio = 0.15,
  bg = transp("white")
)


#######################
# PERMANOVA analysis #
#######################
# PERMANOVA across all pools using first 100 PCs on raw r,g,b
pools.PCA <- prcomp(rgb_uv_pool5DT1PR[,c(10:28833)], scale.=FALSE)
fviz_eig(pools.PCA, addlabels = TRUE)
pools.varexp <- pools.PCA$sdev^2 / sum(pools.PCA$sdev^2)
sum(pools.varexp[1:100])
pools.PCs <- data.frame(pools.PCA$x[,1:100])
pools.PCs$Pools <- rgb_uv_pool5DT1PR$StrPool
pools.PCs <- pools.PCs %>% # Move distance column to the end 
  relocate(Pools, .before = 1)  


#PERMANOVA across groups using first 100 PCs on raw r,g,b
pools.PC.df <- data.frame(pools.PCs[,2:101])
permanova <- adonis2(pools.PC.df ~ Pools, data = pools.PCs, permutations=999, method = "euc")


###############################################################
# Pairwise mean pool color pattern versus geographic distance # 
###############################################################
#Are pools within streams (tributaries) significantly different from each other?

#Aripo Intro 
# dropping Pools F & G with < 5 males)
Aripo_df <- rgb_uv_pool5DT1PR %>%
  filter(Stream == "Aripo",
         !PoolLetter %in% c("F", "G")) %>%
  droplevels()

#PERMANOVA
# Color variables only
Aripo.colors <- Aripo_df[, 10:28833]

# PERMANOVA
Ari.permanova <- adonis2(
  Aripo.colors ~ PoolLetter,
  data = Aripo_df,
  permutations = 9999,
  method = "euclidean")

Ari.permanova

# Pairwise PERMANOVA
Ari.dist <- dist(Aripo.colors, method = "euclidean")

Ari.ppm <- RVAideMemoire::pairwise.perm.manova(
  Ari.dist,
  Aripo_df$PoolLetter,
  nperm = 9999)

Ari.ppm

# Mantel Test
# Here, I am calculating the column means for all color channels across all sampling points 
Aripo_means <- Aripo_df %>% #Aripo_df is my dataframe with all my color values and my grouping variable I care about (PoolLetter)
  group_by(PoolLetter) %>%
  summarise(across(starts_with("r_") | starts_with("g_") | starts_with("b_") | starts_with("uv_"), ~mean(., na.rm = TRUE)),
            .groups = "drop") #Now I have a mean color vector for each of my groups
Ari.means <- Aripo_means[,c(2:28825)]
Ari.means <- as.matrix(Ari.means)
Ari.meansDist <- dist(Ari.means, method = "euclidean")
Ari.meansDist <- as.matrix(Ari.meansDist)
rownames(Ari.meansDist) <- c("A","B","C","D","E","H","I")
colnames(Ari.meansDist) <- c("A","B","C","D","E","H","I")

# Use lower.tri to extract the lower triangle values
AriDist_lt <- Ari.meansDist[lower.tri(Ari.meansDist)]
AriDist_long <- c()  # Empty vector to store comparison names
# Loop through to generate comparison names for each pair in the lower triangle
for (i in 1:(nrow(Ari.meansDist)-1)) {
  for (j in (i+1):ncol(Ari.meansDist)) {
    # Combine the names of the comparisons
    AriDist_long <- c(AriDist_long, paste(rownames(Ari.meansDist)[i], colnames(Ari.meansDist)[j], sep = "_"))
  }
}

# Check the results
AriDist_long  # Comparison names for the lower triangle
AriDist_lt  # Distance values for the lower triangle

AriComp <-data.frame(cbind(AriDist_long,AriDist_lt))
colnames(AriComp) <- c("PoolPair","MeanDist")

Aripo_means$distance <- NA  # Initialize with NA values
# Apply conditions and assign values based on PoolLetter
Aripo_means$distance[Aripo_means$PoolLetter == "A"] <- 5.2
Aripo_means$distance[Aripo_means$PoolLetter == "B"] <- 51.85
Aripo_means$distance[Aripo_means$PoolLetter == "C"] <- 84.75
Aripo_means$distance[Aripo_means$PoolLetter == "D"] <- 125.55
Aripo_means$distance[Aripo_means$PoolLetter == "E"] <- 179.4
Aripo_means$distance[Aripo_means$PoolLetter == "H"] <- 482.75
Aripo_means$distance[Aripo_means$PoolLetter == "I"] <- 507.5

Aripo_means <- Aripo_means %>% # Move distance column to the end temporarily
  relocate(distance, .after = 1)  # Move distance to the 9th position

Ari.gdist <- Aripo_means[, 2, drop = FALSE]
Ari.gdist <- dist(Ari.gdist,method = "euclidean")
Ari.gdist.vec <- as.vector(Ari.gdist)

Ari.dist <- data.frame(Ari.gdist.vec,AriComp)
Ari.dist$MeanDist <- as.numeric(Ari.dist$MeanDist)

Fig2.1 <- ggplot(Ari.dist, aes(x = Ari.gdist.vec, y = MeanDist)) +
  geom_point() + geom_smooth(method= "lm", col="black", se=FALSE) +
  labs(x = "Pairwise Geographic Distance", y = "Pairwise Differences in Mean Color Pattern") +
  theme_classic() + scale_x_continuous(limits = c(0, 600)) +
  scale_y_continuous(limits = c(0, 8.5)) 

ggsave(filename = here("figures", "Fig2.1.png"),
       plot = Fig2.1,
       width = 8,
       height = 6,
       dpi = 300)

Aripo.mantel <- mantel(Ari.gdist, Ari.meansDist,method = "pearson", permutations = 999)
Aripo.mantel

#Naranjo
Nar_df <- subset(rgb_uv_pool5DT1PR, Stream=="Nar")

#PERMANOVA
# Color variables only
Nar.colors <- Nar_df[, 10:28833]

# PERMANOVA
Nar.permanova <- adonis2(
  Nar.colors ~ PoolLetter,
  data = Nar_df,
  permutations = 9999,
  method = "euclidean"
)

# Pairwise PERMANOVA
Nar.dist <- dist(Nar.colors, method = "euclidean")

Nar.ppm <- RVAideMemoire::pairwise.perm.manova(
  Nar.dist,
  Nar_df$PoolLetter,
  nperm = 9999
)

Nar.ppm


#Mantel
Nar_means <- Nar_df %>%
  group_by(PoolLetter) %>%
  summarise(across(starts_with("r_") | starts_with("g_") | starts_with("b_") | starts_with("uv"), ~mean(., na.rm = TRUE)),
            .groups = "drop")
Nar.means <- Nar_means[,c(2:28825)]
Nar.means <- as.matrix(Nar.means)
Nar.meansDist <- dist(Nar.means, method = "euclidean")
Nar.meansDist <- as.matrix(Nar.meansDist)
rownames(Nar.meansDist) <- c("A","B","C","D","E","F","G","H","I")
colnames(Nar.meansDist) <- c("A","B","C","D","E","F","G","H","I")

# Use lower.tri to extract the lower triangle values
NarDist_lt <- Nar.meansDist[lower.tri(Nar.meansDist)]

NarDist_long <- c()  # Empty vector to store comparison names
# Loop through to generate comparison names for each pair in the lower triangle
for (i in 1:(nrow(Nar.meansDist)-1)) {
  for (j in (i+1):ncol(Nar.meansDist)) {
    # Combine the names of the comparisons
    NarDist_long <- c(NarDist_long, paste(rownames(Nar.meansDist)[i], colnames(Nar.meansDist)[j], sep = "_"))
  }
}

# Check the results
NarDist_long  # Comparison names for the lower triangle
NarDist_lt  # Distance values for the lower triangle

NarComp <-data.frame(cbind(NarDist_long,NarDist_lt))
colnames(NarComp) <- c("PoolPair","MeanDist")


Nar_means$distance <- NA  # Initialize with NA values
# Apply conditions and assign values based on PoolLetter
Nar_means$distance[Nar_means$PoolLetter == "A"] <- 12.3
Nar_means$distance[Nar_means$PoolLetter == "B"] <- 61.25
Nar_means$distance[Nar_means$PoolLetter == "C"] <- 118.25
Nar_means$distance[Nar_means$PoolLetter == "D"] <- 217.5
Nar_means$distance[Nar_means$PoolLetter == "E"] <- 257.5
Nar_means$distance[Nar_means$PoolLetter == "F"] <- 296.1
Nar_means$distance[Nar_means$PoolLetter == "G"] <- 321.4
Nar_means$distance[Nar_means$PoolLetter == "H"] <- 343.5
Nar_means$distance[Nar_means$PoolLetter == "I"] <- 382.4

Nar_means <- Nar_means %>% # Move distance column to the end temporarily
  relocate(distance, .after = 1)  # Move distance to the 9th position

Nar.gdist <- Nar_means[, 2, drop = FALSE]
Nar.gdist <- dist(Nar.gdist,method = "euclidean")
Nar.gdist.vec <- as.vector(Nar.gdist)

Nar.dist <- data.frame(Nar.gdist.vec,NarComp)
Nar.dist$MeanDist <- as.numeric(Nar.dist$MeanDist)

Fig2.2 <- ggplot(Nar.dist, aes(x = Nar.gdist.vec, y = MeanDist)) +
  geom_point() + geom_smooth(method= "lm", col="black", se=FALSE) +
  labs(x = "Pairwise Geographic Distance", y = "Pairwise Differences in Mean Color Pattern") +
  theme_classic() + scale_x_continuous(limits = c(0, 600)) +
  scale_y_continuous(limits = c(0, 8.5)) 

ggsave(filename = here("figures", "Fig2.2.png"),
       plot = Fig2.2,
       width = 8,
       height = 6,
       dpi = 300)
#Mantel test
Nar.mantel <- mantel(Nar.gdist, Nar.meansDist,method = "pearson", permutations = 999)
Nar.mantel



#Mainstream
Main_df <- subset(rgb_uv_pool5DT1PR, Stream=="Temp")

#PERMANOVA
# Color variables only
Main.colors <- Main_df[, 10:28833]

# PERMANOVA
Main.permanova <- adonis2(
  Main.colors ~ PoolLetter,
  data = Main_df,
  permutations = 9999,
  method = "euclidean"
)

# Pairwise PERMANOVA
Main.dist <- dist(Main.colors, method = "euclidean")

Main.ppm <- RVAideMemoire::pairwise.perm.manova(
  Main.dist,
  Main_df$PoolLetter,
  nperm = 9999
)

Main.ppm

#Mantel
Main_means <- Main_df %>%
  group_by(PoolLetter) %>%
  summarise(across(starts_with("r_") | starts_with("g_") | starts_with("b_") | starts_with("uv_"), ~mean(., na.rm = TRUE)),
            .groups = "drop")
Main.means <- Main_means[,c(2:28825)]
Main.means <- as.matrix(Main.means)
Main.meansDist <- dist(Main.means, method = "euclidean")
Main.meansDist <- as.matrix(Main.meansDist)
rownames(Main.meansDist) <- c("A","B","C","D","E","F","G","H","I")
colnames(Main.meansDist) <- c("A","B","C","D","E","F","G","H","I")

# Use lower.tri to extract the lower triangle values
MainDist_lt <- Main.meansDist[lower.tri(Main.meansDist)]

MainDist_long <- c()  # Empty vector to store comparison names
# Loop through to generate comparison names for each pair in the lower triangle
for (i in 1:(nrow(Main.meansDist)-1)) {
  for (j in (i+1):ncol(Main.meansDist)) {
    # Combine the names of the comparisons
    MainDist_long <- c(MainDist_long, paste(rownames(Main.meansDist)[i], colnames(Main.meansDist)[j], sep = "_"))
  }
}

# Check the results
MainDist_long  # Comparison names for the lower triangle
MainDist_lt  # Distance values for the lower triangle

MainComp <-data.frame(cbind(MainDist_long,MainDist_lt))
colnames(MainComp) <- c("PoolPair","MeanDist")


Main_means$distance <- NA  # Initialize with NA values

# Apply conditions and assign values based on PoolLetter
Main_means$distance[Main_means$PoolLetter == "A"] <- 22.7
Main_means$distance[Main_means$PoolLetter == "B"] <- 60.8
Main_means$distance[Main_means$PoolLetter == "C"] <- 100.05
Main_means$distance[Main_means$PoolLetter == "D"] <- 128.9
Main_means$distance[Main_means$PoolLetter == "E"] <- 189.4
Main_means$distance[Main_means$PoolLetter == "F"] <- 262.95
Main_means$distance[Main_means$PoolLetter == "G"] <- 340.05
Main_means$distance[Main_means$PoolLetter == "H"] <- 378.3
Main_means$distance[Main_means$PoolLetter == "I"] <- 467.85

Main_means <- Main_means %>% # Move distance column to the end Mainorarily
  relocate(distance, .after = 1)  # Move distance to the 9th position

Main.gdist <- Main_means[, 2, drop = FALSE]
Main.gdist <- dist(Main.gdist,method = "euclidean")

Main.gdist.vec <- as.vector(Main.gdist)

Main.dist <- data.frame(Main.gdist.vec,MainComp)
Main.dist$MeanDist <- as.numeric(Main.dist$MeanDist)

Fig2.3 <- ggplot(Main.dist, aes(x = Main.gdist.vec, y = MeanDist)) +
  geom_point() + geom_smooth(method= "lm", col="black", se=FALSE) +
  labs(x = "Pairwise Geographic Distance", y = "Pairwise Differences in Mean Color Pattern") +
  theme_classic() + scale_x_continuous(limits = c(0, 600)) +
  scale_y_continuous(limits = c(0, 8.5)) 

ggsave(filename = here("figures", "Fig2.3.png"),
       plot = Fig2.3,
       width = 8,
       height = 6,
       dpi = 300)
#Mantel test
Main.mantel <- mantel(Main.gdist, Main.meansDist,method = "pearson", permutations = 999)



###################################################################
### Pulling out orange, black, and uv reflectance for each male ###
###################################################################
rgb_uv_comp <- rgb_uv_pool5DT1PR %>%
  pivot_longer(
    cols = matches("^(r|g|b|uv)_"),
    names_to = c("channel", "point"),
    names_pattern = "^(r|g|b|uv)_(.+)$"
  ) %>%
  pivot_wider(names_from = channel, values_from = value)

# Function to convert RGB to HSV
rgb_to_hsv <- function(r, g, b) {
  # Find min and max of RGB values
  cmax <- max(r, g, b)
  cmin <- min(r, g, b)
  delta <- cmax - cmin
  # Calculate Hue (H)
  if (delta == 0) {
    h <- 0  # Undefined, but we can set to 0
  } else if (cmax == r) {
    h <- (60 * ((g - b) / delta) + 360) %% 360
  } else if (cmax == g) {
    h <- (60 * ((b - r) / delta) + 120) %% 360
  } else {
    h <- (60 * ((r - g) / delta) + 240) %% 360
  }
  # Calculate Saturation (S)
  s <- if (cmax == 0) {
    0
  } else {
    delta / cmax
  }
  # Calculate Value (V)
  v <- cmax
  # Return the HSV values as a named vector
  return(c(H = h, S = s, V = v))
}


# Apply rgb_to_hsv function row-wise and add columns for H, S, and V
fish_hue <- rgb_uv_comp %>%
  mutate(
    hsv_list = pmap(list(r, g, b), rgb_to_hsv),
    hue_value = map_dbl(hsv_list, ~ .x["H"]),
    saturation = map_dbl(hsv_list, ~ .x["S"]),
    value = map_dbl(hsv_list, ~ .x["V"])
  )


# Make sure hsv values are in the range [0,1]
fish_hue <- fish_hue %>%
  mutate(
    hsv_color = hsv(hue_value / 360, saturation, value)
  )

# Apply rgb_to_hsv function row-wise and add columns for H, S, and V
fish_hue$hue <- fish_hue$hue_value/360

hue_xy <- merge(fish_hue, xy5DT1PR, by = "point", all.x = TRUE)

hue_body <- hue_xy %>%
  filter(!(x >= 788 & x <= 880 & y >= 358 & y <= 430))

#Create orange threshold
orange_pixels <- hue_body %>%
  filter(hue >= 0, hue <= 0.1875,
         saturation >= 0.2,
         value >= 0.4)

#Calculate total number of pixels meeting orange threshold
orange_counts <- orange_pixels %>%
  group_by(ID) %>%
  summarise(orange_area = n())

#Compute total area of orange on each male
orange_counts$totalpts <- 6151
orange_counts$TotalOr <- orange_counts$orange_area/orange_counts$totalpts

#Create black threhold
black_pixels <- hue_body %>%
  filter(hue >= 0,
         saturation >= 0.08,
         value >= 0, value <= .3)

#Calculate total number of pixels meeting black threhold
black_counts <- black_pixels %>%
  group_by(ID) %>%
  summarise(black_area = n())

#Compute total black area on each male
black_counts$totalpts <- 6151
black_counts$TotalBlk <- black_counts$black_area/black_counts$totalpts


#Create summed value of UV
fish_hue_summed <- hue_body %>%
  group_by(ID) %>%
  summarise(
    uv_sum = sum(uv, na.rm = TRUE)
  )

# Figure S8: How orange and black are being pulled out
fish_id <- c("TR32MUP0TEMP.tif")
fish1 <- hue_xy %>% filter(ID == fish_id)

fish.black <- black_pixels %>% filter(ID == fish_id)
fish.orange <- orange_pixels %>% filter(ID == fish_id)
fish1 <- hue_xy %>% filter(ID == fish_id)

FigS8 <- ggplot(fish.black, aes(x = x, y = y)) + 
  geom_point(data = fish1, aes(x = x, y = y), color = "grey") +
  geom_point(data = fish.orange, aes(x = x, y = y), color = "orange", alpha=0.5) +
  geom_point(data = fish.black, aes(x = x, y = y), color = "black") +
  theme_classic()

ggsave(filename = here("figures", "FigS8.png"),
       plot = FigS8,
       width = 8,
       height = 6,
       dpi = 300)
