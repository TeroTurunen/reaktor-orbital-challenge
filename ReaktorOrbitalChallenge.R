######################################################################
# Reaktor orbital challenge https://reaktor.com/orbital-challenge/
# Submitter:  tero.turunen@gmail.com
# LinkedIn:   https://fi.linkedin.com/in/teroturunen
######################################################################

# Load required libraries
if(!require(lpSolveAPI)){
  install.packages("lpSolveAPI")
}
library(lpSolveAPI)

# Constants 
EARTH_RADIUS <- 6371
ROUTE_START <- 1

# Distance calculation in km between two points on earth
distance <- function(long1, lat1, long2, lat2) {
  long1 <- (long1 * pi / 180)
  long2 <- (long2 * pi / 180)
  lat1 <- (lat1 * pi / 180)
  lat2 <- (lat2 * pi / 180)

  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = EARTH_RADIUS * c

  return(d)
}

# Line of sight calculation in km on earth 
lineofsight <- function(altitude){
  return (EARTH_RADIUS * acos(EARTH_RADIUS / (EARTH_RADIUS + altitude)))
}

# Load challenge file
file = 'C:/Users/mutte/Documents/Projects/Reaktor/reaktor_orbital_challenge.txt'
data <- read.csv(file , stringsAsFactors = FALSE, skip = 0, header=FALSE, col.names = c('ID','latitude','longitude','altitude')) 

# Extract seed
seed <- data[1, 1]
data <- data[-1,]

# Extract route information from dataset
rows <- nrow(data)
route_start_lat <- data[rows - 1, 2]
route_start_long <-  data[rows - 1, 3]
route_end_lat <- data[rows - 1, 4]
route_end_long <- as.numeric(data[rows, 1])

# Remove route information from dataset
data <- data[1:(rows - 2), ]

# Bind route back to dataset with altitude 0
data <- rbind(c('START', route_start_lat, route_start_long, 0 ), data)
data <- rbind(data, c('END', route_end_lat, route_end_long, 0 ) )

# Extract altitude and satelliteid as vectors and remove from dataset
satID <- data$ID
alt <- as.numeric(data$altitude)
data <- data[, c('latitude','longitude')]

# Create matrix of all index combinations 
satmatrix <- merge(1:rows, 1:rows, all=TRUE)

# Create distance vector between each satellite 
dist <- vector()
for (s in 1:nrow(satmatrix)){
    sat_dist <- distance(as.numeric(data$longitude[satmatrix[s,1]]), as.numeric(data$latitude[satmatrix[s,1]]), as.numeric(data$longitude[satmatrix[s,2]]), as.numeric(data$latitude[satmatrix[s,2]]) )
    dist <- append(dist, sat_dist)
}
# Create maximum line of sight vector between each satellite 
los <- vector()
for (s in 1:nrow(satmatrix)){
  sat_los <- lineofsight(alt[satmatrix[s, 1]]) + lineofsight(alt[satmatrix[s, 2]])
  los <- append(los, sat_los)
}

# Bind distance vector to satmatrix
satmatrix <- cbind(satmatrix, dist)
# Bind line of sight vector to satmatrix
satmatrix <- cbind(satmatrix, los)
# Remove self references
satmatrix <- satmatrix[!(satmatrix[,1] == satmatrix[,2]),] 
# Remove combinations not in line of sight
satmatrix <- satmatrix[satmatrix$dist <= satmatrix$los,]
# Reset row names
rownames(satmatrix) <- 1:nrow(satmatrix)
# Rename columns
names(satmatrix)[names(satmatrix)=='x'] <- 'sat_from'
names(satmatrix)[names(satmatrix)=='y'] <- 'sat_to'

#####################################
# Create mixed integer linear program to solve minimum hops and sequence      
#####################################

ROUTE_END <- rows

# Create cost vector. Hop cost is 1 and sequence 0
cost <- rep(1, nrow(satmatrix))
cost <- append(cost, rep(0, rows))

# Create model object
lpmodel <- make.lp(0, length(cost), verbose = "neutral")

# Set objective function to minimize hops and limit solver timeout to 60 seconds
set.objfn(lpmodel, cost)
lp.control(lpmodel, sense='min')
lp.control(lpmodel, timeout=60)

# Set model constraints
for (c in 1:ROUTE_END){
  if(c==ROUTE_START){
    # Do not accept incoming connections
    indc <- as.numeric(rownames(satmatrix[satmatrix$sat_to==c,]))
    costs <- rep(1, length(indc))
    add.constraint(lpmodel, costs, "=", 0, indc)
 
    # Set route start constraint, must start from here   
    indc <- as.numeric(rownames(satmatrix[satmatrix$sat_from==c,]))
    costs <- rep(1, length(indc))
    add.constraint(lpmodel, costs, "=", 1, indc)
  }
  else if(c==ROUTE_END){
    # Must have one incoming connection 
    indc <- as.numeric(rownames(satmatrix[satmatrix$sat_to==c,]))
    costs <- rep(1, length(indc))
    add.constraint(lpmodel, costs, "=", 1, indc)

    # No output connection allowed
    indc <- as.numeric(rownames(satmatrix[satmatrix$sat_from==c,]))
    costs <- rep(1, length(indc))
    add.constraint(lpmodel, costs, "=", 0, indc)
  }
  else{
    # Satellite must have both in and out link or none 
    indc <- as.numeric(rownames(satmatrix[satmatrix$sat_from==c,]))
    costs <- rep(1, length(indc))

    indc_to <- as.numeric(rownames(satmatrix[satmatrix$sat_to==c,])) 
    indc <- append(indc, indc_to)
    costs <- append(costs, rep(-1, length(indc_to)))
    add.constraint(lpmodel, costs, "=", 0, indc)
  }
}

# Create sequence vector for satellites and ground stations
satseq <- (nrow(satmatrix)+1):length(cost)

# Set sequence contraints
for(s in satseq){
  add.constraint(lpmodel, 1, "<=", rows, s)
  add.constraint(lpmodel, 1, ">=", 0, s)
}

# Subtour elimination (MTZ constraint)
for (n in 1:nrow(satmatrix)){
  if((satmatrix[n,]$sat_to > 1) & (satmatrix[n,]$sat_from > 1)){
  indc <- c(n, satseq[satmatrix[n,]$sat_to], satseq[satmatrix[n,]$sat_from] )
    costs <- c(rows, -1, 1)
    add.constraint(lpmodel, costs, "<=", (rows -1), indc)
  }
}

# Set satellite hop decision variables as binary
for (n in 1:(nrow(satmatrix))){
  set.type(lpmodel, n, c("binary"))
}

# Set sequence decision variables as integer
for (n in nrow(satmatrix):(nrow(satmatrix) + rows) ){
  set.type(lpmodel, n, c("integer"))
}

# Solve model
solve(lpmodel)

# Create result for satellite decision variables 
sat_result <- cbind.data.frame(
  sapply(satmatrix$sat_to, as.numeric),
  sapply(get.variables(lpmodel)[1:nrow(satmatrix)], as.numeric)
)
colnames(sat_result) <- c('satIndex', 'decision')

# Create result for sequence variables
seq_result <- cbind.data.frame(
  sapply(satID, as.character),
  sapply(get.variables(lpmodel)[satseq], as.numeric),
  stringsAsFactors = FALSE
)
colnames(seq_result) <- c('satID', 'satSeq')

# Select only satellites with incoming link and order resultset by sequence
seq_result <- seq_result[satID[sat_result[sat_result$decision==1,]$satIndex], ]
seq_result <- seq_result[order(seq_result$satSeq), ]

#####################################
# Print results
#####################################
# Print seed
# cat(seed)
# Print minimum hops including ground stations
#cat(get.objective(lpmodel))

# Print response as comma separated list
cat(paste(seq_result[seq_result$satID != 'END', c('satID') ], collapse=", "))


