##Overlays and extracts the abundance (freq of cell occurrence) of each pairs

merge.names<-function(hh){
  aa<-length(hh)
  bb<-rbind(as.data.frame(hh[[1]]),as.data.frame(hh[[2]]))
  if(aa==2){
  	return(bb)
  } else for(i in 2:aa){bb<-rbind(bb,as.data.frame(hh[[i]]))};
  return(bb)
}

species.search<-function(spec1,spec2){
#### Loading useful packages####
  my_packages<-c('raster', 'rgdal', 'biomod2', 'sp', 'maptools', 'SDMTools', 'rgbif', 'dismo', 'rgeos', 'XML')
  lapply(my_packages, require, character.only=T)

#   source("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals&Meetings/SIGC2013/sigccode/bs.functions.new.r")
#   source("c:/Users/sean maher/Dropbox/SIGC Resources Dropbox/Effects of species interactions on plants meta-analyses/bs.functions.new.r")

  worldclim_data <- getData('worldclim', var='bio', res=10)
  
  ## A list of length equal to length(species_names)
  species_records <- vector("list", length(spec1))
  ### Downloading records from GBIF using gbif() and storing them in species_records
  
  ## Criteria:
  ## Only georeferenced records
  ## Only records from 1 Jan 1970 onwards
  for (i in seq(along=species_records)){
  	species_records[[i]] <- gbif(unlist(lapply(strsplit(spec1, " "),
  	                                               function(x) x[1]))[i], 
                                paste(unlist(lapply(strsplit(spec1, " "),
                                                   function(x) x[2]))[i], "*", sep = ""), 
                                geo=TRUE, removeZeros = TRUE, args=c('startdate=1970-01-01'), nrecs = 10000, end = 10000)#, nrecs = 500, end = 500

  ## Data cleaning
  ## Removing records with > 20000 m
    species_records[[i]] <- species_records[[i]][species_records[[i]]$coordUncertaintyM <= 20000 | is.na(species_records[[i]]$coordUncertaintyM),]
    ## Adding a cellID field identifying which worldclim raster cell each gbif record falls into
    species_records[[i]]$cellID <- factor(cellFromXY(worldclim_data[[1]], species_records[[i]][c("lon", "lat")])) 
    ## Remove duplicated records that fall within the same cellID
    species_records[[i]] <- species_records[[i]][!duplicated(species_records[[i]]$cellID),]
    }
    
    if(length(species_records)==1){SP1 <- species_records;
    	                             SP1P <- as.data.frame(SP1)
} else SP1P <- merge.names(species_records) #for species with two species names (i.e. species name has changed)
#### Species Distribution 2 Data #### 
  ### Creating useful objects
    
  ## A list of length equal to length(species_names)
  species_records2 <- vector("list", length(spec2))
  ### Downloading records from GBIF using gbif() and storing them in species_records
  ## Criteria:
  ## Only georeferenced records
  ## Only records from 1 Jan 1970 onwards
  for (i in seq(along=species_records2)){
  	species_records2[[i]] <- gbif(unlist(lapply(strsplit(spec2, " "), function(x) x[1]))[i],
  	  paste(unlist(lapply(strsplit(spec2, " "), function(x) x[2]))[i], "*", sep = ""),
  	  geo=TRUE, removeZeros = TRUE, args=c('startdate=1970-01-01'), nrecs = 10000, end = 10000)#, nrecs = 500, end=500
  ### Data cleaning
  ## Removing records with > 20000 m uncertainty
  species_records2[[i]] <- species_records2[[i]][species_records2[[i]]$coordUncertaintyM <= 20000 | is.na(species_records2[[i]]$coordUncertaintyM),]
  ## Adding a cellID field identifying which worldclim raster cell each gbif record falls into
  species_records2[[i]]$cellID <- factor(cellFromXY(worldclim_data[[1]], species_records2[[i]][c("lon", "lat")])) 
  ## Remove duplicated records that fall within the same cellID
     species_records2[[i]] <- species_records2[[i]][!duplicated(species_records2[[i]]$cellID),]
     }
     
     if(length(species_records2)==1){SP2 <- species_records2;
     	                             SP2P <- as.data.frame(SP2)
} else SP2P <- merge.names(species_records2) #for species with two species names (i.e. species name has changed)

  
  #take out NAs
  ind<-which((is.na(SP1P$lon)==F)&(is.na(SP1P$lat)==F))
  SP1P<-SP1P[ind,]
  ind<-which((is.na(SP2P$lon)==F)&(is.na(SP2P$lat)==F))
  SP2P<-SP2P[ind,]
  
  ####################
  ####PROCESSING######
  ####################
  
  #plot points
  data(wrld_simpl)
  #   plot(wrld_simpl, xlim=c(-130,90), ylim=c(-50,70),axes=TRUE, col='light yellow')
  #   points(SP1P$lon, SP1P$lat, col='orange', pch=20, cex=0.75)
  #   points(SP2P$lon, SP2P$lat, col='blue', pch=20, cex=0.75)
  #   
  #make a 40km radius bubble around points for "species range", not sure how else to do this....
   
   proJ<-CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
   coordinates(SP1P) <- ~lon+lat
   coordinates(SP2P) <- ~lon+lat
   proj4string(SP1P) <- proJproj4string(SP2P) <- proJ
   circsSP1 <- circles(SP1P, d=40000, lonlat=TRUE)
   polSP1 <- gUnaryUnion(circsSP1@polygons)
   #   proj4string(polSP1) <- proj
   circsSP2 <- circles(SP2P, d=40000, lonlat=TRUE)
   polSP2 <- gUnaryUnion(circsSP2@polygons)
   
   #crop the bubbles to continents not oceans: 
   rangeSP1 <- crop(x=polSP1,wrld_simpl)
   rangeSP2 <- crop(x=polSP2,wrld_simpl)
   rangeSP12<-gUnion(rangeSP1, rangeSP2)

#   plot(wrld_simpl, xlim=c(-130,90), ylim=c(-50,70),axes=TRUE, col='light yellow')
#   plot(rangeSP1,col=5,add=T)
#   plot(rangeSP2,col=4,add=T)
#   plot(rangeSP12,col=3,add=T)
# 
  #crop predictors to the extent of the ranges of sp 1+sp 2
  predictors.crop <- mask(x=worldclim_data,rangeSP12)
  #   plot(predictors.crop$bio1)
  
  #get lat long for all points in both species ranges
  latlons<-rasterToPoints(predictors.crop,spatial=T)
  cells<-cellFromXY(predictors.crop, latlons)
  allcellsAB<-cbind(cells,rasterToPoints(predictors.crop,spatial=F))
  #   head(allcellsAB)

  #get cell numbers for AB, A0, 0B, 00
    spApres<-extract(predictors.crop,SP1P,cellnumbers=T)
    spBpres<-extract(predictors.crop,SP2P,cellnumbers=T)
    
    cAB<-intersect(spApres[,"cells"],spBpres[,"cells"])
    cA0<-setdiff(spApres[,"cells"],spBpres[,"cells"])
    c0B<-setdiff(spBpres[,"cells"],spApres[,"cells"])
    c00<-setdiff(allcellsAB[,"cells"],union(spApres[,"cells"],spBpres[,"cells"]))
    
   #get data for all those cell numbers
   pAB<-cbind(extract(predictors.crop,cAB),cAB,xyFromCell(predictors.crop,cAB))
   pA0<-cbind(extract(predictors.crop,cA0),cA0,xyFromCell(predictors.crop,cA0))
   p0B<-cbind(extract(predictors.crop,c0B),c0B,xyFromCell(predictors.crop,c0B))
   p00<-cbind(extract(predictors.crop,c00),c00,xyFromCell(predictors.crop,c00))
   
   return(list(pAB=pAB,
              pA0=pA0,
              p0B=p0B,
              p00=p00))
              }

species.check<-function(spec1,spec2){
	aaa<-species.search(spec1,spec2)
	return(c(dim(aaa[[1]])[1],
	      dim(aaa[[2]])[1],
	      dim(aaa[[3]])[1],
	      dim(aaa[[4]])[1]))
	}
	
	