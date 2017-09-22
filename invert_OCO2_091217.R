#-- This is necessary to read ncdf3/ncdf4-hdf5 files
library(ncdf4)

#- This is necessary for openmp use in R
library(foreach)
library(doSNOW)

#-  This does a generic bind of arrays
library(abind)

#-- Utility functions

grid2transcom = function(mat,model.grid.x=1,model.grid.y=1,transcom=TRUE)
  {
  dist.calc = function(lat1,lat2,lon1,lon2)
  {
    a = (sin((lat2-lat1)/90*0.5*pi/2))^2 + (cos(lat1/90*pi/2)*cos(lat2/90*pi/2))*
               (sin((lon2-lon1)/90*0.5*pi/2))^2
    c = 2*atan2(sqrt(a),sqrt(1-a))
    d = 6371 * c	
 	return(d)
  }

  if(model.grid.x==1 & model.grid.y==1){
   transfil = nc_open("/discover/nobackup/aschuh/data/misc/iregions.nc")

   regions = ncvar_get(transfil,"transcom_regions")

   nc_close(transfil)
   }else{
    transfil = nc_open("/discover/nobackup/aschuh/data/transcom/transcom_regions_2x2.5.nc")

    regions = ncvar_get(transfil,"region")

    nc_close(transfil)
  }

 #-- Areas for grid cells

  grid.y=model.grid.y
  grid.x=model.grid.x

  if(model.grid.x==1 & model.grid.y==1){
	tabl.transcom = cbind(c(seq(-90+0.5*grid.y,
                      90-0.5*grid.y,by=grid.y)),
                      c(seq(-90+0.5*grid.y,
                      90-0.5*grid.y,by=grid.y)),                   
                       rep(0,180),rep(grid.x,180))
  }else{
	tabl.transcom = cbind(c(-90+0.5*(0.5*grid.y),seq(-90+0.5*grid.y,
                      90-1.5*grid.y,by=grid.y)+0.5*grid.y,90-0.5*(0.5*grid.y)),
                      c(-90+0.5*(0.5*grid.y),seq(-90+0.5*grid.y,
                      90-1.5*grid.y,by=grid.y)+0.5*grid.y,90-0.5*(0.5*grid.y)),                   
                       rep(0,91),rep(grid.x,91))
                  } 
                  
  dist.lat = apply(tabl.transcom,1,FUN=function(X){ dist.calc(X[1],X[2],X[3],X[4]) } )                 
                 
  dist.lat.mat = matrix(rep(dist.lat,360),nrow=360,ncol=180,byrow=T)

  #sapply(1:56,FUN=function(x){sum(aggregate(as.vector(neetmp2[,,x] * 111 * dist.lat.mat * 1000^2),list(as.vector(regions)),sum)$x)})
  if(transcom)
  {
  transcom.agg      = aggregate(as.vector(mat * 111 * dist.lat.mat * 1000^2),list(as.vector(regions)),sum)
  transcom.agg.vect = as.vector(unlist(transcom.agg))
  return(transcom.agg.vect)
  }else{
  	return(mat * 111 * dist.lat.mat * 1000^2)
}
}

grid2ocean30 = function(mat,model.grid.x=1,model.grid.y=1)
  {
  dist.calc = function(lat1,lat2,lon1,lon2)
  {
    a = (sin((lat2-lat1)/90*0.5*pi/2))^2 + (cos(lat1/90*pi/2)*cos(lat2/90*pi/2))*
               (sin((lon2-lon1)/90*0.5*pi/2))^2
    c = 2*atan2(sqrt(a),sqrt(1-a))
    d = 6371 * c	
 	return(d)
  }

  if(model.grid.x==1 & model.grid.y==1){
   transfil = nc_open("/discover/nobackup/aschuh/data/misc/iregions.nc")

   regions = ncvar_get(transfil,"transcom_regions")

   nc_close(transfil)
   }else{
    transfil = nc_open("/discover/nobackup/aschuh/data/transcom/transcom_regions_2x2.5.nc")

    regions = ncvar_get(transfil,"region")
    o_regions = ncvar_get(transfil,"ocean_regions")

    nc_close(transfil)
  }

 #-- Areas for grid cells

  grid.y=model.grid.y
  grid.x=model.grid.x

  if(model.grid.x==1 & model.grid.y==1){
	tabl.transcom = cbind(c(seq(-90+0.5*grid.y,
                      90-0.5*grid.y,by=grid.y)),
                      c(seq(-90+0.5*grid.y,
                      90-0.5*grid.y,by=grid.y)),                   
                       rep(0,180),rep(grid.x,180))
  }else{
	tabl.transcom = cbind(c(-90+0.5*(0.5*grid.y),seq(-90+0.5*grid.y,
                      90-1.5*grid.y,by=grid.y)+0.5*grid.y,90-0.5*(0.5*grid.y)),
                      c(-90+0.5*(0.5*grid.y),seq(-90+0.5*grid.y,
                      90-1.5*grid.y,by=grid.y)+0.5*grid.y,90-0.5*(0.5*grid.y)),                   
                       rep(0,91),rep(grid.x,91))
                  } 
                  
  dist.lat = apply(tabl.transcom,1,FUN=function(X){ dist.calc(X[1],X[2],X[3],X[4]) } )                 
                 
  dist.lat.mat = matrix(rep(dist.lat,360),nrow=360,ncol=180,byrow=T)

  #sapply(1:56,FUN=function(x){sum(aggregate(as.vector(neetmp2[,,x] * 111 * dist.lat.mat * 1000^2),list(as.vector(regions)),sum)$x)})
  transcom.agg      = aggregate(as.vector(mat * 111 * dist.lat.mat * 1000^2),list(as.vector(o_regions)),sum)
  transcom.agg.vect = as.vector(unlist(transcom.agg))
  return(transcom.agg.vect)
}


pad = function (x, width, fill = " ", left = TRUE) 
{
    xneg = FALSE
    nc = nchar(x)
    if (width > nc) {
        if (left) {
            str = paste(paste(rep(fill, width - nc), collapse = ""), 
                x, sep = "")
        }
        else {
            str = paste(x, paste(rep(fill, width - nc), collapse = ""), 
                sep = "")
        }
    }
    else {
        str = x
    }
    return(str)
}

###############################################################
#--  Pull in OCO-2 files w/ ancillary variables we can use
#--  These should match observation by observation to output files
###############################################################  

#lite_loc = "/discover/nobackup/aschuh/run/lite_wbias_20160826_standardized/inv_weighted/"
#lite_loc = "/discover/nobackup/aschuh/run/lite_wbias_20161205_standardized/inv_weight_insitu/"
#lite_loc = "/discover/nobackup/aschuh/run/lite_wbias_20161208_standardized/inv_weight_insitu/"
lite_loc = "/discover/nobackup/aschuh/run/lite_wbias_20170815_standardized/inv_weight_insitu/"

  
fls.orig.2014 = sort(list.files(lite_loc,
     recursive=FALSE,pattern="2014.4x5",full.names=TRUE))
     
fls.orig.2015 = sort(list.files(lite_loc,
     recursive=FALSE,pattern="2015.4x5",full.names=TRUE))  

fls.orig.2016 = sort(list.files(lite_loc,
     recursive=FALSE,pattern="2016.4x5",full.names=TRUE))
          
fls.orig = c(fls.orig.2014,fls.orig.2015,fls.orig.2016)
     
for(i in 1:(length(fls.orig)))
#for(i in 1:528)
{
	cat(i)
	con = nc_open(fls.orig[i])
	type_in = ncvar_get(con,"type")
	median_in = ncvar_get(con,"xco2_median")
	time_in =  ncvar_get(con,"time")
	lat = ncvar_get(con,"latitude")
	lon =  ncvar_get(con,"longitude")
	logdws = ncvar_get(con,"xco2_logdws")
	co2_grad_del = ncvar_get(con,"xco2_co2_grad_del")	
	dp = ncvar_get(con,"xco2_dp")	
	s31 = ncvar_get(con,"xco2_s31")
	sza = ncvar_get(con,"xco2_sza")	
    airmass = ncvar_get(con,"xco2_airmass")		
    #balbedo = ncvar_get(con,"xco2_balbedo")
    xco2 = ncvar_get(con,"xco2")
    xco2_sd = ncvar_get(con,"xco2_uncertainty")
    xco2_raw = ncvar_get(con,"xco2_raw")
    
	nc_close(con)
	if(i==1)
	{
		type_orig = type_in
		xco2_median = median_in
		timev = time_in
		latv = lat
		lonv = lon
	   logdwsv = logdws
	   co2_grad_delv = co2_grad_del
    	dpv = dp
    	s31v = s31
    	szav = sza
       airmassv = airmass
       #balbedov = balbedo
      xco2_sdv = xco2_sd
      xco2v = xco2
      xco2_rawv = xco2_raw
	}else{
		type_orig=c(type_orig,type_in)
		xco2_median = c(xco2_median,median_in)
		timev = c(timev,time_in)
		latv = c(latv,lat)
		lonv = c(lonv,lon)		
		s31v=c(s31v,s31)
		szav = c(szav,sza)
	 logdwsv = c(logdwsv,logdws)
	   co2_grad_delv = c(co2_grad_delv,co2_grad_del)
    	dpv = c(dpv,dp)
    	airmassv = c(airmassv,airmass)
      #balbedov = c(balbedov,balbedo)    	
    	 xco2_sdv = c( xco2_sdv,xco2_sd)
    	xco2v = c(xco2v,xco2)
       xco2_rawv = c(xco2_rawv,xco2_raw)
    }
}

############################
#-  End of OCO-2 data pull
############################


##################################################################
#-- TCCON data was put into the "regular" xco2 vector (xco2v)
#-- This drops those data into the raw xco2 vector to repace NAs
##################################################################

xco2_rawv[type_orig==200] = xco2v[type_orig==200]

###########################################################################
#-- I'm putting a centering term when I create these 3 vars in original OCO2
#-- input files so I have to put this code in to
#-- uncenter bias terms, need to remove this in lite_plus_tccon*.R script
###########################################################################

logdwsv = logdwsv - 2.9

co2_grad_delv = co2_grad_delv +8.4

dpv = dpv +1.4

#s31v = s31v - 0.2

################################################################################
#-- INITIAL PULL OF OBSPACK DATA
#--  Output is coming out in concentrations and needs to come out in ppm hence ....
##################################################################################

scalar_noaa =10^6

#-- "Bad" indicator for obspack data to be assimilated

con = nc_open("/discover/nobackup/aschuh/run/obspack_co2_1_NRT_v3.3_prepped_inputs_2017-04-26_links/total_flask_input.nc")
mdm = ncvar_get(con,"mdm")
noaa_lat = ncvar_get(con,"latitude")
noaa_lon = ncvar_get(con,"longitude")
noaa_time = ncvar_get(con,"time")
realflask = ncvar_get(con,"value")*scalar_noaa

id = ncvar_get(con,"obspack_id")

#-- If mdm doesn't exist, we shouldn't assimilate
good_obspack = !is.na(mdm) 
nc_close(con)


#--  Join OCO2 and NOAA CO2 obs/errors

#co2v = c(xco2v,realflask)

#co2v_sd = c(xco2_sdv,mdm)

############################################################
#-- Fixed contributions
#-- These are NOT being optimized, only fed into CO2 vector
############################################################

#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_082616/total/oco2_total_fossil.nc")
#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_111716/total/oco2_total_fossil.nc")
con = nc_open("/discover/nobackup/aschuh/run.v11_oco2_080117/fossil/oco2/oco2_total_0_0_0_0_0_0.nc")
fossil = ncvar_get(con,"xco2_model") *1e6
nc_close(con)

#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_082616/total/oco2_total_bg.nc")
#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_111716/total/oco2_total_bg.nc")
con = nc_open("/discover/nobackup/aschuh/run.v11_oco2_080117/bg/oco2/oco2_total_0_0_0_0_0_0.nc")
bg = ncvar_get(con,"xco2_model") *1e6
nc_close(con)

#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_111716/total/oco2_total_fires.nc")
#fires = ncvar_get(con,"xco2_model") *1e6
#nc_close(con)

#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_082616/total/oco2_total_bg_400.nc")
#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_111716/total/oco2_total_bg_400.nc")
con = nc_open("/discover/nobackup/aschuh/run.v11_oco2_080117/bg400/oco2/oco2_total_0_0_0_0_0_0.nc")
bg400 = ncvar_get(con,"xco2_model") *1e6
nc_close(con)

#fossil[fossil > 1e+8] = 400
#oceans[oceans > 1e+8] = 400
#bg[bg > 1e+8] = mean(bg[bg<1e+8])
#bg400[bg400 > 1e+8] = 400

oco2.fixed = list(fossil=fossil,bg=bg,bg400=bg400,fires=NULL)


#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_082616/total/obspack_total_fossil.nc")
#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_111716/total/obspack_total_fossil.nc")
con = nc_open("/discover/nobackup/aschuh/run.v11_oco2_080117/fossil/obspack/obspack_total_0_0_0_0_0.nc")
fossil = ncvar_get(con,"flask") *scalar_noaa
nc_close(con)

#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_082616/total/obspack_total_bg.nc")
#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_111716/total/obspack_total_bg.nc")
con = nc_open("/discover/nobackup/aschuh/run.v11_oco2_080117/bg/obspack/obspack_total_0_0_0_0_0_0.nc")
bg = ncvar_get(con,"flask") *scalar_noaa
nc_close(con)

#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_111716/total/obspack_total_fires.nc")
#fires = ncvar_get(con,"flask") 
#nc_close(con)

#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_082616/total/obspack_total_bg_400.nc")
#con = nc_open("/discover/nobackup/aschuh/GEOS-CHEM_output/oco2_111716/total/obspack_total_bg_400.nc")
con = nc_open("/discover/nobackup/aschuh/run.v11_oco2_080117/bg400/obspack/obspack_total_0_0_0_0_0_0.nc")
bg400 = ncvar_get(con,"flask") * scalar_noaa
nc_close(con)

fossil[fossil > 1e+8] = 400
#oceans[oceans > 1e+8] = 400
#bg[bg > 1e+8 ] = mean(bg[bg!=0])
#bg400[bg400 > 1e+8 ] = 400

obspack.fixed = list(fossil=fossil[good_obspack],bg=bg[good_obspack],bg400=bg400[good_obspack],fires=NULL)

############################################################
############################################################
#-- CREATING JACOBIANS FOR OCEAN, LAND AND CROPS
############################################################
############################################################

###############################################
#-- Build up ocean pulse matrix (Jacobian) for OCO-2
#-- observation operator.
###############################################

#-- This is for OCO2 ------------

#fls.ocn = list.files("/discover/nobackup/aschuh/run.v11_oco2_080117/ocean/oco2/",recursive=FALSE,pattern="oco2_total",full.names=TRUE)

#fls.ocn.short = list.files("/discover/nobackup/aschuh/run.v11_oco2_080117/ocean/oco2/",recursive=FALSE,pattern="oco2_total",full.names=FALSE)

#registerDoSNOW(makeCluster(4, type = "SOCK"))

#datmat.ocean = foreach(i=1:length(fls.ocn),.combine=rbind,.packages=c("ncdf4")) %dopar%
#  {
#    cat(i)
#    idx = as.vector(as.numeric(unlist(strsplit(substring(fls.ocn.short[i],12,26),"_"))))
#	
#	con = nc_open(fls.ocn[i])
#    dat = ncvar_get(con,"xco2_model")
#   nc_close(con)
#  print(paste("idx:",idx," :",length(c(idx,dat))))
#  return(c(idx,dat))
#   } 
       
#save(datmat.ocean,file="/discover/nobackup/aschuh/datmat_ocean.091317.rda")
#load("/discover/nobackup/aschuh/datmat_ocean.080916.rda")
load("/discover/nobackup/aschuh/datmat_ocean.091317.rda")

ind.obs = datmat.ocean > 1e+10 | is.na(datmat.ocean)

###############################################################
#- don't think I need this, it only fills the very 
#- end of the record where we're missing data but not using it
################################################################

datmat.ocean[,-(1:6)] = datmat.ocean[,-(1:6)] * 1e6

print(paste("replacing",sum(ind.obs),"obs with 400",sep=" "))

datmat.ocean[ind.obs] = 400

#-- This essentially flips the ocean signal around the background (reverses sign)

datmat.ocean[,7:(dim(datmat.ocean)[2])] =  t(apply(datmat.ocean[,7:(dim(datmat.ocean)[2])],1,
      FUN=function(x){oco2.fixed$bg400-x}))
      	  
#datmat.ocean[,5:(dim(datmat.ocean)[2])] =  apply(datmat.ocean[,5:(dim(datmat.ocean)[2])],1,FUN=function(x){bg400 - (x-bg400)})

#-- End of OCO2 PULL ------------

##########################################################
#--  Reads output "pulses for ocean" in terms of *NOAA* obspack
#--  observation operator, then constructs Jacobian as above
##########################################################

#fls.ocn = list.files("/discover/nobackup/aschuh/run.v11_oco2_080117/ocean/obspack/",recursive=FALSE,pattern="obspack_total",full.names=TRUE)

#fls.ocn.short = list.files("/discover/nobackup/aschuh/run.v11_oco2_080117/ocean/obspack/",recursive=FALSE,pattern="obspack_total",full.names=FALSE)


#registerDoSNOW(makeCluster(12, type = "SOCK"))

#datmat.ocean.obspack = foreach(i=1:length(fls.ocn),.combine=rbind,.packages=c("ncdf4")) %dopar%
#     {
#    cat(i)
#    idx = as.vector(as.numeric(unlist(strsplit(substring(fls.ocn.short[i],15,29),"_"))))
#	
#	con = nc_open(fls.ocn[i])
#    dat = ncvar_get(con,"flask")[good_obspack] * 1e6 #*scalar_noaa
#   nc_close(con)
#  print(paste("idx:",idx," :",length(c(idx,dat))))
#  return(c(idx,dat))
#   } 
       
#save(datmat.ocean,file="/discover/nobackup/aschuh/datmat_ocean.080916.rda")
#load("/discover/nobackup/aschuh/datmat_ocean.080916.rda")
#save(datmat.ocean.obspack,file="/discover/nobackup/aschuh/datmat_ocean_obspack.091317.rda")
load("/discover/nobackup/aschuh/datmat_ocean_obspack.091317.rda")
#ind.obs = datmat.ocean.obspack > 1e+10 | is.na(datmat.ocean.obspack)

#**********************
#datmat.ocean.obspack[ind.obs] = 400
#**********************

#datmat.ocean.obspack[,7:(dim(datmat.ocean.obspack)[2])] =  t(apply(datmat.ocean.obspack[,7:(dim(datmat.ocean.obspack)[2])],1,
#      FUN=function(x){obspack.fixed$bg400 + obspack.fixed$bg400-x}))

datmat.ocean.obspack[,7:(dim(datmat.ocean.obspack)[2])] =  t(apply(datmat.ocean.obspack[,7:(dim(datmat.ocean.obspack)[2])],1,
      FUN=function(x){(obspack.fixed$bg400-x)}))
       
#-- End of OBSPACK PULL ------------

##########################################################
#--  Reads output "pulses for land" in terms of OCO2 XCO2
#--  observation operator, then saves Jacobian at end
##########################################################

fls = list.files("/discover/nobackup/aschuh/run.v11_oco2_080117/land/oco2/",recursive=FALSE,pattern="oco2_total",full.names=TRUE)

#fls.short =  list.files("/discover/nobackup/aschuh/run.v11_oco2_080117/land/oco2/",recursive=FALSE,pattern="oco2_total",full.names=FALSE)

#######################################
#-- RUN ONCE, THEN LOAD FROM THEN ON
#######################################

con = nc_open(fls[1])
tim = ncvar_get(con,"time")
dat = ncvar_get(con,"xco2_model")
xco2 = ncvar_get(con,"xco2_observed")   	
xco2_sd = ncvar_get(con,"xco2_uncertainty")
nc_close(con)

library(lubridate)
ind_2014_obs_oco2 = tim < time_length(interval(ymd("1985-01-01"), ymd("2015-01-01")),"hour")
ind_2015_obs_oco2 = tim >= time_length(interval(ymd("1985-01-01"), ymd("2015-01-01")),"hour") & tim < time_length(interval(ymd("1985-01-01"), ymd("2016-01-01")),"hour")
ind_2016_obs_oco2 = tim >= time_length(interval(ymd("1985-01-01"), ymd("2016-01-01")),"hour") & tim < time_length(interval(ymd("1985-01-01"), ymd("2017-01-01")),"hour")


#registerDoSNOW(makeCluster(12, type = "SOCK"))

#datmat = 
#foreach(i=1:length(fls),.packages=c("ncdf4"),.combine=rbind) %dopar%
#    {
#   cat(i)
# idx = as.vector(as.numeric(unlist(strsplit(substring(fls.short[i],12,26),"_"))))
#	con = nc_open(fls[i])
#  dat = ncvar_get(con,"xco2_model")
#   nc_close(con)
#  return(c(idx,dat[ind_2014_obs_oco2]))
#  } 
       
#save(datmat,file="/discover/nobackup/aschuh/datmat_2014.091317.rda")

#datmat = 
#foreach(i=1:length(fls),.packages=c("ncdf4"),.combine=rbind) %dopar%
#    {
#   cat(i)
# idx = as.vector(as.numeric(unlist(strsplit(substring(fls.short[i],12,26),"_"))))
#	con = nc_open(fls[i])
#  dat = ncvar_get(con,"xco2_model")
#   nc_close(con)
#  return(c(idx,dat[ind_2015_obs_oco2]))
#  } 
       
#save(datmat,file="/discover/nobackup/aschuh/datmat_2015.091317.rda")

#datmat = 
#foreach(i=1:length(fls),.packages=c("ncdf4"),.combine=rbind) %dopar%
#    {
#   cat(i)
# idx = as.vector(as.numeric(unlist(strsplit(substring(fls.short[i],12,26),"_"))))
#	con = nc_open(fls[i])
#  dat = ncvar_get(con,"xco2_model")
#   nc_close(con)
#  return(c(idx,dat[ind_2016_obs_oco2]))
#  } 
       
#save(datmat,file="/discover/nobackup/aschuh/datmat_2016.091317.rda")


##########################################
#--------------
##########################################

#load("/discover/nobackup/aschuh/datmat.091516.rda")       
#load("/discover/nobackup/aschuh/datmat.120116.rda")
#datmat.bio = datmat

#load("/discover/nobackup/aschuh/datmat_2014.091317.rda")
#datmat_oco2_land_2014 = datmat

load("/discover/nobackup/aschuh/datmat_2015.091317.rda")
datmat_oco2_land_2015 = datmat
rm(datmat)
datmat.bio.oco2 = datmat_oco2_land_2015
rm(datmat_oco2_land_2015)

datmat.bio.oco2[,-(1:6)] = datmat.bio.oco2[,-(1:6)] * 1e6

datmat.bio.oco2[,7:(dim(datmat.bio.oco2)[2])] =  t(apply(datmat.bio.oco2[,7:(dim(datmat.bio.oco2)[2])],1,
      FUN=function(x){(oco2.fixed$bg400[ind_2015_obs_oco2]-x)}))
      
ind.obs = datmat.bio.oco2 > 1e+10 | is.na(datmat.bio.oco2)

#**********************
print(paste("replacing",sum(ind.obs),"obs w/ 399.5"))
datmat.bio.oco2[ind.obs] = 399.5
#**********************

rm(ind.obs)
gc()

##########################################################

#########################
#--  Build up land pulse matrix for OBSPACK
#########################

#fls = list.files("/discover/nobackup/aschuh/run.v11_oco2_080117/land/obspack/",recursive=FALSE,pattern="obspack_total",full.names=TRUE)

#fls.short =  list.files("/discover/nobackup/aschuh/run.v11_oco2_080117/land/obspack/",recursive=FALSE,pattern="obspack_total",full.names=FALSE)

#######################################
#-- RUN ONCE, THEN LOAD FROM THEN ON
#######################################

tim = noaa_time[good_obspack]

library(lubridate)
ind_2014_obs_obspack = tim < time_length(interval(ymd("1970-01-01"), ymd("2015-01-01")),"second")
ind_2015_obs_obspack = tim >= time_length(interval(ymd("1970-01-01"), ymd("2015-01-01")),"second") & tim < time_length(interval(ymd("1970-01-01"), ymd("2016-01-01")),"second")
ind_2016_obs_obspack = tim >= time_length(interval(ymd("1970-01-01"), ymd("2016-01-01")),"second") & tim < time_length(interval(ymd("1970-01-01"), ymd("2017-01-01")),"second")




#registerDoSNOW(makeCluster(28, type = "SOCK"))

#datmat.obspack =
#foreach(i=1:length(fls),.export=c("fls","fls.short"),.packages=c("ncdf4"),.combine=rbind) %dopar%
#     {
#    cat(i)
# idx = as.vector(as.numeric(unlist(strsplit(substring(fls.short[i],15,30),"_"))))
#
#	con = nc_open(fls[i])
#  dat = ncvar_get (con,"flask")[good_obspack][ind_2014_obs_obspack]
#   nc_close(con)
#  return(c(idx,dat))
#} 

#save(datmat.obspack,file="/discover/nobackup/aschuh/datmat.obspack_2014.091317.rda")

#datmat.obspack =
#foreach(i=1:length(fls),.packages=c("ncdf4"),.combine=rbind) %dopar%
#     {
#    cat(i)
# idx = as.vector(as.numeric(unlist(strsplit(substring(fls.short[i],15,30),"_"))))
#
#	con = nc_open(fls[i])
#  dat = ncvar_get (con,"flask")[good_obspack][ind_2015_obs_obspack]
#   nc_close(con)
#  return(c(idx,dat))
#} 

#save(datmat.obspack,file="/discover/nobackup/aschuh/datmat.obspack_2015.091317.rda")

#datmat.obspack =
#foreach(i=1:length(fls),.packages=c("ncdf4"),.combine=rbind) %dopar%
#     {
#    cat(i)
# idx = as.vector(as.numeric(unlist(strsplit(substring(fls.short[i],15,30),"_"))))
#
#	con = nc_open(fls[i])
#  dat = ncvar_get (con,"flask")[good_obspack][ind_2016_obs_obspack]
#   nc_close(con)
#  return(c(idx,dat))
#} 

#save(datmat.obspack,file="/discover/nobackup/aschuh/datmat.obspack_2016.091317.rda")
       

##########################################
#--------------
##########################################

#load("/discover/nobackup/aschuh/datmat.obspack.091516.updated.rda")       
#load("/discover/nobackup/aschuh/datmat.obspack.120116.rda") 

load("/discover/nobackup/aschuh/datmat.obspack_2015.091317.rda")
datmat.bio.obspack = datmat.obspack
rm(datmat.obspack)
gc()

datmat.bio.obspack[,7:(dim(datmat.bio.obspack)[2])] = datmat.bio.obspack[,7:(dim(datmat.bio.obspack)[2])] * scalar_noaa

datmat.bio.obspack[,7:(dim(datmat.bio.obspack)[2])] =  t(apply(datmat.bio.obspack[,7:(dim(datmat.bio.obspack)[2])],1,
      FUN=function(x){(obspack.fixed$bg400[ind_2015_obs_obspack]-x)}))
      
      
ind.obs = datmat.bio.obspack > 1e+10 | is.na(datmat.bio.obspack) 

#**********************
print(paste("replacing",sum(ind.obs),"obs w/ 399.5"))
datmat.bio.obspack[ind.obs] = 399.5
#**********************

rm(ind.obs)
gc()

############################################################
############################################################
#-- END OF CREATING JACOBIANS FOR OCEAN, LAND AND CROPS
############################################################
############################################################

###########################################################################
#-- Builds up "bias" pulse by tying in covariates from native OCO-2 files
###########################################################################

#-- biastype=1 is simple dataset dependent IC offset
#-- biastype=2 is more elaborate "full" bias correction

biastype = 2

#-- Give OBSPACK TYPE=300 identifier
type = c(type_orig,rep(300,length(obspack.fixed$bg)))

land_ind = type %in% c(110,111,112,113)

land_glint = type == 111
land_nadir = type == 110
land_trans = type == 113
land_target = type == 112

ocean_ind = type %in% c(101)

tccon_ind = type %in% c(200)

land_bias_logdws = rep(0,length(type))
land_bias_logdws[land_ind] = logdwsv[land_ind] 
land_bias_logdws[land_ind] = land_bias_logdws[land_ind]  -mean(land_bias_logdws[land_ind])

land_bias_co2_grad_del = rep(0,length(type))
land_bias_co2_grad_del[land_ind] = co2_grad_delv[land_ind] 
land_bias_co2_grad_del[land_ind] = land_bias_co2_grad_del[land_ind] -mean(land_bias_co2_grad_del[land_ind])

land_bias_s31 = rep(0,length(type))
land_bias_s31[land_ind] = s31v[land_ind] 
land_bias_s31[land_ind] = land_bias_s31[land_ind] # -mean(land_bias_s31[land_ind])

#land_bias_sza = rep(0,length(type))
#land_bias_sza[land_glint] = szav[land_glint] 

#land_bias_sza = land_bias_sza - 40
#land_bias_sza[land_bias_sza <=0 ] = 0

land_bias_dp = rep(0,length(type))
land_bias_dp[land_ind] = dpv[land_ind] 
print(paste("removing",mean(land_bias_dp[land_ind]),"from","land dp"))
land_bias_dp[land_ind] = land_bias_dp[land_ind]-mean(land_bias_dp[land_ind])

bias_airmass = rep(0,length(type))
#bias_airmass[!land_nadir & !tccon_ind] = airmassv[!land_nadir & !tccon_ind]
#bias_airmass[!land_nadir & !tccon_ind] = bias_airmass[!land_nadir & !tccon_ind] -mean(bias_airmass[!land_nadir & !tccon_ind])
biased_airmass_ind = land_glint | land_trans | land_target | ocean_ind
bias_airmass[biased_airmass_ind] = airmassv[biased_airmass_ind]
bias_airmass[biased_airmass_ind] = bias_airmass[biased_airmass_ind] -mean(bias_airmass[biased_airmass_ind])


land_bias_constant_nadir = rep(0,length(type))
#land_bias_constant[land_nadir] = (1/0.9955)
#land_bias_constant[land_glint] = (1/0.997)
#land_bias_constant[land_trans] = (1/0.997)
#land_bias_constant[land_target] = (1/0.99625)
land_bias_constant_nadir[land_nadir] = 1

land_bias_constant_glint = rep(0,length(type))
land_bias_constant_glint[land_glint | land_target] = 1

ocean_bias_co2_grad_del = rep(0,length(type))
ocean_bias_co2_grad_del[ocean_ind] = co2_grad_delv[ocean_ind]
print(paste("removing",mean(ocean_bias_co2_grad_del[ocean_ind]),"from","ocean co2graddel"))
ocean_bias_co2_grad_del[ocean_ind] = ocean_bias_co2_grad_del[ocean_ind] - mean(ocean_bias_co2_grad_del[ocean_ind])


ocean_bias_constant = rep(0,length(type))
ocean_bias_constant[ocean_ind] = 1

ocean_bias_dp = rep(0,length(type))
ocean_bias_dp[ocean_ind] = dpv[ocean_ind]
print(paste("removing",mean(ocean_bias_dp[ocean_ind]),"from","ocean dp"))

ocean_bias_dp[ocean_ind] = ocean_bias_dp[ocean_ind] - mean(ocean_bias_dp[ocean_ind])

overall_IC_bias = rep(1,length(type))

if(biastype==1)
{
 datmat.bias =  matrix(overall_IC_bias,nrow=1)

datmat.bias = cbind(matrix(rep(1000,5),nrow=1),datmat.bias)


}
if(biastype==2)
{
 datmat.bias =  rbind(land_bias_constant_nadir,land_bias_constant_glint,land_bias_logdws,land_bias_co2_grad_del,land_bias_dp,
                                land_bias_s31,
                                #land_bias_sza,
                                ocean_bias_co2_grad_del,ocean_bias_dp,ocean_bias_constant) #,xco2_rawv)


datmat.bias = cbind(cbind(rep(1000,dim(datmat.bias)[1]),rep(1000,dim(datmat.bias)[1]),rep(1000,dim(datmat.bias)[1]),rep(1000,dim(datmat.bias)[1]),rep(1000,dim(datmat.bias)[1])),datmat.bias)
}

#fixed_bias = rep(0,length(type))
#fixed_bias[land_glint] =  xco2_rawv[land_glint]*(0.003/0.997)
#fixed_bias[land_nadir] =  xco2_rawv[land_nadir]*(0.0045/0.9955)
#fixed_bias[land_trans] =  xco2_rawv[land_trans]*(0.00375/0.99625)
#fixed_bias[land_target] =  xco2_rawv[land_target]*(0.003/0.997)
#fixed_bias[ocean_ind] =  xco2_rawv[ocean_ind]*(0.003/0.997)

#gpp = apply(datmat[-(1:4),datmat[2,]==1 & datmat[3,]==0]-400,1,sum)
#resp = apply(datmat[-(1:4),datmat[2,]==0 & datmat[3,]==0]-400,1,sum)

#ind =  datmat[,2]==0 & !(datmat[,3] %in% c(11,12,13) )

#-- This worked for just priors
#ind =  datmat[,2]==0 
#datmat = datmat[ind,]
###---------------------------


#datmat[datmat[,1]==1,-(1:4)] = -datmat[datmat[,1]==1,-(1:4)]

#resp = apply(datmat[datmat[,1]==0,-(1:4)],2,sum,na.rm=TRUE)
#gpp = apply(datmat[datmat[,1]==1,-(1:4)],2,sum,na.rm=TRUE)

resp_xco2 = apply(datmat.bio.oco2[datmat.bio.oco2[,2]==0 & datmat.bio.oco2[,4]==0 & datmat.bio.oco2[,5] != 25,-(1:6)],2,sum,na.rm=TRUE)
resp_noaa = apply(datmat.bio.obspack[datmat.bio.obspack[,2]==0 & datmat.bio.obspack[,4]==0 & datmat.bio.obspack[,5] != 25,-(1:6)],2,sum,na.rm=TRUE)

crops_xco2 = apply(datmat.bio.oco2[datmat.bio.oco2[,2]==0 & datmat.bio.oco2[,4]==0 & datmat.bio.oco2[,5] == 25,-(1:6)],2,sum,na.rm=TRUE)
crops_noaa = apply(datmat.bio.obspack[datmat.bio.obspack[,2]==0 & datmat.bio.obspack[,4]==0 & datmat.bio.obspack[,5] == 25,-(1:6)],2,sum,na.rm=TRUE)

gpp_xco2 = apply(datmat.bio.oco2[datmat.bio.oco2[,2]==1 & datmat.bio.oco2[,4]==0 & datmat.bio.oco2[,5] != 25,-(1:6)],2,sum,na.rm=TRUE)
gpp_noaa = apply(datmat.bio.obspack[datmat.bio.obspack[,2]==1 & datmat.bio.obspack[,4]==0 & datmat.bio.obspack[,5] != 25,-(1:6)],2,sum,na.rm=TRUE)


#ocean.flux = apply(apply(datmat.ocean[datmat.ocean[,2]==0,-(1:4)],1,FUN=function(x){x-bg400}),1,sum,na.rm=TRUE) 
#ocean.flux_xco2 =  apply(apply(datmat.ocean[datmat.ocean[,4]==0,-(1:6)],1,FUN=function(x){x-oco2.fixed$bg400}),1,sum,na.rm=TRUE)
ocean.flux_xco2 =  cbind(datmat.ocean[datmat.ocean[,4]==0,1:6],datmat.ocean[datmat.ocean[,4]==0,ind_2015_obs_oco2])
ocean.flux_xco2 =  apply(ocean.flux_xco2[,-(1:6)],2,sum,na.rm=TRUE)

#ocean.flux_noaa =  apply(apply(datmat.ocean.obspack[datmat.ocean.obspack[ind_2014_obs_oco2,3]==0,-(1:6)],1,FUN=function(x){x-obspack.fixed$bg400[ind_2014_obs_oco2]}),1,sum,na.rm=TRUE)
ocean.flux_noaa =  cbind(datmat.ocean.obspack[datmat.ocean.obspack[,4]==0,1:6],datmat.ocean.obspack[datmat.ocean.obspack[,4]==0,ind_2015_obs_obspack])
ocean.flux_noaa =  apply(ocean.flux_noaa[,-(1:6)],2,sum,na.rm=TRUE)

#resp[resp>1e+8] = 399
#resp[resp>1e+8] = 399

nee = c(resp_xco2,resp_noaa) - c(gpp_xco2,gpp_noaa)

#nee.leftover = apply(datmat[datmat[,2]!=0 & ind,-(1:4)],2,sum,na.rm=TRUE)

#-- THIS IS FINAL MODEL A PRIORI PREDICTION OF XCO2 FROM MODEL ******

#xco2_model_total = nee+ c(crops,crops.obspack) + ( c(oco2.fixed$fossil,obspack.fixed$fossil)-c(oco2.fixed$bg400,obspack.fixed$bg400) )+ c(ocean.flux_xco2,ocean.flux_noaa) + c(oco2.fixed$bg,obspack.fixed$bg)  #+ fixed_bias

xco2_model_total = nee  + 
#c(crops,crops.obspack) + 
( c(oco2.fixed$fossil,obspack.fixed$fossil)-c(oco2.fixed$bg400,obspack.fixed$bg400) ) +
#      ( c(oco2.fixed$fires,obspack.fixed$fires)-c(oco2.fixed$bg400,obspack.fixed$bg400) ) +
         c(ocean.flux_xco2,ocean.flux_noaa) + c(oco2.fixed$bg,obspack.fixed$bg)

#xco2_fixed_flux =  (crops) + (fossil-bg400) + bg


#######################################################################
#---- GOOOD DATA IS ONLY 1:175203 (175204:188125 from Dec 2015 is screwed up)
#######################################################################


#-- THIS IS DIAGNOSTIC PLOT OF MODEL XCO2 VS OBSERVED XCO2

library(zoo)

#indind = 1:length(xco2_model_total) %in% c(1:134466,(142638+1):(142638+67298))

indind = 1:length(xco2_model_total) %in%  c(1:155312,(161090+1):(161090+67125))

#xco2_model_total = xco2_model_total - adjj  - 0.5

aa = filter(xco2_model_total[indind],rep(1/500,500))

#plot((nee+(fossil-399.5)+(oceans-399.5))[seq(1,274000,by=100)],pch=16,cex=0.25,col="red")

#bb = rollapply(xco2_rawv+ fixed_bias,500,FUN=mean,na.rm=TRUE)

bb = rollapply(c(xco2v,realflask[good_obspack])[indind],500,FUN=mean,na.rm=TRUE)

#points(xco2[seq(1,274000,by=100)],pch=16,cex=0.25,col="blue")

#adjj = (1:length(crops))/length(crops) * (0.23*11)

rng = range(c(aa,bb),na.rm=TRUE)

plot(aa,col="red",type="l",ylim=c(390,418))

points(bb,col="blue",type="l")



#sdd = apply(datmat[datmat[,2]==0,-(1:5)],1,sd,na.rm=TRUE)

#indx = datmat[datmat[,2]==0,1:5][order(sdd),]

#--
#-- These are all the pulses that are exactly equal after 150000 samples, probably those = bg400 pulse
#ind_none = datmat[,132600] == datmat[1,132600]

bias=FALSE

if(bias){
	datmat = rbind(cbind(datmat.bio,datmat.bio.obspack[,-c(1:5)]),cbind(datmat.ocean,datmat.ocean.obspack[,-c(1:5)]),datmat.bias)
	}else{
       datmat = rbind(cbind(datmat.bio,datmat.bio.obspack[,-c(1:5)]),cbind(datmat.ocean,datmat.ocean.obspack[,-c(1:5)]))
      }
      

##########################################################################################
#-  This is the INVERSION CODE.  Basic solution for a Bayesian regression with Gaussian 
#-  prior for coefficients.  This can be set up as a function, or stepped through for
#-  single applications.
##########################################################################################

invert = function(data_types,error_floor_perc,land_prior_sd,ocean_prior_sd)
{ 

#land_prior_sd = 0.3;ocean_prior_sd = 0.1;error_floor_perc = 0;data_types = c(101,110,111,300,200)

#ocean_prior_sd = 0.1

#error_floor_perc = 0

#data_types = c(101,110,111,300,200)

sigma = matrix(0,ncol=dim(datmat)[1],nrow=dim(datmat)[1])

for(i in 1:(dim(sigma)[1]))
{
	if(datmat[i,2]==0) {
	    if(i < 3697){
	   	   sigma[i,i] =   land_prior_sd^2    #0.001^2
	    }
	    if(i > 3696 & i < 3846){
	     	sigma[i,i] =   ocean_prior_sd^2    #0.001^2
	   	  }
	   	if(i >= 3847){
	     # if( i == 3790 | i == 3793 )
	   	  #{
	   	 # 	  sigma[i,i] =   0.00001^2    #0.001^2
         #}else{
	      sigma[i,i] =   land_prior_sd^2    #0.001^2
	    	}
	   	 }else{
	   	 	
	    if(i < 3697){
	   	   sigma[i,i] =   land_prior_sd^2    #0.001^2
	    }
	    if(i > 3696 & i < 3847){
	     	sigma[i,i] =   ocean_prior_sd^2    #0.001^2
	   	  }
	   	if(i >= 3847){
	   	
	   	  #if( i == 3790 | i == 3793 )
	   	 #   {
	   	 #    sigma[i,i] =   0.00001^2    #0.001^2
         #  }else{
	     	sigma[i,i] =   0.3^2    #0.001^2
	   	 }
	 }
}

#ind = rep(TRUE,length(xco2_sd))
#ind = type != 101 & type != 200
#ind =  type != 101 & 1:(dim(datmat)[2]-4) %in% 1:175203 & abs(xco2_model_total) < 1000
#-- 101 is ocean glint

#data_types = c(300)

########################################
#-- Subsetting to about March 1, 2016
########################################

#ind = (type %in% data_types)  & 1:(dim(datmat)[2]-5) %in% c(1:112974,120061:166790) & abs(xco2_model_total) < 1000

#ind = (type %in% data_types)  & 1:(dim(datmat)[2]-5) %in% c(1:134466,(142638+1):(142638+67298)) & abs(xco2_model_total) < 1000 & abs(xco2_model_total) > 0

ind = (type %in% data_types)  & 1:(dim(datmat)[2]-5) %in% c(1:155312,(161090+1):(161090+67125)) & abs(xco2_model_total) < 1000 & abs(xco2_model_total) > 0


#ind = (type %in% data_types)  & 1:(dim(datmat)[2]-4) %in% 1:110000 & abs(xco2_model_total) < 1000

X = datmat[,-(1:5)] 

#X_rest = X[,!ind]
#X = X[,ind]

if(bias){
	
	 bias_len = dim(datmat)[1]-(3696 + 150)
	  #XXX=matrix(rep(c(oco2.fixed$bg400,obspack.fixed$bg400)[ind],dim(datmat)[1]-bias_len),nrow=dim(datmat)[1]-bias_len,byrow=TRUE)
	  #X[1:(dim(datmat)[1]-bias_len),] = X[1:(dim(datmat)[1]-bias_len),]-XXX
	  
	  #XXX=matrix(rep(c(oco2.fixed$bg400,obspack.fixed$bg400)[!ind],dim(datmat)[1]-bias_len),nrow=dim(datmat)[1]-bias_len,byrow=TRUE)
	  #X_rest[1:(dim(datmat)[1]-bias_len),] = X_rest[1:(dim(datmat)[1]-bias_len),]-XXX
	  	
	  #XXX=matrix(rep(c(oco2.fixed$bg400,obspack.fixed$bg400),dim(datmat)[1]-bias_len),nrow=dim(datmat)[1]-bias_len,byrow=TRUE)
	  
	  #X[1:(dim(datmat)[1]-bias_len),] = X[1:(dim(datmat)[1]-bias_len),]-XXX
	  	  	  
	  X[1:(dim(datmat)[1]-bias_len),] = sweep(X[1:(dim(datmat)[1]-bias_len),],2,c(oco2.fixed$bg400,obspack.fixed$bg400),'-')	 
	   	  
	  }else{
      #XXX=matrix(rep(c(oco2.fixed$bg400,obspack.fixed$bg400)[ind],dim(datmat)[1]),nrow=dim(datmat)[1],byrow=TRUE)
      #X = X-XXX
      
      #XXX=matrix(rep(c(oco2.fixed$bg400,obspack.fixed$bg400)[!ind],dim(datmat)[1]),nrow=dim(datmat)[1],byrow=TRUE)
      #X_rest = X_rest-XXX    
      
      #XXX=matrix(rep(c(oco2.fixed$bg400,obspack.fixed$bg400),dim(datmat)[1]),nrow=dim(datmat)[1],byrow=TRUE)
      #X = X-XXX  
      
      X = sweep(X,2,c(oco2.fixed$bg400,obspack.fixed$bg400),'-')	
 }
 
#X_full = X
X_rest = X[,!ind]
X = X[,ind]
 
#XXX=matrix(rep(bg400[ind],dim(datmat)[1]-1),nrow=dim(datmat)[1]-1,byrow=TRUE)

#-- Adjustment is to include bias terms

#X[1:(dim(datmat)[1]-1),] = X[1:(dim(datmat)[1]-1),]-XXX
#X[1:(dim(datmat)[1]-7),] = X[1:(dim(datmat)[1]-7),]-XXX


#-- Here I'm pirating xco2_sdv by adding mdm on to end
co2_sdv_full = c(xco2_sdv,mdm[good_obspack])
co2_sdv = co2_sdv_full[ind]

#-- WATCH HERE, ADDING error FLOOR*******
#xco2_sdv[type==101] = xco2_sdv[type==101]  + 0.5*error_floor_perc
#xco2_sdv[type %in% c(110,111,112)] = xco2_sdv[type %in% c(110,111,112)]  + 0.7*error_floor_perc
#xco2_sd_final = xco2_sdv[ind] 

#R_vector = xco2_sd_final^-1 
#R2_vector = xco2_sd_final^-2

R_vector = co2_sdv^-1 
R2_vector = co2_sdv^-2

R = matrix(rep(R_vector,(dim(X)[1])),byrow=TRUE,ncol=length(R_vector))
R2 = matrix(rep(R2_vector,(dim(X)[1])),byrow=TRUE,ncol=length(R2_vector))

XR = X*R
XR2 = X*R2

#rm(R)
#rm(R2)

tX_X = tcrossprod(XR)

#rm(XR)

part1 = solve(sigma) + tX_X

#rm(tX_X)

#-- This is the posterior covariance matrix of the state

part2 = solve(part1)

#K = part2 %*% XR2

#-- RUNNING W/ RAW XCO2 NOW !!!!!
if(bias &  biastype!=1){
	  #part3 = matrix((xco2_rawv-xco2_model_total)[ind],ncol=1)
	  part3 = matrix((c(xco2_rawv,realflask[good_obspack])-xco2_model_total)[ind],ncol=1)
}else{
       #part3 = matrix((xco2v-xco2_model_total)[ind],ncol=1)
       part3 = matrix((c(xco2v,realflask[good_obspack])-xco2_model_total)[ind],ncol=1)
       
       #-- FOR S31
       #part3 = matrix((c(xco2v + land_bias_s31[1:length(xco2v)],realflask[good_obspack])-xco2_model_total)[ind],ncol=1)
       
       #-- FOR CLIM ERRORS FROM TRANSPORT
       load("/discover/nobackup/aschuh/clim.err.rda")
       
       mons = cumsum(c(0,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31)
       resid = (timev - min(timev))/(3600*24)+6
       m_ind = sapply(resid,FUN=function(x){max(grep(TRUE,x > mons))}) %% 12
       m_ind[m_ind==0] = 12
       
       #-- 
       lats = seq(-90,90,length=91)
       lat_ind = sapply(latv,FUN=function(x){max(grep(TRUE,x > lats))}) 
       
       resid_transport = clim.err[cbind(lat_ind,m_ind)]
       part3 = matrix((c(resid_transport,realflask[good_obspack]))[ind],ncol=1)
       
}
#part3 = matrix((xco2_rawv-xco2_model_total)[ind],ncol=1)
#part3 = matrix(( (xco2_rawv + fixed_bias)-xco2_model_total)[ind],ncol=1)

part4 = XR2  %*% part3

#rm(XR2)

part5 = part2 %*% part4

#rm(part4)
#rm(part3)

adj = t(part5) %*% X 

adj = as.vector(adj)

adj_rest = t(part5) %*% X_rest

adj_rest = as.vector(adj_rest)

full_adj = rep(0,length(xco2_model_total))

full_adj[ind] = adj

full_adj[!ind] = adj_rest

#XXX=matrix(rep(bg400[ind],528),nrow=528,byrow=TRUE)

#adj2 = t(part5) %*% (X)*2.135

#adj2 = as.vector(adj2)

#-- RUNNING W/ RAW XCO2 NOW !!!!!

if(bias & biastype!=1){
	resids.prior.bias = (c(xco2_rawv,realflask[good_obspack]) - xco2_model_total)[ind]
	resids.post.bias = (c(xco2_rawv,realflask[good_obspack])- xco2_model_total)[ind] - adj
	}else{
      resids.prior = (c(xco2v,realflask[good_obspack]) - xco2_model_total)[ind]
      resids.post = (c(xco2v,realflask[good_obspack])- xco2_model_total)[ind] - adj
     }
     
xco2.post = xco2_model_total[ind] + adj #  - 2*part5[length(part5)]

full.xco2.post = xco2_model_total + full_adj

outdat = cbind(time=c(timev,noaa_time[good_obspack]),lon=c(lonv,noaa_lon[good_obspack]),lat=c(latv,noaa_lat[good_obspack]),type=type,xco2.prior =xco2_model_total ,xco2.post = full.xco2.post, assim_TF=ind,obs=c(xco2v,realflask[good_obspack]) )

 #test = t(sweep(X_full[,type==200],1,as.vector(part5),'*')) %*% X_full
 #test2 = sweep(test,MARGIN=2,co2_sdv_full^-2,'*')
 #rm(test)
 #rm(X_full)
 
return(list(sigma=sigma,part2=part2,part5=part5,outdat=outdat))
}

###################################
#-- END OF INVERSION CODE
###################################

#--- S, the influence matrix
#--   S = diag(R2_vector[1:10000])%*%t(X[,1:10000])%*%part2%*%(X[,1:10000])


#-- THIS IS IF WE WANT TO RUN OVER A VARIETY OF INVERSION CASES

#data.sets = list(c(101),c(110),c(111,112),c(101,110),c(110,111,112),c(101,111,112),c(101,110,111,112),c(200))

data.sets = list(c(101,110),c(101,110,300,200),c(101),c(110),c(300),c(200,300),c(200),c(111),c(101,300),c(110,300),c(111,300))

error.perc = c(0,0.5,1,2)

land.pr.sdv = c(0.1,0.2,0.3)

ocean.pr.sdv = c(0.05,0.1,0.2)

gridrun = expand.grid(1:11,1,3,1)

#gridrun = rbind(gridrun[56:66,],gridrun[-(56:66),])

#a = foreach(i=1:(dim(gridrun)[1])) %dopar% {

#error.perc = c(0,0.5,1,2)

land.pr.sdv = seq(0.05,0.3,by=0.025)

ocean.pr.sdv = seq(0.025,0.325,by=0.03)

gridrun = expand.grid(5,1,1:11,1:11)

for(i in 1:(dim(gridrun)[1]))
{
	print(i)
	x = gridrun[i,]
	print(paste("x:",x))
	
	out = invert(data.sets[[as.numeric(x[1])]],
	              error.perc[as.numeric(x[2])],
	              land.pr.sdv[as.numeric(x[3])],
	              ocean.pr.sdv[as.numeric(x[4])])
	
	save(out,file=paste("/discover/nobackup/aschuh/oco2_inversion_suite_S31/oco2_inversion_",paste(x,collapse="_"),".rda",sep=""))
	
}





plot.dat = cbind(resids.prior=resids.prior,resids.post=resids.post,
                          time=timev[ind],lat=latv[ind])
                          
#save(plot.dat,file="/discover/nobackup/aschuh/plot.dat.052516.rda")                          

calc.objs = c("tX_X","part1","part2","part3","part4","part5","XR","XR2","X","R","R2","XXX")

#-- Plotting 


#hist(resids.prior,breaks=100,col="red",xlim=c(-10,10),freq=FALSE,ylim=c(0,0.5),main="Residuals (Model-OCO2)",xlab="CO2 (ppm)")
#points(density(resids.post),col="blue",lwd=4,type="l")
#legend(-10,0.4,c("Prior residuals","Posterior residuals"),col=c("red","blue"),lty=c(1,1),lwd=c(4,4))

#-- Comparison of resids w/bias vs w/o bias
d1 = density(resids.prior)
d2 = density(resids.post)
d3 = density(resids.prior.bias)
d4 = density(resids.post.bias)

plot(d1$x,d1$y,col="red",lty=4,xlim=c(-5,5),ylim=c(0,0.6),type="l",xlab="resids",ylab="density")
points(d2$x,d2$y,col="blue",lty=4,xlim=c(-5,5),type="l")
points(d3$x,d3$y,col="red",lty=1,xlim=c(-5,5),type="l")
points(d4$x,d4$y,col="blue",lty=1,xlim=c(-5,5),type="l")
grid()
legend(1.5,0.5,c(paste("BC Prior mn=",round(mean(resids.prior),2)," sd=",round(sd(resids.prior),2),sep=""),
paste("BC Post mn=",round(mean(resids.post),2)," sd=",round(sd(resids.post),2),sep=""),
paste("RAW Prior mn=",round(mean(resids.prior.bias),2)," sd=",round(sd(resids.prior.bias),2),sep=""),
paste("RAW Post mn=",round(mean(resids.post.bias),3)," sd=",round(sd(resids.post.bias),2),sep="")),col=c("red","blue","red","blue"),lty=c(4,4,1,1),cex=0.7)



#####################################################################
#-- Plot of model residuals, oriented towards bias terms right now
#####################################################################

aa = filter(xco2_model_total[ind],rep(1/500,500))

aa2 = filter(xco2_model_total[ind] + as.vector(adj -adjB),rep(1/500,500))

#plot((nee+(fossil-399.5)+(oceans-399.5))[seq(1,274000,by=100)],pch=16,cex=0.25,col="red")

bb = filter(xco2_rawv[ind] - as.vector(adjB),rep(1/500,500))

#bb = filter((xco2_rawv + fixed_bias)[ind],rep(1/500,500))

plot(aa,col="pink",type="l",ylim=c(394,404))

points(aa2,col="red",type="l",ylim=c(394,401))

points(bb,col="blue",type="l")

points(filter(xco2.post,rep(1/500,500)),col="green",type="l")


#--- 

adjB = as.vector(t(part5[3847:(dim(part5)[1]),1]) %*% X[3847:(dim(part5)[1]),])
adjB_O = as.vector(t(part5[3852:3854,1]) %*% X[3852:3854,])
adjB_L = as.vector(t(part5[3847:3851,1]) %*% X[3847:3851,])
#adjB_airmass = as.vector(t(part5[3855,1]) %*% X[3855,])
adjO = as.vector(t(part5[3697:3846,1]) %*% X[3697:3846,])
adjL = as.vector(t(part5[1:3696,1]) %*% X[1:3696,])

#-- Plotting bias correction

colv = rev(rainbow(6))[match(type,sort(unique(type)))]

#plot((xco2v - xco2_rawv)[ind][seq(1,161090,by=100)],-adjB[seq(1,161090,by=100)],ylim=c(-1,4),xlim=c(-1,4),bg=colv[ind][seq(1,176000,by=100)],cex=0.75,pch=21,xlab="OCO-2 (BC - RAW)",ylab="Inversion estimated BC")

#-- With S31 posthoc, already pulled out 0.2 from s31v

#ttype =101

#plot((xco2v + 7*(s31v)- xco2_rawv)[ind][type_orig[ind] == ttype][seq(1,161090,by=40)],-adjB[type_orig[ind] == ttype][seq(1,161090,by=40)],ylim=c(-1,4),xlim=c(-1,4),bg=colv[ind][type_orig[ind] == ttype][seq(1,161090,by=40)],cex=0.75,pch=21,xlab="OCO-2 (BC - RAW)",ylab="Inversion estimated BC")

plot((xco2v + 7*(s31v)- xco2_rawv)[ind][seq(1,161090,by=40)],-adjB[seq(1,161090,by=40)],ylim=c(-1,4),xlim=c(-1,4),bg=colv[ind][seq(1,161090,by=40)],cex=0.75,pch=21,xlab="OCO-2 (BC - RAW)",ylab="Inversion estimated BC")

legend(3,1,c("OG","LN","LG","LTarg","LTran","TCCON"),pt.bg=rev(rainbow(6)),pch=21,bg="white")

lines(c(-5,5),c(-5,5),lwd=3,col="grey")

grid()

############
#############
################
require(maps)
w = map("world",plot=FALSE)

grd = expand.grid(lon = seq(-177.5,179,by=5),lat=seq(-89,89,by=4) )


spat.dat.full = data.frame(TYPE=type[ind],LAT=c(latv,noaa_lat[good_obspack])[ind],LON=c(lonv,noaa_lon[good_obspack])[ind],MYBIAS=-adjB,OCO2BIAS=(xco2v + 7*(s31v)- xco2_rawv)[ind])

spat.dat.full$DIFF = spat.dat.full$MYBIAS - spat.dat.full$OCO2BIAS

#spat.dat.full$DIFF = spat.dat.full$DIFF - mean(spat.dat.full$DIFF,na.rm=TRUE)

spat.dat = spat.dat.full[!is.na(spat.dat.full$DIFF) & spat.dat.full$TYPE == 110,]


glbl.n <- t(matrix(.Fortran("agrid2d",
                  PACKAGE = "biocycle", n = as.integer(length(spat.dat$LON)),
                  x = as.double(spat.dat$LON),
                  y = as.double(spat.dat$LAT),
                   val = as.double(rep(1,length(spat.dat$DIFF))), szx = as.integer(72),
                  szy = as.integer(45), grdx = as.double(seq(-180,
                    180, by = 5)), grdy = as.double(seq(-90,
                    90, by = 4)),
                  array = matrix(as.double(0), ncol = 45,
                    nrow = 72))$array, ncol = 72,
                  byrow = T))
                  
glbl <- t(matrix(.Fortran("agrid2d",
                  PACKAGE = "biocycle", n = as.integer(length(spat.dat$LON)),
                  x = as.double(spat.dat$LON),
                  y = as.double(spat.dat$LAT),
                   val = as.double(spat.dat$DIFF), szx = as.integer(72),
                  szy = as.integer(45), grdx = as.double(seq(-180,
                    180, by = 5)), grdy = as.double(seq(-90,
                    90, by = 4)),
                  array = matrix(as.double(0), ncol = 45,
                    nrow = 72))$array, ncol = 72,
                  byrow = T))

res = glbl/glbl.n

grd$LN = as.vector(res)

spat.dat = spat.dat.full[!is.na(spat.dat.full$DIFF) & spat.dat.full$TYPE == 111,]


glbl.n <- t(matrix(.Fortran("agrid2d",
                  PACKAGE = "biocycle", n = as.integer(length(spat.dat$LON)),
                  x = as.double(spat.dat$LON),
                  y = as.double(spat.dat$LAT),
                   val = as.double(rep(1,length(spat.dat$DIFF))), szx = as.integer(72),
                  szy = as.integer(45), grdx = as.double(seq(-180,
                    180, by = 5)), grdy = as.double(seq(-90,
                    90, by = 4)),
                  array = matrix(as.double(0), ncol = 45,
                    nrow = 72))$array, ncol = 72,
                  byrow = T))
                  
glbl <- t(matrix(.Fortran("agrid2d",
                  PACKAGE = "biocycle", n = as.integer(length(spat.dat$LON)),
                  x = as.double(spat.dat$LON),
                  y = as.double(spat.dat$LAT),
                   val = as.double(spat.dat$DIFF), szx = as.integer(72),
                  szy = as.integer(45), grdx = as.double(seq(-180,
                    180, by = 5)), grdy = as.double(seq(-90,
                    90, by = 4)),
                  array = matrix(as.double(0), ncol = 45,
                    nrow = 72))$array, ncol = 72,
                  byrow = T))

res = glbl/glbl.n

grd$LG = as.vector(res)

spat.dat = spat.dat.full[!is.na(spat.dat.full$DIFF) & spat.dat.full$TYPE == 101,]


glbl.n <- t(matrix(.Fortran("agrid2d",
                  PACKAGE = "biocycle", n = as.integer(length(spat.dat$LON)),
                  x = as.double(spat.dat$LON),
                  y = as.double(spat.dat$LAT),
                   val = as.double(rep(1,length(spat.dat$DIFF))), szx = as.integer(72),
                  szy = as.integer(45), grdx = as.double(seq(-180,
                    180, by = 5)), grdy = as.double(seq(-90,
                    90, by = 4)),
                  array = matrix(as.double(0), ncol = 45,
                    nrow = 72))$array, ncol = 72,
                  byrow = T))
                  
glbl <- t(matrix(.Fortran("agrid2d",
                  PACKAGE = "biocycle", n = as.integer(length(spat.dat$LON)),
                  x = as.double(spat.dat$LON),
                  y = as.double(spat.dat$LAT),
                   val = as.double(spat.dat$DIFF), szx = as.integer(72),
                  szy = as.integer(45), grdx = as.double(seq(-180,
                    180, by = 5)), grdy = as.double(seq(-90,
                    90, by = 4)),
                  array = matrix(as.double(0), ncol = 45,
                    nrow = 72))$array, ncol = 72,
                  byrow = T))

res = glbl/glbl.n

grd$OG = as.vector(res)



levelplot(OG ~ lon + lat,data=grd,col.regions=my.col(20),
     at=seq(-max(abs(grd$OG),na.rm=TRUE),max(abs(grd$OG),na.rm=TRUE),length=20),
      panel = function(x,y,z,...)
                      {
                        panel.contourplot(x,y,z,...)
                        llines(w$x,w$y,col="black",lty=4)
                      })

levelplot(LN ~ lon + lat,data=grd,col.regions=my.col(20),
     at=seq(-max(abs(c(grd$LN,grd$LG)),na.rm=TRUE),max(abs(c(grd$LN,grd$LG)),na.rm=TRUE),length=20),
      panel = function(x,y,z,...)
                      {
                        panel.contourplot(x,y,z,...)
                        llines(w$x,w$y,col="black",lty=4)
                      })
                      
levelplot(LG ~ lon + lat,data=grd,col.regions=my.col(20),
     at=seq(-max(abs(c(grd$LN,grd$LG)),na.rm=TRUE),max(abs(c(grd$LN,grd$LG)),na.rm=TRUE),length=20),
      panel = function(x,y,z,...)
                      {
                        panel.contourplot(x,y,z,...)
                        llines(w$x,w$y,col="black",lty=4)
                      })






plot(rollapply(adjB[seq(1,176000,by=100)],500,FUN=mean,na.rm=TRUE),pch=0.4,col="black",ylim=c(-2,2))
points(rollapply(adjB_L[seq(1,176000,by=100)],500,FUN=mean,na.rm=TRUE),type="l",col="green",lty=4,ylim=c(-2,2))
points(rollapply(adjB_O[seq(1,176000,by=100)],500,FUN=mean,na.rm=TRUE),type="l",col="blue",lty=4,ylim=c(-2,2))
points(rollapply(adjB_airmass[seq(1,176000,by=100)],500,FUN=mean,na.rm=TRUE),type="l",col="red",lty=4,ylim=c(-2,2))
points(rollapply(adjO[seq(1,176000,by=100)],500,FUN=mean,na.rm=TRUE),pch=0.4,col="blue",ylim=c(-2,2))
points(rollapply(adjL[seq(1,176000,by=100)],500,FUN=mean,na.rm=TRUE),pch=0.4,col="green",ylim=c(-2,2))
points(rollapply((adjB+adjO+adjL)[seq(1,176000,by=100)],500,FUN=mean,na.rm=TRUE),pch=0.4,col="grey",ylim=c(-2,2))
#points(rollapply(((xco2_rawv + fixed_bias)-xco2_model_total)[ind],500,FUN=mean,na.rm=TRUE),pch=0.4,col="orange")


plot(rollapply(((xco2v)-xco2_model_total)[ind][seq(1,50000,by=50)],500,FUN=mean,na.rm=TRUE),pch=0.4,col="orange",ylim=c(-3,3))
points(rollapply((adjO+adjL)[seq(1,50000,by=50)],500,FUN=mean,na.rm=TRUE),pch=0.4,col="grey",ylim=c(-2,2))
points(rollapply(adjO[seq(1,50000,by=50)],500,FUN=mean,na.rm=TRUE),pch=0.4,col="blue",ylim=c(-2,2))
points(rollapply(adjL[seq(1,50000,by=50)],500,FUN=mean,na.rm=TRUE),pch=0.4,col="green",ylim=c(-2,2))
points(rollapply(adjB[seq(1,50000,by=50)],500,FUN=mean,na.rm=TRUE),pch=0.4,col="magenta",ylim=c(-2,2))
points(rollapply((adjB+adjL+adjO)[seq(1,50000,by=50)],500,FUN=mean,na.rm=TRUE),lwd=2,type="l",col="black",lty=4,ylim=c(-2,2))

#-- Sort by reduction in variance for scalar

cbind(datmat[,1:5],diag(part2)/diag(sigma))[rev(order(diag(part2)/diag(sigma))),]


######################################
######################################
#--   CALCULATE POSTERIOR FLUXES
######################################
######################################

coefs = rbind(rep(1,365)
,sin((1:365)/365*2*pi*1)
,sin((1:365)/365*2*pi*2)
,sin((1:365)/365*2*pi*3)
,cos((1:365)/365*2*pi*1)
,cos((1:365)/365*2*pi*2)
,cos((1:365)/365*2*pi*3)
,rep(1,365)
,sin((1:365)/365*2*pi*1)
,sin((1:365)/365*2*pi*2)
,sin((1:365)/365*2*pi*3)
,cos((1:365)/365*2*pi*1)
,cos((1:365)/365*2*pi*2)
,cos((1:365)/365*2*pi*3))



#coefs.ocean = rbind(rep(1,365)
#,sin((1:365)/365*2*pi*1)
#,cos((1:365)/365*2*pi*1))

######

fil = nc_open("/discover/nobackup/aschuh/data/sib_harmonic/sib4_1deg_0014.lu.nc.orig")

ref = ncvar_get(fil,"pft_ref")

area = ncvar_get(fil,"pft_area")

assim = ncvar_get(fil,"assim")

resp = ncvar_get(fil,"resp_balanced")

nc_close(fil)

fil = nc_open("/discover/nobackup/aschuh/data/sib_harmonic/sib4_1deg_0014.regional.g.nc")

resp_crop = ncvar_get(fil,"resp_crop")

nc_close(fil)

######


fil = nc_open("/discover/nobackup/aschuh/data/sib_harmonic/regions.nc")

regions = ncvar_get(fil,"transcom_regions")

grid_cell_area = ncvar_get(fil,"grid_cell_area")
o_regions = ncvar_get(fil,"ocean_regions")
xform = ncvar_get(fil,"xform")

#-- Not sure why we're using 210:239 instead of 211:240
o_xform = xform[12:22,210:239]

nc_close(fil)

#registerDoMC(12)

registerDoSNOW(makeCluster(12, type = "SOCK"))

#for(i in 1:11)
acomb <- function(...) abind(..., along=5)

#






#-- YOU HAVE TO RUN THIS 

retr_harm = function(time)
{
	c(1,
	   sin(time/365*2*pi*1),
	   sin(time/365*2*pi*2),
	   sin(time/365*2*pi*3),
	   cos(time/365*2*pi*1),
	   cos(time/365*2*pi*2),
	   cos(time/365*2*pi*3),
	   1,
	   sin(time/365*2*pi*1),
	   sin(time/365*2*pi*2),
	   sin(time/365*2*pi*3),
	   cos(time/365*2*pi*1),
	   cos(time/365*2*pi*2),
	   cos(time/365*2*pi*3))	
}


fil = nc_open("/discover/nobackup/aschuh/data/sib_harmonic/sib4_flux_20140906.nc")

ref = ncvar_get(fil,"pft_ref")

area = ncvar_get(fil,"area")

nc_close(fil)


fil = nc_open("/discover/nobackup/aschuh/data/sib_harmonic/regions.nc")

regions = ncvar_get(fil,"transcom_regions")

grid_cell_area = ncvar_get(fil,"grid_cell_area")
o_regions = ncvar_get(fil,"ocean_regions")
xform = ncvar_get(fil,"xform")

#-- Not sure why we're using 210:239 instead of 211:240
o_xform = xform[12:22,210:239]

nc_close(fil)

#-- Containers for fluxes

resp_tot = array(0,dim=c(360,180,10,2))
assim_tot = array(0,dim=c(360,180,10,2))

mon_days = cumsum(c(1,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31))

yrs = c(rep(2014,4),rep(2015,12),rep(2016,8))

mths = c(9:12,1:12,1:8)



gen_1x1 = function(fil)
{
load(fil)
load("/discover/nobackup/aschuh/oco2_inversion_postproc_092916/datmat_index.092716.rda")

#out =  cbind(datmat[,1:5],part5,diag(sigma)^0.5,diag(part2)^0.5)
out =  cbind(datmat_index[,1:5],out$part5,diag(out$sigma)^0.5,diag(out$part2)^0.5)
out[3697:3846,5] = out[3697:3846,5] + 100

registerDoSNOW(makeCluster(12, type = "SOCK"))

reg3D = abind(regions,regions,regions,regions,regions,regions,regions,regions,regions,regions,along=3)

solution_land = array(0,dim=c(360,180,10,7,2))

for(i in 1:11)
{
	print(paste("i:",i))
	for(j in 1:24)
	{
		cat(j)
		part5_sol = out[out[,4] == j & out[,5] == i,]
		for(k in 1:7)
		{
			solution_land[,,,k,1][ref == j & reg3D ==i] = part5_sol[k,6]
			solution_land[,,,k,2][ref == j & reg3D ==i] = part5_sol[7+k,6]			
		}
		
	}
}		
		
solution_land2 = aaply(solution_land,c(4,5),.fun=function(x){x*area})
solution_land2 = aperm(solution_land2,c(3,4,5,1,2))
		
#-- TIME LOOP

require(snow)

for(k in 1:24)
{
	fls.tmp = sort(list.files("/discover/nobackup/aschuh/data/sib_harmonic/",pattern=paste("flux_",yrs[k],pad(mths[k],2,"0"),sep=""),full.names=TRUE))
	
	numb.files = length(fls.tmp)


acomb5 <- function(...) abind(..., along=5)


aaa =  foreach(m=1:numb.files,.combine=acomb5,.packages=c("ncdf4","abind"))  %dopar%
 {
  				conn = nc_open(fls.tmp[m])
  				resp = ncvar_get(conn,"resp_bal")
  				assim = ncvar_get(conn,"assim")  	
  				resp_crop = ncvar_get(conn,"resp_crop")
  				nc_close(conn)
  				
  				resp_new = array(0,dim=dim(resp))
  				assim_new = array(0,dim=dim(assim))
        
                resp_final = array(0,dim=c(dim(resp)[1:3]))
                assim_final = array(0,dim=c(dim(assim)[1:3]))
                                
                #-- Run starts on 9/6, pulling fluxes from 9/1                
                harm_tmp = sapply(243 + (mon_days[k]-1)+(m-1)+((1:24)/24),retr_harm)
  				
  				resp_prior = apply(resp,c(1,2,3),mean) * area 
  				assim_prior = apply(assim,c(1,2,3),mean) * area
  				
  				#-- solution * harmonic
  				for(w in 1:24)
  				{
  				cat(w)
  				mult.resp = sweep(solution_land2[,,,,1],c(4),harm_tmp[1:7,w],"*")
  				mult.assim = sweep(solution_land2[,,,,2],c(4),harm_tmp[8:14,w],"*")  			
  					
  				 mult.resp = apply(mult.resp,c(1,2,3),sum)
  				 mult.assim = apply(mult.assim,c(1,2,3),sum)
  				   					
  				#mult.resp= aaply(solution_land2[,,,,1],c(4),.fun=function(x){x*harm_tmp[1:7,]})
  				#mult.assim= aaply(solution_land2[,,,,2],c(4),.fun=function(x){x*retr_harm()})  				
  	
  				resp_new[,,,w] = resp[,,,w] *mult.resp
  				assim_new[,,,w] = assim[,,,w] *mult.assim
  				}
  				
  				resp_final = apply(resp_new,c(1,2,3),mean)
   				assim_final = apply(assim_new,c(1,2,3),mean)
   				 				
  				#assim = aaply(mult.assim,c(1,2,3),.fun=function(x){x*assim})
  				  				
  				#resp_tot[,,,1] =  resp_tot[,,,1]  + 1/numb.files*resp_prior
  				#assim_tot[,,,1] =  assim_tot[,,,1]  + 1/numb.files*assim_prior
  				
  				#resp_tot [,,,2] = resp_tot [,,,2] + 1/numb.files*resp_final
  				#assim_tot[,,,2] = assim_tot[,,,2] + 1/numb.files*assim_final	
  				
  				#abind(resp_tot,assim_tot,along=4)  		
  				abind(1/numb.files*resp_prior,1/numb.files*resp_final,
  				         1/numb.files*assim_prior,1/numb.files*assim_final	,along=4	)
  			}
  			
  			
  			
  			
  			
  			return(abind(resp_tot ,assim_tot,along=4))
  			
}


rr = apply(resp_prior,c(1,2),sum) + apply(resp_crop,c(1,2),mean)

aa = apply(assim_prior,c(1,2),sum)

rr2 = apply(resp_final,c(1,2),sum) 
aa2 = apply(assim_final,c(1,2),sum)

plot(apply(rr-aa,2,mean),col="black",pch=16,ylim=c(-2,2))

#-- look for place where the "+" comes from before gpp adjust: aa2

points(apply((rr-aa) + (rr2+aa2),2,mean),col="blue",pch=16)







#return_diag = function(out,resp,assim,regions,area,ref,coefs)
#{
outtt = foreach(i=1:11,.combine='acomb',.export=c("pad","out","mons","resp","assim","regions","area","ref","coefs","abind"),.noexport=c("datmat",calc.objs),.multicombine=TRUE,.verbose=TRUE) %dopar%
  {
  	print(objects())
  	#print(dim(datmat))
  	resp_tot = array(0,dim=c(360,180,12))
    assim_tot = array(0,dim=c(360,180,12))
   
    resp_tot2 = array(0,dim=c(360,180,12))
   assim_tot2 = array(0,dim=c(360,180,12)) 
    
    PTH = array(NA,dim=c(2*7*24*11,5))
    
  	for( j in 1:24)
  	{
  		print(paste("i:",i," j:",j))
  		comp = out[out[,4] == j & out[,5] == i,]
  		ord = order(comp[,1],comp[,2],comp[,3])
  		comp = comp[ord,]
  		
  		#coef2 = apply(coefs,2,FUN=function(x){x*comp[,5]})
  		#coef3 = apply(coefs,2,FUN=function(x){x*comp[,6]})
  		#coef4 = apply(coefs,2,FUN=function(x){x*comp[,7]})
  		 
  		coef2 = apply(coefs,2,FUN=function(x){x*comp[,6]})
  		coef3 = apply(coefs,2,FUN=function(x){x*comp[,7]})
  		coef4 = apply(coefs,2,FUN=function(x){x*comp[,8]})
  		  		  		
  		ref_ind = ref == j
  		
  		resp_tmp = array(0,dim=dim(resp_tot))
  		
  		assim_tmp = array(0,dim=dim(resp_tot))
  		
  		resp_tmp2 = array(0,dim=dim(resp_tot))
  		
  		assim_tmp2 = array(0,dim=dim(resp_tot))

  		#resp_tmp3 = array(0,dim=(dim(resp)[c(1,2,4)]))
  		
  		#assim_tmp3 = array(0,dim=(dim(resp)[c(1,2,4)]))
  		  		  		
  		#resp_tmp4 = array(0,dim=(dim(resp)[c(1,2,4)]))
  		
  		#assim_tmp4 = array(0,dim=(dim(resp)[c(1,2,4)]))
  		
  		for(k in 1:24)
  		{ 
  			area_tmp = area
  			area_tmp[!ref_ind] = 0
  			
  			fls.tmp = sort(list.files("/discover/nobackup/aschuh/data/sib_harmonic/",pattern=paste("flux_",yrs[k],pad(mths[k],2,"0"),sep=""),full.names=TRUE))
  			
  			for(m in 1:length(fls.tmp))
  			{
  				conn = nc_open(fls.tmp[m])
  				resp = ncvar_get(conn,"resp_bal")
  				assim = ncvar_get(conn,"assim")  				
  			}
  			
  			resp_tmp[,,k] = apply(resp[,,,k]* area_tmp,c(1,2),sum)
  			assim_tmp[,,k] = apply(assim[,,,k]* area_tmp,c(1,2),sum)
  			
  			#resp_tmp[,,,k] = resp_tmp[,,,k] * area
  			#assim_tmp[,,,k] = assim_tmp[,,,k] * area
  			
  			resp_tmp[,,k][regions!=i] = 0
  			assim_tmp[,,k][regions!=i] = 0

  			resp_tmp2[,,k] = resp_tmp[,,k] 
  			assim_tmp2[,,k] = assim_tmp[,,k] 
  			  			
  			resp_tmp[,,k] = resp_tmp[,,k] * mean(coef2[1:7,mons[k]:(mons[k+1]-1)])*7
  			assim_tmp[,,k] = assim_tmp[,,k] * mean(coef2[8:14,mons[k]:(mons[k+1]-1)])*7  		
  				
  		}
  		
  		#-- Have flux for PFT *j* and TRANSOM REGION *i*
  		
  		assim_tot = assim_tot + assim_tmp[,,]
  		resp_tot = resp_tot + resp_tmp[,,]  		
  		
    	assim_tot2 = assim_tot2 + assim_tmp2[,,]
  		resp_tot2 = resp_tot2 + resp_tmp2[,,]  				
  		
  		#assim_tot3 = assim_tot3 + assim_tmp3[,,]
  		#resp_tot3 = resp_tot3 + resp_tmp3[,,]  		
  		
  		#assim_tot4 = assim_tot4 + assim_tmp4[,,]^2
  		#resp_tot4 = resp_tot4 + resp_tmp4[,,]^2 
  	}
  abind(assim_tot,assim_tot2,resp_tot,resp_tot2,along=4) #,assim_tot3,resp_tot3,assim_tot4,resp_tot4,along=4)
}

out_ocean = foreach(i=101:130,.combine='acomb',.packages=c("ncdf4"),.export=c("o_regions","pad","out","resp","regions","area","ref","coefs","abind","mons"),.noexport=c("datmat",calc.objs),.multicombine=TRUE,.verbose=TRUE) %dopar%
  {
  	print(objects())
  	#print(dim(datmat))

    ocean_tot = array(0,dim=(dim(resp)[c(1,2,4)]))
    ocean_tot2 = array(0,dim=(dim(resp)[c(1,2,4)]))    
    
  		print(paste("i:",i))
  		comp = out[out[,5] == i,]
  		ord = order(comp[,2])
  		comp = comp[ord,]
  		
  		#coef2 = apply(coefs,2,FUN=function(x){x*comp[,5]})
  		#dim(coef2)
  		#ref_ind = ref == j
  		
  		#ocean_tmp = array(0,dim=(dim(resp)[c(1,2,4)]))
  		ocean_tmp = array(0,c(360,180,22))
  		#ocean_tmp2 = array(0,dim=(dim(resp)[c(1,2,4)]))
  		ocean_tmp2 = array(0,c(360,180,22))
  		
  	    mon.ocean = c(9,10,11,12,1:12,1:12)
  	    
  	    yr.ocean = c(rep(2014,4),rep(2015,12),rep(2016,12))
  	    
  	    
  		for(k in 1:22)
  		{
  			cat(k)
  			#fil.tmp = nc_open(paste("/discover/nobackup/aschuh/data/CT2015_priors/CT2015.prior_",yr.ocean[k],pad(mon.ocean[k],2,fill="0"),".nc",sep=""))
            fil.tmp = nc_open(paste("/discover/nobackup/aschuh/data/CT_NRT_priors_monthly/CT-NRT.v2016-1.unopt.flux1x1.",yr.ocean[k],pad(mon.ocean[k],2,fill="0"),".nc",sep=""))
            #-- NOT OI FLUXES, using Andy's NRT fluxes
  			oi = ncvar_get(fil.tmp,"ocn_flux_opt")
  			nc_close(fil.tmp)  		
           mm=cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30,31))
           tsteps = dim(oi)[3]
           mons2 = c(mons,mons[13]+(mons[-1]-1))
           strt = mons2[mon.ocean[k]]-1

            ocean_indic = o_regions
            ocean_indic[,] = 0 
            ocean_indic[o_regions == (i-100)] = 1
            	
  			#-- Need to cut out non ocean region fluxes at this point ********
  			  			
  			#area_tmp = area
  			#area_tmp[!ref_ind] = 0
  			
  			#-- Apply harmonic
  			
  			harms.tmp = array(NA,dim=c(3,tsteps))
            #harms.tmp[1,]  = rep(1,tsteps) * comp[1,5]
            #harms.tmp[2,]  = sin(strt + (1:tsteps)*1/8*2*pi*1) * comp[2,5]
            #harms.tmp[3,]  = cos(strt + (1:tsteps)*1/8*2*pi*1)       * comp[3,5]    
            harms.tmp[1,]  = rep(1,tsteps) * comp[1,6]
            harms.tmp[2,]  = sin(strt + (1:tsteps)*1/8*2*pi*1) * comp[2,6]
            harms.tmp[3,]  = cos(strt + (1:tsteps)*1/8*2*pi*1)       * comp[3,6]
            
  			#--  Apply grid cell areas
  			
  			ocean_tmp[,,k] = ocean_indic*apply(oi,c(1,2),FUN=function(x){mean(x*harms.tmp[1,] + x*harms.tmp[2,] + x*harms.tmp[3,] )})
  			ocean_tmp2[,,k] = ocean_indic*apply(oi,c(1,2),FUN=function(x){mean(x)}) 			  	  	   
  				
  		}
  		
  		#assim_tot = assim_tot + assim_tmp[,,]
  		#resp_tot = resp_tot + resp_tmp[,,]  		
  		
    	#assim_tot2 = assim_tot2 + assim_tmp2[,,]
  		#resp_tot2 = resp_tot2 + resp_tmp2[,,]  				
  	    #}
  	abind(ocean_tmp,ocean_tmp2,along=4)
  #abind(assim_tot,assim_tot2,resp_tot,resp_tot2,along=4)
}

out_list = list(outtt=outtt,out_ocean=out_ocean)
}



#outtt = return_diag(out=out,resp=resp,assim=assim,regions=regions,area=area,ref=ref,coefs=coefs)


registerDoSEQ()

#levelplot(nee_annual_adj*10^-6*12*3600*24*365,col.regions=rainbow(100),cuts=100)

###############################
#--  LAND
###############################


nee_mats = apply(outtt,c(1,2,3,4),sum)

nee_prior = nee_mats[,,,4] - nee_mats[,,,2] + resp_crop

#nee_post = nee_mats[,,,3] + nee_mats[,,,1] + nee_prior + resp_crop

#-- This is effective land sink
nee_annual_adj = apply(nee_mats[,,,3] + nee_mats[,,,1],c(1,2),mean)

nee_resp_adj = apply(nee_mats[,,,3] ,c(1,2),mean)

nee_gpp_adj = apply(nee_mats[,,,1],c(1,2),mean)

gpp_annual_prior = apply(nee_mats[,,,2],c(1,2),mean)

gpp_annual_post = -apply(nee_mats[,,,1],c(1,2),mean) + gpp_annual_prior

resp_annual_prior = apply(nee_mats[,,,4],c(1,2),mean) + apply(resp_crop,c(1,2),mean)

resp_annual_post = apply(nee_mats[,,,3],c(1,2),mean) + resp_annual_prior

nee_annual_prior= resp_annual_prior  - gpp_annual_prior

nee_annual_post = resp_annual_post  - gpp_annual_post

print(paste("PRIOR NEE:",format(sum(grid2transcom(nee_annual_prior*10^-6*12*3600*24*365)[24:46]),2)))

print(paste("NEE ADJ:",format(sum(grid2transcom(nee_annual_adj*10^-6*12*3600*24*365)[24:46]),2)))

print(paste("POST NEE:",format(sum(grid2transcom(nee_annual_post*10^-6*12*3600*24*365)[24:46]),2)))

print(paste("PRIOR ER:",format(sum(grid2transcom(resp_annual_prior*10^-6*12*3600*24*365)[24:46]),2)))
print(paste("PRIOR GPP:",format(sum(grid2transcom(gpp_annual_prior*10^-6*12*3600*24*365)[24:46]),2)))

print(paste("POST ER:",format(sum(grid2transcom(resp_annual_post*10^-6*12*3600*24*365)[24:46]),2)))
print(paste("POST GPP:",format(sum(grid2transcom(gpp_annual_post*10^-6*12*3600*24*365)[24:46]),2)))

print(paste("NH:",format(sum(grid2transcom(nee_annual_post*10^-6*12*3600*24*365)[24:46][c(1,2,7,8,11)]),2)))
print(paste("SH:",format(sum(grid2transcom(nee_annual_post*10^-6*12*3600*24*365)[24:46][c(3,4,5,6,9,10)]),2)))

###############################
#--  OCEANS
###############################


ocean_mats = apply(out_ocean,c(1,2,3,4),sum)

ocean_prior = ocean_mats[,,,2]

#-- This is effective ocean sink adjustment

ocean_annual_adj = apply(ocean_mats[,,1:12,2] + ocean_mats[,,1:12,1],c(1,2),mean)


ocean_annual_prior = apply(ocean_mats[,,1:12,2],c(1,2),mean)  #  - need ocean sink of -3 PgC

ocean_annual_post = apply(ocean_mats[,,1:12,1],c(1,2),mean) + apply(ocean_mats[,,1:12,2],c(1,2),mean)  #  - need ocean sink of -3 PgC

print(paste("PRIOR NEE:",format(sum(grid2transcom(ocean_annual_prior*12*3600*24*365)[35:46]),2)))

print(paste("OCEAN POST:",format(sum(grid2transcom(ocean_annual_post*12*3600*24*365)[35:46]),2)))

print(paste("NH:",format(sum(grid2transcom(nee_annual_post*10^-6*12*3600*24*365)[24:46][c(1,2,7,8,11)]),2)))
print(paste("SH:",format(sum(grid2transcom(nee_annual_post*10^-6*12*3600*24*365)[24:46][c(3,4,5,6,9,10)]),2)))


######################################
######################################
#--   END OF CALCULATE POSTERIOR FLUXES
######################################
######################################



#-- Plotting

#-- bar plots

return_t3 = function(outtt)
{
nee_mats = apply(outtt,c(1,2,3,4),sum)
nee_prior = nee_mats[,,,4] - nee_mats[,,,2] + resp_crop
nee_annual_adj = apply(nee_mats[,,,3] + nee_mats[,,,1],c(1,2),mean)
nee_resp_adj = apply(nee_mats[,,,3] ,c(1,2),mean)
nee_gpp_adj = apply(nee_mats[,,,1],c(1,2),mean)
gpp_annual_prior = apply(nee_mats[,,,2],c(1,2),mean)
gpp_annual_post = -apply(nee_mats[,,,1],c(1,2),mean) + gpp_annual_prior
resp_annual_prior = apply(nee_mats[,,,4],c(1,2),mean) + apply(resp_crop,c(1,2),mean)
resp_annual_post = apply(nee_mats[,,,3],c(1,2),mean) + resp_annual_prior
nee_annual_prior= resp_annual_prior  - gpp_annual_prior
nee_annual_post = resp_annual_post  - gpp_annual_post
t3 = grid2transcom(nee_annual_post*10^-6*12*3600*24*365)[24:46]
t3 = c(t3[1:11],sum(t3[12:23]))
return(t3)
}

load("/discover/nobackup/aschuh/outtt.all_data.rda")
t3 = return_t3(outtt)
t3_tot = t3

load("/discover/nobackup/aschuh/outtt.no_land_glint.rda")
t3=return_t3(outtt)
t3_tot = rbind(t3_tot,t3)

load("/discover/nobackup/aschuh/outtt.no_ocean_glint.rda")
t3=return_t3(outtt)
t3_tot = rbind(t3_tot,t3)

load("/discover/nobackup/aschuh/outtt.no_glint.rda")
t3=return_t3(outtt)
t3_tot = rbind(t3_tot,t3)

load("/discover/nobackup/aschuh/outtt.no_ocean_glint_0.35.rda")
t3=return_t3(outtt)
t3_tot = rbind(t3_tot,t3)

load("/discover/nobackup/aschuh/outtt.no_ocean_glint_0.7.rda")
t3=return_t3(outtt)
t3_tot = rbind(t3_tot,t3)

t3_tot = cbind(t3_tot,apply(t3_tot[,c(1,2,7,8,11)],c(1),sum),apply(t3_tot[,c(3,4,5,6,9,10)],c(1),sum),apply(t3_tot,c(1),sum))

t3_tot2 = as.matrix(t3_tot)

dimnames(t3_tot2) = list(c("All Data","No Land GL","No Ocean GL","No GL","No Ocean GL + 0.35","No Ocean GL + 0.7")
,c("Boreal NA","NA Temp","SA Trop","SA Temp","N Africa","S Africa","Boreal AS","Temp AS","Trop AS","Aust","Europe","Extra","NH","SH","Global"))



barplot(t3_tot2/1e+15,beside=TRUE, cex.names=1.2,cex.axis = 1.5,cex.lab=1.4,ylab="PgC Flux",las=2,legend.text=dimnames(t3_tot2)[[1]],col=c("black","blue","green","red","cyan","orange"))

grid()

barplot(t3_tot2/1e+15,beside=TRUE, cex.axis = 1.5,ylab="PgC Flux",cex.lab = 1.4,cex.names=1.2,las=2,legend.text=dimnames(t3_tot2)[[1]],col=c("black","blue","green","red","cyan","orange"),add=TRUE)


library(RColorBrewer)

  #-- Load world map
  require(maps)
  mm = map("world",plot=FALSE)

  constant = 10^-6*12*3600*24*365

  #-- Build grid for plotting
  ll.grid = expand.grid(lon=seq(-179.5,179.5,by=1),lat=seq(-89.5,89.5,1))

  ll.grid$nee.prior = as.vector(nee_annual_prior*constant )
  ll.grid$nee.post = as.vector(nee_annual_post*constant )
  ll.grid$nee.post [ll.grid$nee.post > 500] = 500
  ll.grid$nee.post [ll.grid$nee.post < -500] = -500
    
  ll.grid$crops = as.vector(apply(nee_mats[,,,4] - nee_mats[,,,2],c(1,2),mean)*constant )
  
  ll.grid$crops.resp = as.vector(apply(resp_crop,c(1,2),mean)*constant )
  
  ll.grid$crops.resp[ll.grid$crops.resp > 1000] = 1000
  
  ll.grid$nee.adj =  ll.grid$nee.post  -  ll.grid$nee.prior
  ll.grid$nee.adj [ll.grid$nee.adj > 500] = 500
  ll.grid$nee.adj [ll.grid$nee.adj < -500] = -500
  
  
 p1 = levelplot(crops ~ lon + lat,data=ll.grid,
              #col.regions=c(rev(brewer.pal(9,"Blues")),brewer.pal(9,"Greys")),
              #col.regions=rev(brewer.pal(11,"PuOr")),
              col.regions=rev(brewer.pal(9,"Blues")),
              ylab="",xlab="",
              # main="2010 Annual NEE",
              #at=seq(-340,340,length=12),
              cuts=8,
              scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
              par.settings=list(layout.heights=list(left.padding=-3,right.padding=-3)),
              panel = function(x,y,z,...)
                      {
                        panel.contourplot(x,y,z,...)
                        llines(mm$x,mm$y,col="black",lty=4)
                      })
                      
p2 = levelplot(crops.resp ~ lon + lat,data=ll.grid,
              #col.regions=c(rev(brewer.pal(9,"Blues")),brewer.pal(9,"Greys")),
              col.regions=(brewer.pal(9,"Reds")),
                            ylab="",xlab="",
              # main="2010 Annual NEE",
              #at=seq(-340,340,length=12),
              cuts=8,
              scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
              par.settings=list(layout.heights=list(left.padding=-3,right.padding=-3)),
              panel = function(x,y,z,...)
                      {
                        panel.contourplot(x,y,z,...)
                        llines(mm$x,mm$y,col="black",lty=4)
                      })             
                      

p3 = levelplot(nee.post - nee.prior ~ lon + lat,data=ll.grid,
              #col.regions=c(rev(brewer.pal(9,"Blues")),brewer.pal(9,"Greys")),
              #col.regions=rev(brewer.pal(9,"PuOr")),
               #            ylab="",xlab="",
              # main="2010 Annual NEE",
              col.regions=my.col(50),
              cuts=49,
              scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
              par.settings=list(layout.heights=list(left.padding=-3,right.padding=-3)),
              panel = function(x,y,z,...)
                      {
                        panel.contourplot(x,y,z,...)
                        llines(mm$x,mm$y,col="black",lty=4)
                      })  
                      
p4 = levelplot(nee.post ~ lon + lat,data=ll.grid,
              #col.regions=c(rev(brewer.pal(9,"Blues")),brewer.pal(9,"Greys")),
              col.regions=rev(brewer.pal(9,"PuOr")),
                           ylab="",xlab="",
              # main="2010 Annual NEE",
              #at=seq(-340,340,length=12),
              cuts=8,
              scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
              par.settings=list(layout.heights=list(left.padding=-3,right.padding=-3)),
              panel = function(x,y,z,...)
                      {
                        panel.contourplot(x,y,z,...)
                        llines(mm$x,mm$y,col="black",lty=4)
                      })             
                      
                      
print(p1, split=c(1,1,2,2), more=TRUE)
print(p2, split=c(2,1,2,2), more=TRUE)
print(p3, split=c(1,2,2,2), more=TRUE)
print(p4, split=c(2,2,2,2), more=FALSE)


#-- Plotting seasonal cycles

flx = apply(outtt[,,,,4],c(3,4),FUN=function(x){grid2transcom(x*10^-6*12*3600*24*365)[27]})/1e+15

flx_crop = apply(resp_crop,c(3),FUN=function(x){grid2transcom(x*10^-6*12*3600*24*365)[27]})/1e+15

flx_df = cbind(flx,flx_crop)

flx_df = flx_df[c(9:12,1:8),]

flx_df = flx_df / 12

#plot(-flx[,2]/12,type="b",col="blue",lty=4,ylim=c(-1.5,1.5),ylab="GPP <--  PgC  -->  RESP",xaxt="n",xlab="")
plot(-flx_df[,2],type="b",col="blue",lty=4,ylim=c(-1.5,1.5),ylab="GPP <--  PgC  -->  RESP",xaxt="n",xlab="")

axis(1,side=1,at=1:12,labels=c("Sept","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug"))

#points(-(flx[,2]-flx[,1])/12,type="b",col="blue",lty=1)
points(-(flx_df[,2]-flx_df[,1]),type="b",col="blue",lty=1)

#points((flx_crop+flx[,4])/12,type="b",col="red",lty=4)
points((flx_df[,5]+flx_df[,4]),type="b",col="red",lty=4)

#points((flx_crop+flx[,4]+flx[,3])/12,type="b",col="red",lty=1)
points((flx_df[,5]+flx_df[,4]+flx_df[,3]),type="b",col="red",lty=1)

#points((flx_crop+flx[,4]+flx[,3] - (flx[,2]-flx[,1]) )/12,type="b",col="magenta",lty=1,lwd=3)
points((flx_df[,5]+flx_df[,4]+flx_df[,3] - (flx_df[,2]-flx_df[,1]) ),type="b",col="magenta",lty=1,lwd=3)

grid()

legend(10,-0.55,c("Post Resp","Prior Resp","Post GPP","Prior GPP","Post NEE"),
             col=c("red","red","blue","blue","magenta"),lty=c(1,4,1,4,1),lwd=c(1,1,1,1,3),cex=0.85,bg="cyan")

#-- Create flux matrix, FLUXTYPE X HARMONIC X PFT X TRANREGION X TIME




acomb4 <- function(...) abind(..., along=4)

dim(datmat)[1]

sub = grep(   TRUE,  datmat[,1] == 0  & datmat[,4] == REG   )

flux_ts = foreach(i=sub,.verbose=TRUE,.combine='acomb4') %dopar%
  {
   cat(i)
   FTYPE = datmat[i,1]
   HARM   = datmat[i,1]
   PFT      = datmat[i,3]
   REG     = datmat[i,4]   

    ref_ind = ref == PFT

   area_tmp = area
   area_tmp[!ref_ind] = 0

   flux_tmp = array(0,dim=(dim(resp)[c(1,2,4)]))
    
   if(FTYPE==0)
   {
   	 for(k in 1:12){
   	 	        flux_tmp[,,k] = apply(resp[,,,k]* area_tmp,c(1,2),sum)
   	 	     }
    }
    
   if(FTYPE==1)
   {
   	 for(k in 1:12){
   	 	        flux_tmp[,,k] = apply(resp[,,,k]* area_tmp,c(1,2),sum)
   	 	     }
    }
    
    flux_tmp[regions!=REG] = 0
    
    comp = out[out[,3] == PFT & out[,4] == REG & out[,2] == HARM & out[,1] == FTYPE,]
    #ord = order(comp[,1],comp[,2],comp[,3])
  	#comp = comp[ord,]
    
    #coef2 = apply(coefs[FTYPE*7+(HARM+1)],2,FUN=function(x){x*comp[,5]})
    coef2 = coefs[FTYPE*7+(HARM+1),] * comp[5]
    
  for(k in 1:12)
  		{
  			flux_tmp[,,k]  = flux_tmp[,,k] * mean(coef2[mons[k]:(mons[k+1]-1)])
  				
  		}

   flux_tmp
   
   #abind(assim_tot,assim_tot2,resp_tot,resp_tot2,along=4)
}









#-- BIAS STUFF

 dat = cbind(xco2 = xco2v,xco2r = xco2_rawv,dp = dpv + 1.4,airmass=airmassv,co2_grad_del=co2_grad_delv + 8.4,logdws=logdwsv -2.9,type=type)
 dat = as.data.frame(dat)
 lm(xco2-0.995*xco2r ~ dp + airmass + co2_grad_del+ logdws,data=dat[dat$type==101,])
 
 
 
 
 
 
 
#--  This is for output for Experiments

fls_out = c("gridded_flux_template_schuh_SE.nc4",
"gridded_flux_template_schuh_SEi.nc4",
"gridded_flux_template_schuh_OG.nc4",
"gridded_flux_template_schuh_LN.nc4",
"gridded_flux_template_schuh_IS.nc4",
"gridded_flux_template_schuh_TCi.nc4",
"gridded_flux_template_schuh_TC.nc4",
"gridded_flux_template_schuh_LG.nc4",
"gridded_flux_template_schuh_OGi.nc4",
"gridded_flux_template_schuh_LNi.nc4",
"gridded_flux_template_schuh_LGi.nc4")

fls_out = sapply(fls_out,FUN=function(x){paste("/discover/nobackup/aschuh/for_sean_010417/",x,sep="")})

exp_out = sapply(1:11,FUN=function(x){paste("/discover/nobackup/aschuh/oco2_inversion_suite_120716/oco2_inversion_",x,"_1_3_2.rda",sep="")})


for(w in 1:11)
{
	#load(exp_out[w])
	fil_out = fls_out[w]
	dat_out = gen_1x1(exp_out[w])
	con = nc_open(fil_out,write=TRUE)


nee_mats = apply(dat_out$outtt,c(1,2,3,4),sum)

#-- Dump land prior into file

nee_prior = nee_mats[,,,4] - nee_mats[,,,2] + resp_crop

#ncvar_put(con,varid="land_flux",vals=nee_prior,start=c(1,1,1),count=c(-1,-1,12))

#ncvar_put(con,varid="land_flux",vals=nee_prior,start=c(1,1,13),count=c(-1,-1,12))

#-- Dump land posterior into file

nee_post = nee_mats[,,,3] + nee_mats[,,,1] + nee_prior 

ncvar_put(con,varid="land_flux",vals=nee_post*10^-6*10^-3*12,start=c(1,1,1),count=c(-1,-1,12))

ncvar_put(con,varid="land_flux",vals=nee_post*10^-6*10^-3*12,start=c(1,1,13),count=c(-1,-1,12))




#-- Ocean

ocean_mats = apply(dat_out$out_ocean,c(1,2,3,4),sum)

ocean_prior = ocean_mats[,,,2]

#-- This is effective ocean sink adjustment

ocean_post = ocean_mats[,,,2] + ocean_mats[,,,1]

ncvar_put(con,varid="ocean_flux",vals=ocean_post*10^-3*12,start=c(1,1,1),count=c(-1,-1,22))


#--  Net flux prior and post

net_flux_post = ocean_post*10^-3*12 + abind(nee_post,nee_post[,,1:10])*10^-6*10^-3*12

net_flux_prior = ocean_prior*10^-3*12 + abind(nee_prior,nee_prior[,,1:10])*10^-6*10^-3*12

ncvar_put(con,varid="net_flux_prior",vals=net_flux_prior,start=c(1,1,1),count=c(-1,-1,22))

ncvar_put(con,varid="net_flux",vals=net_flux_post,start=c(1,1,1),count=c(-1,-1,22))


nc_close(con)

}