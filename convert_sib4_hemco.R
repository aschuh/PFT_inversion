################################
#. Word of warning.  This code partitions sib4 gridded 1x1 output into Transcom by PFT files
#. Any sib4 flux associated w/ TranscomReg > 11 (oceans + greeland/antarctica) is "thrown" away, 
#. probably about 1-2% of the fluxes.  Shouldn't affect balance too much since we're "close" to balanced
#. grid cell to grid cell (not exact because of crop redistibution)
#. Since we run the inversion on the processed TR x PFT files, those are the basis for inversion and 
#. posterior fluxes, should be fine. 
################################

library(plyr)
library(ncdf4)

#template = nc_open("/projects/sandbox/sib4_hemco_rewrite/sib4_fluxes/template.nc")
template = nc_open("/discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4/template.nc")

#transfil = nc_open("/projects/sandbox/iregions.nc")
transfil = nc_open("/discover/nobackup/aschuh/data/misc/iregions.nc")
tregions = ncvar_get(transfil,"transcom_regions")
oregions = ncvar_get(transfil,"ocean_regions")
nc_close(transfil)


fls = list.files("/discover/nobackup/aschuh/data/sib_harmonic",pattern="sib4_flux",full.names=TRUE)
fls.short = list.files("/discover/nobackup/aschuh/data/sib_harmonic",pattern="sib4_flux",full.names=FALSE)

labs = as.vector(unlist(sapply(fls.short,FUN=function(x){strsplit(strsplit(x,"_")[[1]][3],"\\.")[[1]][1]})))

#for(i in 1:length(fls))
library(abind)
library(foreach)
library(doSNOW)

registerDoSNOW(makeCluster(14, type = "SOCK"))

#foreach(i in 2:3)
foreach(i=1:length(fls),.packages=c("ncdf4","plyr","abind"))  %dopar%
{
    	print(paste("Working on file :",fls[i]))
	
	#ncnew = nc_create("/projects/sandbox/sib4_hemco_rewrite/sib4_fluxes/test.nc",vars=varASSIM,force_v4=TRUE,verbose=FALSE)
	
	dat = load.ncdf(fls[i])
	
	for(j in 1:24)
	{
	  print(paste("Working on PFT:",j))
		
	  tmpgpp = dat$assim 
	  tmpresp = dat$resp.bal 
	  tmpgpp[dat$pft.ref != j]	= 0
	  tmpresp[dat$pft.ref != j]	= 0		  
	  		
	  tmp.gpp = aaply(tmpgpp[,,,],.margins=c(4),
       .fun=function(x)
       {
     	return(x*dat$area)
       })

	  tmp.resp = aaply(tmpresp[,,,],.margins=c(4),
       .fun=function(x)
       {
     	return(x*dat$area)
       })
       	  
      tmp.gpp = aperm(apply(tmp.gpp,c(1,2,3),sum),c(2,3,1))
      
      tmp.resp = aperm(apply(tmp.resp,c(1,2,3),sum),c(2,3,1))


      full.gpp = array(0,dim=c(360,180,24,11))
      
      full.resp = array(0,dim=c(360,180,24,11))      

      full.gpp = abind(tmp.gpp,tmp.gpp,tmp.gpp,tmp.gpp,tmp.gpp,
                  tmp.gpp,tmp.gpp,tmp.gpp,tmp.gpp,tmp.gpp,tmp.gpp,along=4)

      full.resp = abind(tmp.resp,tmp.resp,tmp.resp,tmp.resp,tmp.resp,tmp.resp,
                   tmp.resp,tmp.resp,tmp.resp,tmp.resp,tmp.resp,along=4)
                          
      for(k in 1:11)
      {
      	print(paste("...working on transcom region:",k))
      	
      	full.gpp[,,,k] = aperm(aaply(full.gpp[,,,k],c(3),
          	   .fun=function(x){x[tregions!=k]=0;x})
          	     ,c(2,3,1)) 
          	     
      	full.resp[,,,k] = aperm(aaply(full.resp[,,,k],c(3),
          	   .fun=function(x){x[tregions!=k]=0;x})
          	     ,c(2,3,1)) 
          	     
      }


#      for(k in 1:11)
#      {
#      	full.gpp[,,,k] = aperm(aaply(tmp.gpp,c(3),.fun=function(x){x[tregions==k]=0;x}),c(2,3,1))
#      	full.resp[,,,k] = aperm(aaply(tmp.resp,c(3),.fun=function(x){x[tregions==k]=0;x}),c(2,3,1))
#      }
      
      if( j ==1){
           pft.dim = ncdim_def( "PFT", "count", 1:25, unlim=FALSE)
      
           tregion.dim = ncdim_def( "TRANSCOM", "count", 1:11, unlim=FALSE)
         
           t <- ncdim_def( "time", "hours since 2014-01-01", 0:23, unlim=TRUE)
          }

      newdim = list(template$dim$lon,template$dim$lat,t,pft.dim)
            
      varASSIM <- ncvar_def(name="ASSIM", 
          units="umol m-2 sec-1", 
          dim=newdim,
          -1, 
		longname="ASSIM Flux", 
		compression=5,
		prec="float")
		
	  varRESP <- ncvar_def(name="RESP.BAL", 
          units="umol m-2 sec-1", 
          dim=newdim,
          -1, 
		longname="RESP balanced flux", 
		compression=5,
		prec="float")
           
      if(j == 1)
      {
      	 ncnew = list()
      	 
      	 for(k in 1:11)
      	 { 
      	  print(paste("...WRITING PFT ",j," for transcom region:",k))   	 
      	  
      	   ncnew[[k]] = nc_create(paste("/discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4/sib4_fluxes_TR_",k,"_",labs[i],".nc",sep=""),
      	   vars=list(varRESP,varASSIM),force_v4=TRUE,verbose=FALSE)
      	   
      	   ncvar_put( ncnew[[k]], varASSIM,full.gpp[,,,k] ,start=c(1,1,1,j),count=c(-1,-1,-1,1))
     
          ncvar_put( ncnew[[k]], varRESP,full.resp[,,,k] ,start=c(1,1,1,j),count=c(-1,-1,-1,1))
          
         }
        }else{
        
        ncnew = list()
              	 
      	 for(k in 1:11)
      	 {
      	 	        
         print(paste("...WRITING PFT ",j," for transcom region:",k))   
         	      	  
      	  ncnew[[k]] = nc_open(paste("/discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4/sib4_fluxes_TR_",k,"_",labs[i],".nc",sep=""),write=TRUE)

      	  ncvar_put( ncnew[[k]], varASSIM,full.gpp[,,,k] ,start=c(1,1,1,j),count=c(-1,-1,-1,1))
     
         ncvar_put( ncnew[[k]], varRESP,full.resp[,,,k] ,start=c(1,1,1,j),count=c(-1,-1,-1,1))
        }
                	 
        }
      
     #ncvar_put( ncnew, varASSIM,tmp.gpp ,start=c(1,1,1,j),count=c(-1,-1,-1,1))
     
     #ncvar_put( ncnew, varRESP,tmp.resp ,start=c(1,1,1,j),count=c(-1,-1,-1,1))
     
	}
	
	#-- Add CROPS separately
	 
	 #require(abind)
	  
	 tmp.crop = dat$resp.crop

    zeros = array(0,dim=c(360,180,24))

    for(k in 1:11)
    {
       tc.tmp.crop = tmp.crop

       tc.tmp.crop[tregions!=k] = 0
       
       ncvar_put( ncnew[[k]], varASSIM,zeros,start=c(1,1,1,25),count=c(-1,-1,-1,1))
     
       ncvar_put( ncnew[[k]], varRESP,tc.tmp.crop ,start=c(1,1,1,25),count=c(-1,-1,-1,1))
      
       nc_close(ncnew[[k]])
    }
   }

#---

fls = list.files("/discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4",full.names=TRUE,pattern="_2016")

fls2 = list.files("/discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4",full.names=TRUE,pattern="_2017")

flss = c(fls,fls2)

for(i in 1:length(flss))
{
	require(ncdf4)
	print(paste("working on ",i," of ",length(flss),sep=""))
	ff = nc_open(flss[i])
	if(sum(names(ff$dim)=="Time") > 0)
	{
	  nc_close(ff)
	  system(paste("/usr/local/other/SLES11.1/nco/4.4.4/intel-12.1.0.233/bin/ncrename -d Time,time ",flss[i]))
	  ff = nc_open(flss[i])
	 }
	 if(!ncatt_get(ff,"time",attname="calendar")$hasatt)
	 {
	  system(paste("/usr/local/other/SLES11.1/nco/4.4.4/intel-12.1.0.233/bin/ncatted -a calendar,time,o,c,'standard'  ",flss[i],sep=""))
	  }
	 nc_close(ff)
}


###########################
#-- I have to add symlinks to "fake" Sept 16 - Jan 2017 flux data
###########################

unq_dates= gsub("-","",seq(as.Date("2016/9/1"), as.Date("2017/1/1"), "days"))

new_dates = sapply(unq_dates,FUN=function(x){paste(as.numeric(as.character(substring(x,1,4)))-1,substring(x,5,6),substring(x,7,8),sep="")})

for(i in 1:length(unq_dates))
{	
	cat(i)
	for(j in 1:11)
	{
	 if(substring(unq_dates[i],1,4)==2016){yr.target = 2015;yr.link = 2016}
	 if(substring(unq_dates[i],1,4)==2017){yr.target = 2015;yr.link = 2017}
	 	 
	  system(paste("ln -s /discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4/",yr.target,"/sib4_fluxes_TR_",j,"_",new_dates[i],".nc  ",
	           "/discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4/",yr.link,"/sib4_fluxes_TR_",j,"_",unq_dates[i],".nc",sep=""))
	}	
}


############################
#. Create global fluxes
############################

unq_datetags_fils=sapply(sort(seq(as.Date("2014/9/1"), as.Date("2017/3/1"), "days")),FUN=function(x){paste(substring(x,1,4),substring(x,6,7),substring(x,9,10),sep="")})

fls = list.files("/discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4",full.names=TRUE,pattern="TR",recursive=TRUE)

template = load.ncdf("/discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4/template.nc")

#for(i in 1:length(fls))
library(abind)
library(foreach)
library(doSNOW)

registerDoSNOW(makeCluster(28, type = "SOCK"))

CROPS_TF = FALSE

#foreach(i in 2:3)
for(i in 200:length(unq_datetags_fils))

#foreach(i=839:845,.packages=c("ncdf4","plyr","abind"))  %dopar%

#foreach(i=469:487,.packages=c("ncdf4","plyr","abind"))  %dopar%
foreach(i=1:length(unq_datetags_fils),.packages=c("ncdf4","plyr","abind"))  %dopar%
{
	fls.tmp = fls[grep(unq_datetags_fils[i],fls)]
	#fls.tmp
	for(j in 1:length(fls.tmp))
	#foreach(j=1:length(fls.tmp)) %dopar%
	{
		fls.tmp
		print(paste("opening",fls.tmp[j]))
		dat = load.ncdf(fls.tmp[j])
		if(CROPS_TF){dat$resp.bal = dat$resp.bal[,,,25];dat$assim = dat$assim[,,,25]}
		
		if(j ==1){resp.bal = dat$resp.bal;assim = dat$assim}else{
			resp.bal = resp.bal + dat$resp.bal;assim = assim + dat$assim
		}
	}
	
	assim = apply(assim,c(1,2,3),sum)
	resp.bal = apply(resp.bal,c(1,2,3),sum)

    outfil = 	paste("/discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4/global/sib4_fluxes_global_",unq_datetags_fils[i],".nc",sep="")	
	      	   
   print(paste("opening",outfil))            
     t <- ncdim_def( "time", "hours since 2014-01-01", 0:23, unlim=TRUE)

      newdim = list(template$dim$lon,template$dim$lat,t)
            
      varASSIM <- ncvar_def(name="ASSIM", 
          units="umol m-2 sec-1", 
          dim=newdim,
          -1, 
		longname="ASSIM Flux", 
		compression=5,
		prec="float")
		
	  varRESP <- ncvar_def(name="RESP.BAL", 
          units="umol m-2 sec-1", 
          dim=newdim,
          -1, 
		longname="RESP balanced flux", 
		compression=5,
		prec="float")

	hand = nc_create(outfil,
      	   vars=list(varRESP,varASSIM),force_v4=TRUE,verbose=FALSE)
      	   		
	   ncvar_put( hand, varASSIM,assim ,start=c(1,1,1),count=c(-1,-1,-1))
     
      ncvar_put( hand, varRESP,resp.bal ,start=c(1,1,1),count=c(-1,-1,-1))
      
     nc_close(hand)
}