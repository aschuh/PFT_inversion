##################################################
#--  Inversion script for PFT-based Harmonic scheme
##################################################


read.input.geos <- function(filename) {
  return(readLines(filename))
}

write.input.geos <- function(input.geos.data,
                             filename) {
  writeLines(input.geos.data,con=filename)
}

# remove trailing white space
trim <- function(x) {
    return(gsub("[ \t]*$","",x))
  }

rtrim <- trim

# remove leading white space
ltrim <- function(x) {
    return(gsub("^[ \t]*","",x))
  }


locate.key.input.geos <- function(input.geos.data,key,domain=NULL) {

  # "domain" is a synonym for "menu".  To look for "Output file name"
  # within the ObsPack menu (which is capitalized within the
  # input.geos file), you should use domain="OBSPACK MENU".
  
  # parentheses cause problems for grep, we need to escape them
  key <- gsub('\\(','\\\\(',gsub('\\)','\\\\)',key))

  if(is.null(domain)) {
    #lx <- grep(sprintf("^%s[[:blank:]]*:",key),input.geos.data)
    lx < - grep(paste("*",key,"*",sep=""),input.geos.data)
  } else {
    lx.dstart <- grep(sprintf('^%%%%%%[[:blank:]]*%s[[:blank:]]*%%%%%%[[:blank:]]:',domain),input.geos.data)
    if(length(lx.dstart) == 0) {
      stop(sprintf("No such domain \"%s\" in input.geos data (\"input.geos.data\" here).",domain))
    }
    if(length(lx.dstart) > 1) {
      stop(sprintf("%d matches to domain \"%s\" in input.geos data (\"input.geos.data\" here).",length(lx.dstart),domain))
    }
    # Now we have the location of the domain start.  Search forward
    # for domain end (the next "%%% * MENU %%%:" text or end-of-file)
    nlines.tot <- length(input.geos.data)
    lx.dend <- lx.dstart+1
    while(TRUE) {
    	#print(grepl("^%%% .* MENU %%%.*:$",input.geos.data[lx.dend:lx.dend]))
      #if(grepl("^%%% .* MENU %%%.*:.*$",input.geos.data[lx.dend:lx.dend])) {
      if(grepl("*MENU*:*",input.geos.data[lx.dend:lx.dend])) {  	
        break
      }
      lx.dend <- lx.dend+1
    }
    #lx <- lx.dstart-1+grep(sprintf("^%s[[:blank:]]*:",key),input.geos.data[lx.dstart:lx.dend])
    lx <- lx.dstart-1+grep(paste("*",key,"*",sep=""),input.geos.data[lx.dstart:lx.dend])
    
  }
  if(length(lx)==0) {
    stop(sprintf("No such key \"%s\" in input.geos data (\"input.geos.data\" here).",key))
  }

  if(length(lx)>1) {
    stop(sprintf("%d matches to key \"%s\" in input.geos data (\"input.geos.data\" here).",
                 length(lx),key))
  }


  return(lx)
}

set.value.input.geos <- function(input.geos.data,key,val) {
  lx <- locate.key.input.geos(input.geos.data,key)

  # parentheses cause problems for grep, we need to escape them
  key <- gsub('\\(','\\\\(',gsub('\\)','\\\\)',key))

  # this is fancier that it needs to be, to preserve spacing before colons in input.geos 
  input.geos.data[lx] <- paste(gsub(sprintf("(^%s[[:blank:]]*:).*",key),"\\1",input.geos.data[lx]),val)
  return(input.geos.data)
}

set.value.input.geos2 <- function(input.geos.data,key,val) {

   # parentheses cause problems for grep, we need to escape them
  key <- gsub('\\(','\\\\(',gsub('\\)','\\\\)',key))


  #lx <- grep(input.geos.data,key)
  matches = sum(grepl(glob2rx(key),input.geos.data))
  if(matches==1)
  {
  lx = grep(glob2rx(key),input.geos.data)
  }else{
  	print(paste("problem w/ keyword ",matches," matches",sep=""))
  	stop()
  	}
  	
  input.geos.data[lx] <- paste(strsplit(input.geos.data[lx],":")[[1]][1],":  ",val,sep="")
  return(input.geos.data)
}

get.value.input.geos <- function(input.geos.data,key,domain=NULL) {
  lx <- locate.key.input.geos(input.geos.data,key,domain=domain)

  # parentheses cause problems for grep, we need to escape them
  key <- gsub('\\(','\\\\(',gsub('\\)','\\\\)',key))

  return(trim(ltrim(gsub(sprintf("^%s[[:blank:]]*:",key),"",input.geos.data[lx]))))
}
 

dir.exists2 = function (paths) 
{
    x = base::file.info(paths)$isdir
    !is.na(x) & x
}

create_dirs = function(baserundir,baseoutdir)
{
	
    mk_subdir=function(dir,x){
	  if(!dir.exists2(paste(dir,x,sep="/"))){dir.create(paste(dir,x,sep="/"))}
    }
	
	mk_subdir(baserundir,"")
	mk_subdir(baserundir,"land")
	mk_subdir(baserundir,"ocean")
	mk_subdir(baserundir,"fossil")
	mk_subdir(baserundir,"bg")
	mk_subdir(baserundir,"fires")
		
	subdirs = apply(expand.grid(c("land","ocean","fossil","fires","bg"),
	              c("CO2","oco2","obspack","input.geoses","hemco.configs","restarts","logs","betas")),
	                 1,FUN=function(x){paste(x,collapse="/")})			
	
	sapply(subdirs,FUN=function(x){mk_subdir(baserundir,x)})
}

replace_file_string = function(text_vector,line_num,file)
{
	text_vector[line_num]= paste(strsplit(singles_file[line_num],":")[[1]][1],": ",file,sep="")	
}

fill.pft.input.hemco=function(gppresp=0,sincos=0,harmonic=0,pft=1,tr=1,dates,input.hemco.file)
{
	input.data = read.input.geos(input.hemco.file)
	input.data = set.value.input.geos2(input.geos.data=input.data,key="GPPRESP",val=gppresp) 
     input.data = set.value.input.geos2(input.geos.data=input.data,key="SINCOS",val=sincos)
     set.value.input.geos2(input.data,key="PFT",val=pft)
     set.value.input.geos2(input.data,key="TR",val=tr)
     set.value.input.geos2(input.data,key="HARMONIC",val=harmonic)
     return(input.data)
               }

fill.pft.input.geos=function(gppresp=0,sincos=0,harmonic=0,pft=1,tr=1,input.geos.file)
{
	input.data = read.input.geos(input.hemco.file)
	input.data =set.value.input.geos2(input.data,key="Start YYYYMMDD, HHMMSS",gppresp)
     input.data =set.value.input.geos2(input.data,key="End   YYYYMMDD, HHMMSS",sincos)
     input.data =set.value.input.geos2(input.data,key="Input restart file",pft)
     input.data =set.value.input.geos2(input.data,key="MERRA2",tr)
     input.data =set.value.input.geos2(input.data,key="Run directory",harmonic)
     
     return(input.data)
               }

create_input_geoses = function(baserundir="/discover/nobackup/aschuh/run.v11_oco2_080117",
                   baseoutdir = "/discover/nobackup/aschuh/run.v11_oco2_080117/",
           template_input_geos="/discover/nobackup/aschuh/inversion_template_files/input.geoses/input.geos.template")
{
	#-- FOR LAND
     combos.grid = expand.grid(gppresp=0:1,sincos=0:1,harmonic=0:3,pft=1:25,tr=1:11,year=2014:2016)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     input.data.orig = readLines(template_input_geos)
   	 fluxtype = "land"

    	input.data =set.value.input.geos2(input.geos.data=input.data.orig,key="*Run*directory*",paste(baserundir,"/",fluxtype,sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*=>*GEOS-FP*","../MERRA/MERRA.output/GEOS_4x5.d/GEOS_FP/YYYY/MM/")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Transport*Timestep*","30")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Convect*Timestep*","30")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Sat*track*info*","/discover/nobackup/aschuh/run/lite_wbias_20170815_standardized/inv_weight_insitu/oco2.MMDDYYYY.4x5.nc")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*OBSPACK*infile*","/discover/nobackup/aschuh/run/obspack_co2_1_NRT_v3.3_prepped_inputs_2017-04-26_links/flask_input.YYYYMMDD00.nc")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Avg*timeser*file*",paste(baseoutdir,fluxtype,"/CO2/ts_1h_avg.YYYYMMDDhh",sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Input*restart*file*","/discover/nobackup/aschuh/inversion_template_files/CO2/restart.2014090600.nc")
    	
      #-- Going to need to add options for 4x5 vs 2x25 including i/j/l indicies out ND50	
      
      for(i in 1:(dim(combos.grid)[1]))
      {
 
    	  combo_tag = paste(combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep="")
    	
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Sat*output*file*name*",paste(baseoutdir,fluxtype,"/oco2/oco2.YYYYMMDD_",combo_tag,".nc",sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*OBSPACK*outfile*",paste(baseoutdir,fluxtype,"/obspack/obspack.YYYYMMDD00_",combo_tag,".nc",sep=""))
     	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Start*YYYYMMDD*HHMMSS*","20140901 000000")
       	input.data =set.value.input.geos2(input.geos.data=input.data,key="*End*YYYYMMDD*HHMMSS*","20170101 000000")
      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*HEMCO*Input*file*",
      	paste(baserundir,"/",fluxtype,"/hemco.configs/HEMCO_Config_",combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep=""))
      	
      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Output*restart*file*",paste(baseoutdir,fluxtype,"/restarts/GEOSChem_restart.",combo_tag,sep=""))
    	
       print(paste("writing ",i," of ",dim(combos.grid)[1]," input.geoses for land..."))
       write.input.geos(input.geos.data=input.data,
                             paste(baserundir,"/",fluxtype,"/input.geoses/input.geos_",combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep=""))
       }
    

	#-- FOR OCEAN
     combos.grid = expand.grid(gppresp=0,sincos=0:1,harmonic=0:2,pft=1,tr=1:30,year=2014:2016)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     input.data.orig = readLines(template_input_geos)
   	fluxtype = "ocean"
   	   	
	input.data =set.value.input.geos2(input.geos.data=input.data.orig,key="*Run*directory*",paste(baserundir,"/",fluxtype,sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*=>*GEOS-FP*","../MERRA/MERRA.output/GEOS_4x5.d/GEOS_FP/YYYY/MM/")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Transport*Timestep*","30")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Convect*Timestep*","30")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Sat*track*info*","/discover/nobackup/aschuh/run/lite_wbias_20170815_standardized/inv_weight_insitu/oco2.MMDDYYYY.4x5.nc")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*OBSPACK*infile*","/discover/nobackup/aschuh/run/obspack_co2_1_NRT_v3.3_prepped_inputs_2017-04-26_links/flask_input.YYYYMMDD00.nc")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Avg*timeser*file*",paste(baseoutdir,fluxtype,"/CO2/ts_1h_avg.YYYYMMDDhh",sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Input*restart*file*","/discover/nobackup/aschuh/inversion_template_files/CO2/restart.2014090600.nc")
    	    	
      for(i in 1:(dim(combos.grid)[1]))
      {
    	  combo_tag = paste(combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep="")
    	
    	    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Sat*output*file*name*",paste(baseoutdir,fluxtype,"/oco2/oco2.YYYYMMDD_",combo_tag,".nc",sep=""))
    	    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*OBSPACK*outfile*",paste(baseoutdir,fluxtype,"/obspack/obspack.YYYYMMDD00_",combo_tag,".nc",sep=""))
    	    	
      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*HEMCO*Input*file*",
      	 paste(baserundir,"/",fluxtype,"/hemco.configs/HEMCO_Config_",combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep=""))
      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Start*YYYYMMDD*HHMMSS*","20140901 000000")
       	input.data =set.value.input.geos2(input.geos.data=input.data,key="*End*YYYYMMDD*HHMMSS*","20170101 000000")

      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Output*restart*file*",paste(baseoutdir,fluxtype,"/restarts/GEOSChem_restart.",combo_tag,sep=""))
      	    	
       print(paste("writing ",i," of ",dim(combos.grid)[1]," input.geoses for ocean..."))
       write.input.geos(input.geos.data=input.data,
                             paste(baserundir,"/",fluxtype,"/input.geoses/input.geos_",combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep=""))
       }
   		
	#-- FOR FIRES
     combos.grid = expand.grid(gppresp=0,sincos=0:1,harmonic=0:2,pft=1,tr=1:30,year=2014:2016)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     input.data.orig = readLines(template_input_geos)
   	fluxtype = "fires"
   	   	
	input.data =set.value.input.geos2(input.geos.data=input.data.orig,key="*Run*directory*",paste(baserundir,"/",fluxtype,sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*=>*GEOS-FP*","../MERRA/MERRA.output/GEOS_4x5.d/GEOS_FP/YYYY/MM/")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Transport*Timestep*","30")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Convect*Timestep*","30")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Sat*track*info*","/discover/nobackup/aschuh/run/lite_wbias_20170815_standardized/inv_weight_insitu/oco2.MMDDYYYY.4x5.nc")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*OBSPACK*infile*","/discover/nobackup/aschuh/run/obspack_co2_1_NRT_v3.3_prepped_inputs_2017-04-26_links/flask_input.YYYYMMDD00.nc")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Avg*timeser*file*",paste(baseoutdir,fluxtype,"/CO2/ts_1h_avg.YYYYMMDDhh",sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Input*restart*file*","/discover/nobackup/aschuh/inversion_template_files/CO2/restart.2014090600.nc")
    	    	
      for(i in 1:(dim(combos.grid)[1]))
      {
    	  combo_tag = paste(combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep="")
    	
    	    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Sat*output*file*name*",paste(baseoutdir,fluxtype,"/oco2/oco2.YYYYMMDD_",combo_tag,".nc",sep=""))
    	    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*OBSPACK*outfile*",paste(baseoutdir,fluxtype,"/obspack/obspack.YYYYMMDD00_",combo_tag,".nc",sep=""))
    	    	
      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*HEMCO*Input*file*",
      	 paste(baserundir,"/",fluxtype,"/hemco.configs/HEMCO_Config_",combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep=""))
      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Start*YYYYMMDD*HHMMSS*","20140901 000000")
       	input.data =set.value.input.geos2(input.geos.data=input.data,key="*End*YYYYMMDD*HHMMSS*","20170101 000000")

      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Output*restart*file*",paste(baseoutdir,fluxtype,"/restarts/GEOSChem_restart.",combo_tag,sep=""))
      	    	
       print(paste("writing ",i," of ",dim(combos.grid)[1]," input.geoses for ocean..."))
       write.input.geos(input.geos.data=input.data,
                             paste(baserundir,"/",fluxtype,"/input.geoses/input.geos_",combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep=""))
       }
       
      #-- FOR FOSSIL
     combos.grid = expand.grid(gppresp=0,sincos=0,harmonic=0,pft=0,tr=0,year=0)
     input.data.orig = readLines(template_input_geos)
   	fluxtype = "fossil"
   	   	
	input.data =set.value.input.geos2(input.geos.data=input.data.orig,key="*Run*directory*",paste(baserundir,"/",fluxtype,sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*=>*GEOS-FP*","../MERRA/MERRA.output/GEOS_4x5.d/GEOS_FP/YYYY/MM/")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Transport*Timestep*","30")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Convect*Timestep*","30")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Sat*track*info*","/discover/nobackup/aschuh/run/lite_wbias_20170815_standardized/inv_weight_insitu/oco2.MMDDYYYY.4x5.nc")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*OBSPACK*infile*","/discover/nobackup/aschuh/run/obspack_co2_1_NRT_v3.3_prepped_inputs_2017-04-26_links/flask_input.YYYYMMDD00.nc")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Avg*timeser*file*",paste(baseoutdir,fluxtype,"/CO2/ts_1h_avg.YYYYMMDDhh",sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Input*restart*file*","/discover/nobackup/aschuh/inversion_template_files/CO2/restart.2014090600.nc")
    	    	
      for(i in 1:(dim(combos.grid)[1]))
      {
    	  combo_tag = paste(combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep="")
    	
    	    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Sat*output*file*name*",paste(baseoutdir,fluxtype,"/oco2/oco2.YYYYMMDD_",combo_tag,".nc",sep=""))
    	    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*OBSPACK*outfile*",paste(baseoutdir,fluxtype,"/obspack/obspack.YYYYMMDD00_",combo_tag,".nc",sep=""))
    	    	
      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*HEMCO*Input*file*",
      	 paste(baserundir,"/",fluxtype,"/hemco.configs/HEMCO_Config_",combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep=""))
      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Start*YYYYMMDD*HHMMSS*","20140901 000000")
       	input.data =set.value.input.geos2(input.geos.data=input.data,key="*End*YYYYMMDD*HHMMSS*","20170101 000000")

      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Output*restart*file*",paste(baseoutdir,fluxtype,"/restarts/GEOSChem_restart.",combo_tag,sep=""))
      	    	
       print(paste("writing ",i," of ",dim(combos.grid)[1]," input.geoses for fossil..."))
       write.input.geos(input.geos.data=input.data,
                             paste(baserundir,"/",fluxtype,"/input.geoses/input.geos_",combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep=""))
       }
       
   	#-- FOR BG
     combos.grid = expand.grid(gppresp=0,sincos=0,harmonic=0,pft=0,tr=0,year=0)
     input.data.orig = readLines(template_input_geos)
   	 fluxtype = "bg"
   	   	
	 input.data =set.value.input.geos2(input.geos.data=input.data.orig,key="*Run*directory*",paste(baserundir,"/",fluxtype,sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*=>*GEOS-FP*","../MERRA/MERRA.output/GEOS_4x5.d/GEOS_FP/YYYY/MM/")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Transport*Timestep*","30")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Convect*Timestep*","30")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Sat*track*info*","/discover/nobackup/aschuh/run/lite_wbias_20170815_standardized/inv_weight_insitu/oco2.MMDDYYYY.4x5.nc")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*OBSPACK*infile*","/discover/nobackup/aschuh/run/obspack_co2_1_NRT_v3.3_prepped_inputs_2017-04-26_links/flask_input.YYYYMMDD00.nc")
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Avg*timeser*file*",paste(baseoutdir,fluxtype,"/CO2/ts_1h_avg.YYYYMMDDhh",sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Input*restart*file*","/discover/nobackup/aschuh/inversion_template_files/CO2/restart.2014090600.nc.REAL")
    	    	

      for(i in 1:(dim(combos.grid)[1]))
      {
    	combo_tag = paste(combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep="")
    	
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Sat*output*file*name*",paste(baseoutdir,fluxtype,"/oco2/oco2.YYYYMMDD_",combo_tag,".nc",sep=""))
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*OBSPACK*outfile*",paste(baseoutdir,fluxtype,"/obspack/obspack.YYYYMMDD00_",combo_tag,".nc",sep=""))
    	    	
      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*HEMCO*Input*file*",
      	 paste(baserundir,"/",fluxtype,"/hemco.configs/HEMCO_Config_",combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep=""))
      	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Start*YYYYMMDD*HHMMSS*","20140901 000000")
       	input.data =set.value.input.geos2(input.geos.data=input.data,key="*End*YYYYMMDD*HHMMSS*","20170101 000000")
    	
       	input.data =set.value.input.geos2(input.geos.data=input.data,key="*Output*restart*file*",paste(baseoutdir,fluxtype,"/restarts/GEOSChem_restart.",combo_tag,sep=""))
      	    	
      	    	
       print(paste("writing ",i," of ",dim(combos.grid)[1]," input.geoses for fossil..."))
       write.input.geos(input.geos.data=input.data,
                             paste(baserundir,"/",fluxtype,"/input.geoses/input.geos_",combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep=""))
       }
                 		
}

create_hemco_configs = function(baserundir,template_hemco_config="/discover/nobackup/aschuh/inversion_template_files/hemco.configs/HEMCO_Config.rc_template")
{
	
#-- FOR LAND
     combos.grid = expand.grid(gppresp=0:1,sincos=0:1,harmonic=0:3,pft=1:25,tr=1:11,year=2014:2016)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     input.data.orig = readLines(template_hemco_config)
   	 fluxtype = "land"
   	
    for(i in 1:(dim(combos.grid)[1]))
      {
      	
     combo_tag = paste(combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep="")

     input.data =set.value.input.geos2(input.geos.data=input.data.orig,key="*TR:*",combos.grid[i,5])       	
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*PFT:*",combos.grid[i,4])
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*SINCOS:*",combos.grid[i,2])
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*HARMONIC:*",combos.grid[i,3])
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*GPPRESP:*",combos.grid[i,1])

     RESP_TF = combos.grid[i,1] == 0
     GPP_TF = ! RESP_TF

     #-- If HARMONIC ARG is >= 1, then we modulate, otherwise we run base prior
     HARMONIC_TF = combos.grid[i,3] != 0

     if(HARMONIC_TF == TRUE)
       {
     	harm_arg_gpp = "200"
     	harm_arg_resp = "100"
     	}else{
     	harm_arg_gpp = "-"
     	harm_arg_resp = "-"
     	}
     

     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*SIB4_RESP_CO2*",tolower(as.character(RESP_TF)))
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*SIB4_GPP_CO2*",tolower(as.character(GPP_TF)))     
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*CT_OCN_2016*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*CT_FIRE_2016*","false")     
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*BASU_FOSSIL_2016*","false")

     yr_run = combos.grid[i,6]
          
     input.data[grep(glob2rx("*sib4_fluxes*ASSIM*"),input.data)] = paste("0 CO2 /discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4/"
     ,yr_run,"/sib4_fluxes_TR_$TR_$YYYY$MM$DD.nc ASSIM ",(yr_run-1),"-",(yr_run),"/1-12/1-31/0-23 R xy+\"PFT=$PFT\"  umol/m2/s CO2 ",harm_arg_gpp," 1 1",sep="")

     input.data[grep(glob2rx("*sib4_fluxes*RESP.BAL*"),input.data)] = paste("0 CO2 /discover/nobackup/aschuh/data/sib_harmonic/hemco_sib4/"
     ,yr_run,"/sib4_fluxes_TR_$TR_$YYYY$MM$DD.nc RESP.BAL ",(yr_run-1),"-",(yr_run),"/1-12/1-31/0-23 R xy+\"PFT=$PFT\"  umol/m2/s CO2 ",harm_arg_resp," 1 1",sep="")

     input.data =set.value.input.geos2(input.geos.data=input.data,key="*Logfile:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO.log',combo_tag,sep=""))

     input.data =set.value.input.geos2(input.geos.data=input.data,key="*DiagnPrefix:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO_diagnostics',combo_tag,sep=""))
                           	    	        	
     print(paste("writing ",i," of ",dim(combos.grid)[1]," hemco.configs for land..."))
     
     write.input.geos(input.geos.data=input.data,
                             paste(baserundir,"/",fluxtype,"/hemco.configs/HEMCO_Config_",combo_tag,sep=""))
    }
    

	#-- FOR OCEAN
     combos.grid = expand.grid(gppresp=0,sincos=0:1,harmonic=0:2,pft=1,tr=1:30,year=2014:2016)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     input.data.orig = readLines(template_hemco_config)
   	fluxtype = "ocean"
   	   	
	 for(i in 1:(dim(combos.grid)[1]))
      {
           	
    	  combo_tag = paste(combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep="")
    	  
    	input.data =set.value.input.geos2(input.geos.data=input.data.orig,key="*TR:*",combos.grid[i,5])    	
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*PFT:*",combos.grid[i,4])
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*SINCOS:*",combos.grid[i,2])
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*HARMONIC:*",combos.grid[i,3])
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*GPPRESP:*",combos.grid[i,1])


     #-- If HARMONIC ARG is >= 1, then we modulate, otherwise we run base prior
     HARMONIC_TF = combos.grid[i,3] != 0

     if(HARMONIC_TF == TRUE)
       {
     	harm_arg_ocn = "300"
     	}else{
     	harm_arg_ocn = "-"
     	}
     	
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*SIB4_RESP_CO2*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*SIB4_GPP_CO2*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*CT_OCN_2016*","true")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*CT_FIRE_2016*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*BASU_FOSSIL_2016*","false")

     yr_run = combos.grid[i,6]
               
     ocn_mask_ext = as.character(2000+combos.grid[i,5])          
               
     input.data[grep(glob2rx("*ocn_flux_opt*"),input.data)] = paste("0 CO2  /discover/nobackup/aschuh/data/ct2016_optimized/hemco_data/climate/",yr_run,"/CT2016.flux1x1.climate_oceansub.$YYYY$MM$DD.nc ocean_fluxes ",(yr_run-1),"-",(yr_run),"/1-12/1-31/0-23/+90minutes R xy+\"ocean_region=$TR\"  mol/m2/s  CO2  ",harm_arg_ocn," 1 1",sep="")
            	    	        	
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*Logfile:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO.log',combo_tag,sep=""))
          
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*DiagnPrefix:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO_diagnostics',combo_tag,sep=""))
          
     print(paste("writing ",i," of ",dim(combos.grid)[1]," hemco.configs for ocean..."))
     
     write.input.geos(input.geos.data=input.data,
                             paste(baserundir,"/",fluxtype,"/hemco.configs/HEMCO_Config_",combo_tag,sep=""))
    }

	#-- FOR FIRES
     combos.grid = expand.grid(gppresp=0,sincos=0,harmonic=0,pft=0,tr=0,year=0)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     input.data.orig = readLines(template_hemco_config)
   	fluxtype = "fires"
   	   	
	 for(i in 1:(dim(combos.grid)[1]))
      {
           	
    	  combo_tag = paste(combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep="")
    	  
    	input.data =set.value.input.geos2(input.geos.data=input.data.orig,key="*TR:*",combos.grid[i,5])    	
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*PFT:*",combos.grid[i,4])
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*SINCOS:*",combos.grid[i,2])
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*HARMONIC:*",combos.grid[i,3])
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*GPPRESP:*",combos.grid[i,1])

     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*SIB4_RESP_CO2*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*SIB4_GPP_CO2*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*CT_OCN_2016*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*CT_FIRE_2016*","true")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*BASU_FOSSIL_2016*","false")
               
     input.data[grep(glob2rx("*ocn_flux_opt*"),input.data)] = paste("0 CO2  /discover/nobackup/aschuh/data/ct2016_optimized/hemco_data/CT2016.flux1x1.$YYYY$MM$DD.nc ocn_flux_opt ",combos.grid[i,6],"/1-12/1-31/0-23/+90minutes C xy  mol/m2/s  CO2  - 1 1",sep="")
            	    	        	
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*Logfile:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO.log',combo_tag,sep=""))
          
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*DiagnPrefix:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO_diagnostics',combo_tag,sep=""))
          
     print(paste("writing ",i," of ",dim(combos.grid)[1]," hemco.configs for fires..."))
     
     write.input.geos(input.geos.data=input.data,
                             paste(baserundir,"/",fluxtype,"/hemco.configs/HEMCO_Config_",combo_tag,sep=""))
    }



	#-- FOR FOSSIL
     combos.grid = expand.grid(gppresp=0,sincos=0,harmonic=0,pft=0,tr=0,year=0)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     input.data.orig = readLines(template_hemco_config)
   	fluxtype = "fossil"
   	   	
	 for(i in 1:(dim(combos.grid)[1]))
      {
      	
    	  combo_tag = paste(combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep="")
    	  
    	input.data =set.value.input.geos2(input.geos.data=input.data.orig,key="*TR:*",combos.grid[i,5])    	
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*PFT:*",combos.grid[i,4])
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*SINCOS:*",combos.grid[i,2])
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*HARMONIC:*",combos.grid[i,3])
    	input.data =set.value.input.geos2(input.geos.data=input.data,key="*GPPRESP:*",combos.grid[i,1])

     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*SIB4_RESP_CO2*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*SIB4_GPP_CO2*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*CT_OCN_2016*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*CT_FIRE_2016*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*BASU_FOSSIL_2016*","true")
     
     
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*Logfile:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO.log',combo_tag,sep=""))               
            	    	        	
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*Logfile:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO.log',combo_tag,sep=""))
          
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*DiagnPrefix:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO_diagnostics',combo_tag,sep=""))
          
     print(paste("writing ",i," of ",dim(combos.grid)[1]," hemco.configs for fossil..."))
     
     write.input.geos(input.geos.data=input.data,
                             paste(baserundir,"/",fluxtype,"/hemco.configs/HEMCO_Config_",combo_tag,sep=""))
    }
 
	#-- FOR BG
     combos.grid = expand.grid(gppresp=0,sincos=0,harmonic=0,pft=0,tr=0,year=0)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     input.data.orig = readLines(template_hemco_config)
   	 fluxtype = "bg"
   	   	
	 for(i in 1:(dim(combos.grid)[1]))
      {
      	
     combo_tag = paste(combos.grid[i,6],"_",combos.grid[i,1],"_",combos.grid[i,2],"_",combos.grid[i,3],"_",combos.grid[i,4],"_",combos.grid[i,5],sep="")
    	  
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*SIB4_RESP_CO2*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*SIB4_GPP_CO2*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*CT_OCN_2016*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*CT_FIRE_2016*","false")
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*-->*BASU_FOSSIL_2016*","false")
     
     
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*Logfile:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO.log',combo_tag,sep=""))               
            	    	        	
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*Logfile:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO.log',combo_tag,sep=""))
          
     input.data =set.value.input.geos2(input.geos.data=input.data,key="*DiagnPrefix:*",
          paste(baserundir,'/',fluxtype,'/logs/HEMCO_diagnostics',combo_tag,sep=""))
          
     print(paste("writing ",i," of ",dim(combos.grid)[1]," hemco.configs for fossil..."))
     
     write.input.geos(input.geos.data=input.data,
                             paste(baserundir,"/",fluxtype,"/hemco.configs/HEMCO_Config_",combo_tag,sep=""))
    }       	
}

create_launch_scripts = function(baserundir,run.pulse.file="/discover/nobackup/aschuh/inversion_template_files/run.pulse")
{
	 pulse  =  readLines(run.pulse.file)

     lx = grep(glob2rx("*basedir=*"),pulse)

     pulse[lx] = paste("basedir='",baserundir,"'",sep="")
     
	 writeLines(pulse,paste(baserundir,"/run.pulse",sep=""))

     system(paste("ln -s /discover/nobackup/aschuh/Code.v11-01/bin/geos ",baserundir,"/",sep=""))
	 
   #/discover/nobackup/aschuh/run.v9-02_oco2_111716/land/run.pulse 20140906 1 365 0 0 0 01 01 > out_0_0_0_01_01.txt

   	 line_vector1 = vector()
   	 line_vector2 = vector()
   	 line_vector3 = vector()   	 
   	 line_vector4 = vector()
   	 line_vector5 = vector()
   	 
     #-- FOR LAND
     combos.grid = expand.grid(gppresp=0:1,sincos=0:1,harmonic=0:3,pft=1:25,tr=1:11,year=2014:2016)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
   	 fluxtype = "land"
   	   	    	 
    for(i in 1:(dim(combos.grid)[1]))
      {
    	combo_tag = paste(combos.grid[i,6]," ",combos.grid[i,1]," ",combos.grid[i,2]," ",combos.grid[i,3],
    	        " ",combos.grid[i,4]," ",combos.grid[i,5],sep="")
    	        
    	line_vector1[i] = paste(baserundir,"/run.pulse ",fluxtype," ",combo_tag," > out_",gsub(" ","_",combo_tag),"_txt",sep="")
    	
      }

      #-- FOR OCEAN
     combos.grid = expand.grid(gppresp=0,sincos=0:1,harmonic=0:2,pft=1,tr=1:30,year=2014:2016)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     fluxtype = "ocean"

     for(i in 1:(dim(combos.grid)[1]))
      {
    	combo_tag = paste(combos.grid[i,6]," ",combos.grid[i,1]," ",combos.grid[i,2]," ",combos.grid[i,3],
    	        " ",combos.grid[i,4]," ",combos.grid[i,5],sep="")
    	        
    	line_vector2[i] = paste(baserundir,"/run.pulse ",fluxtype," ",combo_tag," > out_",gsub(" ","_",combo_tag),"_txt",sep="")
    	
      }
      
	#-- FOR FOSSIL
     combos.grid = expand.grid(gppresp=0,sincos=0,harmonic=0,pft=0,tr=0,year=0)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     fluxtype = "fossil"
   	   	
	 for(i in 1:(dim(combos.grid)[1]))
      {
    	combo_tag = paste(combos.grid[i,6]," ",combos.grid[i,1]," ",combos.grid[i,2]," ",combos.grid[i,3],
    	        " ",combos.grid[i,4]," ",combos.grid[i,5],sep="")
    	        
    	line_vector3[i] = paste(baserundir,"/run.pulse ",fluxtype," ",combo_tag," > out_",gsub(" ","_",combo_tag),"_txt",sep="")
    	
      }      	
      
      
	#-- FOR FIRES
     combos.grid = expand.grid(gppresp=0,sincos=0,harmonic=0,pft=0,tr=0,year=0)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     fluxtype = "fires"
   	   	
	 for(i in 1:(dim(combos.grid)[1]))
      {
    	combo_tag = paste(combos.grid[i,6]," ",combos.grid[i,1]," ",combos.grid[i,2]," ",combos.grid[i,3],
    	        " ",combos.grid[i,4]," ",combos.grid[i,5],sep="")
    	        
    	line_vector4[i] = paste(baserundir,"/run.pulse ",fluxtype," ",combo_tag," > out_",gsub(" ","_",combo_tag),"_txt",sep="")
    	
      } 
      
      
	#-- FOR BG
     combos.grid = expand.grid(gppresp=0,sincos=0,harmonic=0,pft=0,tr=0,year=0)
     combos.grid = combos.grid[!(combos.grid[,3] == 0 & combos.grid[,2]==1),]
     fluxtype = "bg"
   	   	
	 for(i in 1:(dim(combos.grid)[1]))
      {
    	combo_tag = paste(combos.grid[i,6]," ",combos.grid[i,1]," ",combos.grid[i,2]," ",combos.grid[i,3],
    	        " ",combos.grid[i,4]," ",combos.grid[i,5],sep="")
    	        
    	line_vector5[i] = paste(baserundir,"/run.pulse ",fluxtype," ",combo_tag," > out_",gsub(" ","_",combo_tag),"_txt",sep="")
    	
      }       
      
     writeLines(c(line_vector1,line_vector2,line_vector3,line_vector4,line_vector5),con=paste(baserundir,"/exec.script",sep=""))
      
  }

#-- CREATE RERUN FOR LOST SIMS

create_reruns = function(dir,filename_out){ 
	fls.inputgeos = list.files(dir,pattern="new.exec.script_",full.names=TRUE)
	for(i in 1:length(fls.inputgeos))
	{
		if(i==1) { indat = readLines(fls.inputgeos[i]) }else{
			indat = c(indat,readLines(fls.inputgeos[i]))
		}
	}
	
    ens.attemptrun = sapply(indat,FUN=function(x){gsub("_txt","",gsub("out_","",strsplit(x," > ")[[1]][2]))})
    
    type_attempt = sapply(indat,FUN=function(x){strsplit(x," ")[[1]][2]})

    ran = paste(type_attempt,ens.attemptrun,sep="_")


	fls.output.land = list.files(paste(dir,"/land/oco2/",sep=""),pattern="oco2_total_")
	fls.output.ocean = list.files(paste(dir,"/ocean/oco2/",sep=""),pattern="oco2_total_")
	fls.output.fossil = list.files(paste(dir,"/fossil/oco2/",sep=""),pattern="oco2_total_")
	fls.output.fires = list.files(paste(dir,"/fires/oco2/",sep=""),pattern="oco2_total_")
	fls.output.bg = list.files(paste(dir,"/bg/oco2/",sep=""),pattern="oco2_total_")
			
	fls.output = c(fls.output.land,fls.output.ocean,fls.output.fossil,
	         fls.output.fires,fls.output.bg)
	        
		
	type_ran = c(rep("land",length(fls.output.land)),
	            rep("ocean",length(fls.output.ocean)),
	            rep("fossil",length(fls.output.fossil)),
	            rep("fires",length(fls.output.fires)),
	            rep("bg",length(fls.output.bg)))
	
	ens.run = sapply(fls.output,FUN=function(x){gsub(".nc","",gsub("oco2_total_","",x))})	
	
	ran.done = paste(type_ran,ens.run,sep="_")


	needtorun = ran[!(ran %in% ran.done)]
	
	linestowrite = indat[!(ran %in% ran.done)]
	
	write(linestowrite,file=filename_out)
}

#-- Running

create_dirs(baserundir="/discover/nobackup/aschuh/run.v11_oco2_080117",baseoutdir="/discover/nobackup/aschuh/run.v11_oco2_080117")

create_input_geoses(baserundir="/discover/nobackup/aschuh/run.v11_oco2_080117",
    baseoutdir = "/discover/nobackup/aschuh/run.v11_oco2_080117/",
   template_input_geos="/discover/nobackup/aschuh/inversion_template_files/input.geoses/input.geos.template")

create_hemco_configs("/discover/nobackup/aschuh/run.v11_oco2_080117","/discover/nobackup/aschuh/inversion_template_files/hemco.configs/HEMCO_Config.rc_template")

create_launch_scripts(baserundir="/discover/nobackup/aschuh/run.v11_oco2_080117/",run.pulse.file="/discover/nobackup/aschuh/inversion_template_files/run.pulse")

create_reruns(dir="/discover/nobackup/aschuh/run.v11_oco2_080117",filename="/discover/nobackup/aschuh/run.v11_oco2_080117/dropped_runs")