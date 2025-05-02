#' Create a controller file for caflood simulations
#'
#' @param Input path to folder containing input data
#' @param run_name optional. Complete name describing de simulation
#' @param run_short_name short name used for outputs
#' @param time_start period in seconds corresponding the begging of the simulation
#' @param time_end period in seconds corresponding the ending of the simulation
#' @param max_DT
#' @param min_DT
#' @param alpha_fraction
#' @param max_iterations
#' @param roughness_global adimensional coefficient representing Manning roughness coefficient
#' @param ignore_WD
#' @param tolerance
#' @param slope_tolerance
#' @param boundary_ele
#' @param dem_path name of the digital elevation model file
#' @param rain_path name of rainfall event file
#' @param output_period
#' @param Raster_WD_Tolerance
#' @param Upstream_Reduction
#'
#'
caflood_controller <- function(
    Input,
    run_name = "This is a CAFLOOD simulation",
    run_short_name,
    time_start = 0,
    time_end = 86400,
    max_DT = 60,
    min_DT = 0.01,
    alpha_fraction = 0.1,
    max_iterations = 1000000,
    roughness_global,
    ignore_WD = 0.0001,
    tolerance = 0.0001,
    slope_tolerance = 0.528,
    boundary_ele,
    dem_path,
    dem_name,
    rain_path,
    flow_path,
    output_period = 300,
    Raster_WD_Tolerance = 0.01,
    Upstream_Reduction = 0.5
    ){

  writeLines(
    c(paste0("Simulation Name			,",run_name ),
      paste0("Short Name (for outputs)			,",run_short_name),
      "Version	   			, 1,0,0",
      "Model Type			, WCA2Dv2",
      paste0("Time Start (seconds)		,", time_start),
      paste0("Time End   (seconds)		,",  time_end),
      paste0("Max DT (seconds)		,",max_DT),
      paste0("Min DT (seconds)		,", min_DT),
      paste0("Update DT (seconds)		,", 60),
      paste0("Alpha (Fraction DT 0.0-1.0)	,",alpha_fraction),
      paste0("Max Iterations			, ",max_iterations),
      paste0("Roughness Global		, ",roughness_global),
      paste0("Ignore WD (meter)		, ", ignore_WD),
      paste0("Tolerance (meter)		, ", tolerance),
      paste0("Slope Tolerance (%)		, ", slope_tolerance),
      paste0("Boundary Ele (Hi/Closed-Lo/Open), ", boundary_ele),
      paste0("Elevation ASCII			, ", dem_path),
      paste0("Rain Event CSV			, ", rain_path),
      paste0("Water Level Event CSV		, "),
      paste0("Inflow Event CSV		, ", flow_path),
      paste0("Time Plot CSV			, ", paste0(dem_name,"_WLpoints.csv,"), paste0(dem_name,"_VELpoints.csv")),
      paste0("Raster Grid CSV			, ", "WDraster.csv,", "VELraster.csv"),
      paste0("Output Console			, ", 'true'),
      paste0("Output Period (s)		, ", output_period),
      paste0("Output Computation Time		, ", 'true'),
      paste0("Check Volumes	   		, ", 'true'),
      paste0("Remove Proc Data (No Pre-Proc)	, ", 'true'),
      paste0("Remove Pre-Proc Data		, ", 'true'),
      paste0("Raster VEL Vector Field		, ", 'true'),
      paste0("Raster WD Tolerance (meter)	, ", Raster_WD_Tolerance),
      paste0("Update Peak Every DT		, ", "false"),
      paste0("Expand Domain                   , ", "false"),
      paste0("Ignore Upstream			, ", 'true'),
      paste0("Upstream Reduction (meter)	, ",Upstream_Reduction)),
    paste0(Input,run_short_name,"_","control.csv")
  )

}


#' Create rainfall file
#' @param Event_name name of event for output file
#' @param timesteps time steps of rainfall in seconds
#' @param rainfall rainfall intensity in mm/hr
#'
#'
#' @examples
make_it_rain <- function(Input,Event_name,timesteps,rainfall){
  writeLines(
    c(paste0("Event Name,", Event_name),
      paste(c("Rain Intensity (mm/hr)",rainfall),sep=',',collapse = ","),
      paste(c("Time Stop (seconds)",timesteps),sep=',',collapse = ","),
      "Area (tlx tly brx bry),"),
    paste0(Input,Event_name,"_","rain",".csv")
  )
}


make_it_flow <- function(Input,Event_name,timesteps,water_level,tlx,tly){

  writeLines(
    c(paste0("Event Name,", Event_name),
      paste(c("Inflow (cumecs)",water_level),sep=',',collapse = ","),
      paste(c("Time Stop (seconds)",timesteps),sep=',',collapse = ","),
      paste(c("Zone (tlx tly w h)", c(tlx,tly,1,1)),sep=',',collapse = ",")),
    paste0(Input,Event_name,"_","WaterLevelBC",".csv")
  )

}



#' Create and execute a bash file to run CADDIES
#'
#' @param caflood_path
#' @param Input
#' @param Output
#' @param controller
#'

make_it_run <- function(caflood_path,Input,Output,controller){

  writeLines(
    c("ECHO OFF",
      "call echo iniciando",
      paste("cd",caflood_path),
      paste(".\\caflood.exe /WCA2D",Input,controller,Output),
      "call echo terminou"),
    paste0(tempdir(),"\\run.bat"))

  shell.exec(paste0(tempdir(),"\\run.bat"))
}




#' Create file of Velocity and Water Level control points and snap them to river network
#'
#' @param Input Folder with input files
#' @param Lon
#' @param Lat
#' @param id
#' @param dem
#' @param Event_name
#' @param period
#'
#' @import sp
#' @import whitebox
#' @import raster
#' @import terra
#'
make_control_points <- function(Input,
                                Lon,Lat,id,
                                dem,
                                dem_name,
                                Event_name,
                                period,
                                stream_network
                                ){


    inputfile <- tempdir()

    dem <- raster::projectRaster(dem,crs="EPSG:31983")
    raster::writeRaster(dem,paste0(inputfile,"\\dem.tif"),overwrite=TRUE)
    raster::writeRaster(dem, paste0(Input,dem_name,".asc"),format = "ascii",overwrite=TRUE)
    streams <- st_read(stream_network)

    st_write(streams,paste0(inputfile,"\\streams.shp"),
             delete_dsn = T,
             delete_layer = T)

    wbt_fill_burn(
      dem = "dem.tif",
      streams = "streams.shp",
      output = "dem_filled.tif",
      wd = inputfile
    )

    wbt_d8_flow_accumulation(input = paste0(inputfile,"\\dem_filled.tif"),
                             output = paste0(inputfile,"\\D8FA.tif"))

    wbt_d8_pointer(dem = paste0(inputfile,"\\dem_filled.tif"),
                   output = paste0(inputfile,"\\D8pointer.tif"))

    ppoints <- data.frame(
      Lon = Lon,
      Lat = Lat
    )

    ppointsSP <- SpatialPoints(ppoints, proj4string = CRS("+init=epsg:31982"))
    ppointsSP <- st_as_sf(ppointsSP)

    st_write(ppointsSP,
             paste0(inputfile,"\\pourpoints.shp"),
             delete_dsn = T,
             delete_layer = T)


    threshold <- mean(values(rast(paste0(inputfile,paste0("\\D8FA.tif")))),na.rm=T)

    wbt_extract_streams(flow_accum = paste0(inputfile,"\\D8FA.tif"),
                        output = paste0(inputfile,"\\raster_streams.tif"),
                        threshold = threshold)

    wbt_jenson_snap_pour_points(pour_pts = paste0(inputfile,"\\pourpoints.shp"),
                                streams = paste0(inputfile,"\\raster_streams.tif"),
                                output = paste0(inputfile,"\\snappedpp.shp"),
                                snap_dist = 50) #careful with this! Know the units of your data

    pp <- shapefile(paste0(inputfile,"\\snappedpp.shp"))

    coords <- coordinates(pp)

    ### Write VEL points
    writeLines(
      c(paste0("Time Plot Name		,", Event_name),
        "Physical Variable	, VEL",
        paste(c("Points Name",id),sep = ",",collapse = ','),
        paste(c("Points X Coo", coords[,1]) ,sep = ",",collapse = ','),
        paste(c("Points Y Coo", coords[,2])  ,sep = ",",collapse = ','),
        paste("Period (seconds)", period,sep = ",")),
      paste0(Input,dem_name,"_","VELpoints",".csv")
    )

    ### Write WL points
    writeLines(
      c(paste0("Time Plot Name		,", Event_name),
        "Physical Variable	, WL",
        paste(c("Points Name",id),sep = ",",collapse = ','),
        paste(c("Points X Coo", coords[,1]) ,sep = ",",collapse = ','),
        paste(c("Points Y Coo", coords[,2])  ,sep = ",",collapse = ','),
        paste("Period (seconds)", period,sep = ",")),
      paste0(Input,dem_name,"_","WLpoints",".csv")
    )


    writeLines(
      c(
        "Raster Grid Name	, Teste LiDAR 0.5 m Water Depth Raster Grid",
        paste("Physical Variable","WD",sep = ","),
        paste("Peak","true",sep = ","),
        paste("Period (seconds)","0",sep = ",")),
      paste0(Input,"WDraster",".csv")
    )

    writeLines(
      c(
        "Raster Grid Name	, Teste LiDAR 0.5 m Water Depth Raster Grid",
        paste("Physical Variable","VEL",sep = ","),
        paste("Peak","true",sep = ","),
        paste("Period (seconds)","0",sep = ",")),
      paste0(Input,"VELraster",".csv")

    )

  }




effective_precipitation <- function(P,CN,lambda,SS=0){
  #CN<-CN+(imperviousness)*(98-CN)
  S = 25400/CN-254
  Ia = lambda*S
  Prec.acc <- cumsum(P)-SS
  Prec.eff <- ifelse(Prec.acc <= Ia,0,(Prec.acc-Ia)^2/(Prec.acc-Ia+S))
  Prec.eff.diff <- c(0,diff(Prec.eff,1))
  return(Prec.eff.diff)
}

cafloodr <- function(param,
                     caflood_path,
                     Input,
                     Output,
                     event_name,
                     dem_name,
                     path2,
                     outlet_path = "C:/Projetos/01_DTI-A/02_Input/01_SP/01_Hidro_hora/estacoes.shp",
                     rain_path = "C:/Projetos/01_DTI-A/02_Input/01_SP/01_Hidro_hora/Event_1_SPprec.csv",
                     stream_network,
                     warmup_days = 1,
                     warmup_prec = 0
                     ){

  alpha_fraction = param$alpha_fraction
  roughness_global = param$roughness_global
  slope_tolerance = param$slope_tolerance
  boundary_ele = param$boundary_ele
  Raster_WD_Tolerance = param$Raster_WD_Tolerance
  Upstream_Reduction = param$Upstream_Reduction

  Event_name <- paste0("E_",event_name)
  run_name <- paste(Event_name,dem_name,sep='_')
  run_short_name <- paste(Event_name,dem_name,sep="_")

  dem_path <- paste0(Input,path2,dem_name,".tif")

  stream_network_path <- paste0(stream_network,".shp")

  warmup <- data.frame(
    prec = warmup_prec,
    timesteps = seq(0,warmup_days*24*60*60,600)
  )


  rain <- read.csv(rain_path)
  rain$timesteps <- rain$X*10*60+max(warmup$timesteps)

  rain <- rbind(warmup,rain[,c("prec","timesteps")])

  rain$prec <- na_kalman(rain$prec)

  make_it_rain(Input = Input,
               Event_name = Event_name,
               timesteps = rain$timesteps,
               rainfall  = rain$prec*6)


  outlets <- st_read(outlet_path)
  outlets <- st_transform(outlets,crs = 31983)

  outlets <- outlets |>
    subset(ID %in% c(11,157,1000858,1000886))


  make_control_points(Input = Input,
                      Lon = st_coordinates(outlets)[,1],
                      Lat = st_coordinates(outlets)[,2],
                      id = outlets$ID,
                      dem = raster(dem_path,crs = "EPSG:31983"),
                      dem_name = dem_name,
                      Event_name = Event_name,
                      period = 600,
                      stream_network = stream_network_path
  )

  options(scipen = 999)

  caflood_controller(
    Input,
    run_name = run_name,
    run_short_name,
    time_start = 0,
    time_end = max(rain$timesteps),
    max_DT = 60,
    min_DT = 0.01,
    alpha_fraction = alpha_fraction,
    max_iterations = 100000000,
    roughness_global = roughness_global,
    ignore_WD = 0.0001,
    tolerance = 0.0001,
    slope_tolerance = slope_tolerance,
    boundary_ele = boundary_ele,
    dem_path = paste0(dem_name,".asc"),
    dem_name = dem_name,
    rain_path = paste0(Event_name,"_rain.csv"),
    flow_path = NULL,
    output_period = 300,
    Raster_WD_Tolerance = Raster_WD_Tolerance,
    Upstream_Reduction = Upstream_Reduction
  )

  make_it_run(caflood_path,
              Input,
              Output,
              controller = paste0(run_short_name,"_control.csv"))


}



cafloodr2 <- function(param,
                     caflood_path,
                     Input,
                     Output,
                     event_name,
                     dem_name,
                     path2,
                     outlet_path = "C:/Projetos/01_DTI-A/02_Input/01_SP/01_Hidro_hora/estacoes.shp",
                     rain_path = "C:/Projetos/01_DTI-A/02_Input/01_SP/01_Hidro_hora/Event_1_SPprec.csv",
                     stream_network
){

  alpha_fraction = param$alpha_fraction
  roughness_global = param$roughness_global
  slope_tolerance = param$slope_tolerance
  boundary_ele = param$boundary_ele
  Raster_WD_Tolerance = param$Raster_WD_Tolerance
  Upstream_Reduction = param$Upstream_Reduction

  Event_name <- paste0("E_",event_name)
  run_name <- paste(Event_name,dem_name,sep='_')
  run_short_name <- paste(Event_name,dem_name,sep="_")

  dem_path <- paste0(Input,path2,dem_name,".tif")

  stream_network_path <- paste0(stream_network,".shp")

  rain <- fread(rain_path)
  rain$timesteps <- rain$V1*10*60

  prec <- rain$prec

  flow <- hymodr(param = c(499.99941429 ,  1.42245963 ,  0.38878335,   0.02660172 ,  0.25048186),
                 area = 61.28059440489175,
                 tdelta = 60,
                 e = rep(0, length(prec)),
                 p = prec)$q_tot


  outlets <- st_read(outlet_path)
  outlets <- st_transform(outlets,crs = 31983)

  outlets <- outlets |>
    subset(ID %in% c(11,157,1000858,1000886))


  make_control_points(Input = Input,
                      Lon = st_coordinates(outlets)[,1],
                      Lat = st_coordinates(outlets)[,2],
                      id = outlets$ID,
                      dem = raster(dem_path,crs = "EPSG:31983"),
                      dem_name = dem_name,
                      Event_name = Event_name,
                      period = 600,
                      stream_network = stream_network_path
  )

  coord <- fread("C:/Path/Modelos/Aricanduva_lidar/Input/dem_lidar_10m_FM_WLpoints.csv",
                 skip = 2, header = T)

  make_it_flow(Input = Input,
               Event_name = Event_name,
               timesteps = rain$timesteps,
               water_level  = flow,
               tlx = 3$`1000843`[1],
               tly = coord$`1000843`[2])

  options(scipen = 999)

  caflood_controller(
    Input,
    run_name = run_name,
    run_short_name,
    time_start = 0,
    time_end = max(rain$timesteps),
    max_DT = 60,
    min_DT = 0.01,
    alpha_fraction = alpha_fraction,
    max_iterations = 100000000,
    roughness_global = roughness_global,
    ignore_WD = 0.0001,
    tolerance = 0.0001,
    slope_tolerance = slope_tolerance,
    boundary_ele = boundary_ele,
    dem_path = paste0(dem_name,".asc"),
    dem_name = dem_name,
    rain_path = NULL,
    flow_path = paste0(Event_name,"_WaterLevelBC.csv"),
    output_period = 300,
    Raster_WD_Tolerance = Raster_WD_Tolerance,
    Upstream_Reduction = Upstream_Reduction
  )

  make_it_run(caflood_path,
              Input,
              Output,
              controller = paste0(run_short_name,"_control.csv"))


}


