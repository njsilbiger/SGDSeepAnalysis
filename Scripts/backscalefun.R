

backscale<- function(fit, x.name, y.name, season, site, Modeldata.x,Modeldata.y) {
  R1<-conditional_effects(fit, x.name, resp = y.name, resolution = 100)[[1]] %>%
      mutate(estimate = estimate__*attr(Modeldata.y,"scaled:scale")+attr(Modeldata.y,"scaled:center"),
           lower = lower__*attr(Modeldata.y,"scaled:scale")+attr(Modeldata.y,"scaled:center"),
           upper = upper__*attr(Modeldata.y,"scaled:scale")+attr(Modeldata.y,"scaled:center"),
           season = season, 
           site = site,
  
    x = R1[,x.name]*attr(Modeldata.x,"scaled:scale")+attr(Modeldata.x,"scaled:center"))
  return(R1)
}

b1<-backscale(fit = V_Dry_fit, x.name = "SilicateumolL", y.name = "Proteinaceous", season = "Dry", site = "Varari",Modeldata.x =ModelData$SilicateumolL,Modeldata.y = ModelData$Proteinaceous)

b2<-backscale(fit = V_Wet_fit, x.name = "SilicateumolL", y.name = "Proteinaceous", season = "Wet", site = "Varari",Modeldata.x =ModelData$SilicateumolL,Modeldata.y = ModelData$Proteinaceous)


