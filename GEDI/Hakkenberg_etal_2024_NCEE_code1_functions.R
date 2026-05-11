# code 1 functions ------Hakkenberg et. al. 2024 NCEE; chrishakkenberg@gmail.com

 # mean_morans_p_func --------- 
  mean_morans_p_func <- function(blocking_output,input) {
    morans_p <- vector()
    for (seed in 1:10) {
      tmp <- blocking_output[blocking_output$seed==seed,]
      tmp2 <- input[input$GEDI_shot %in% tmp$GEDI_shot,]
      tmp_dist_mat <- as.matrix(dist(cbind(tmp2$GEDI_utmX, tmp2$GEDI_utmY)))
      tmp_dist_inv <- 1/tmp_dist_mat
      diag(tmp_dist_inv) <- 0
      my_moran <- Moran.I(tmp2$strct_RH_98, tmp_dist_inv) # sig p value means there is autocorrelation
      morans_p <- c(morans_p,my_moran$p.value)
    }
    return(mean(morans_p))
  }

# blocking_across_seeds_func -------

  
  blocking_across_seeds_func <- function (input,my_range,seeds=1:10) {
    footprint_dist_mat <- as.matrix(dist(cbind(input$GEDI_utmX, input$GEDI_utmY)))
    my_list<-list()
    for (seed in seeds){   # start seeds
      input$proximity_index<-1
      set.seed(seed)
      
      for (i in sample(nrow(input))){  # start dist_mat
        if (input[i,"proximity_index"]==1) {
          tmp_index<-footprint_dist_mat[i,][which(footprint_dist_mat[i,] < my_range & footprint_dist_mat[i,]!=0)] 
          if (length(tmp_index)==0){next} else
          {
            input[as.numeric(names(tmp_index)),"proximity_index"] <- 0
          }
        } else {next}
      } # end dist_mat
      
      tmp_output <- as.data.frame(cbind(fire,seed,input[input$proximity_index==1,"GEDI_shot"]))
      tmp_output$seed <- as.integer(tmp_output$seed)
      names(tmp_output) [3] <- "GEDI_shot"
      
      my_list[[seed]] <- tmp_output
    } # end seeds
    output <- do.call(rbind,my_list)
    return(output)
  }
  

my_theme <- function(base_size = 12, base_family = "",mycolor_background="black",plot_background="white",mycolor_elements="white",mycolor_text="white",font="Helvetica") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_blank(),  
      axis.text.x = element_text(family=font,size = base_size*0.6, color = mycolor_text, lineheight = 0.9,angle=25),  
      axis.text.y = element_text(family=font,size = base_size*0.6, color = mycolor_text, lineheight = 0.9),  
      axis.ticks = element_line(color = mycolor_elements, linewidth  =  0.2),  
      axis.title.x = element_text(family=font,size = base_size, color = mycolor_text, margin = ggplot2::margin(0, 10, 0, 0)),  
      axis.title.y = element_text(family=font,size = base_size, color = mycolor_text, angle = 90, margin = ggplot2::margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      legend.background = element_rect(color = plot_background, fill = plot_background),  
      legend.key = element_rect(color = plot_background,  fill = plot_background),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(family=font,size = base_size*0.8, color = mycolor_text),  
      legend.title = element_text(family=font,size = base_size*0.8, face = "bold", hjust = 0, color = mycolor_text),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      panel.background = element_rect(fill = mycolor_background, color  =  NA),  
      panel.border = element_rect(fill = NA, color = mycolor_elements),  
      panel.grid.major = element_line(color = mycolor_background),  
      panel.grid.minor = element_blank(),  
      panel.spacing = unit(0.5, "lines"),   
      plot.background = element_rect(color = plot_background, fill = plot_background),  
      plot.title = element_text(size = base_size*1.2, color = mycolor_text,face="bold",family=font),  
      plot.margin = unit(rep(1, 4), "lines")
    )
}

# ME_wrapper_lin ------

ME_wrapper_lin <- function(tmp_resp,tmp_pred,fire,data){
  skip_to_next  <- FALSE
  tryCatch(re_lme <- lme(tmp_resp ~ tmp_pred, random=~tmp_pred|fire,data=data,method="ML",control =list(msMaxIter = 1000, msMaxEval = 1000)), error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next) {
    re_lme <- lmer(tmp_resp ~ tmp_pred + (tmp_pred | fire), data)
    results <- as.data.frame(summary(re_lme)$coefficients)
  } else {
    results <- as.data.frame(summary(re_lme)$tTable)
  }
    names(results)[1:2]<- c("mCI","se")
    results$lCI<-results[,"mCI"]-results[,"se"]*1.96
    results$hCI<-results[,"mCI"]+results[,"se"]*1.96
    results2 <- results[2,c("lCI","mCI","hCI")]
    results2$sig <- "nonsig"
    row.names(results2) <- NULL
    if ( (results2[,"lCI"]<0 & results2[,"hCI"]<0) | (results2[,"lCI"]>0 & results2[,"hCI"]>0) ) {results2$sig <- "sig"}
  results3 <- cbind(results2[1,],r.squaredGLMM(re_lme),cor(data$tmp_resp,predict(re_lme))^2)
  names(results3) <- c("lCI","mCI","hCI","sig","R2m","R2c","R2alt")
  return(results3)
}

# blocked_lin_func ------

blocked_lin_func <- function(input_mat,spatial_seeds_mat,seeds=10){
  seed<-1;seed_list <- list()
  for (seed in 1:seeds){
    
    tmp_seed <- input_mat[input_mat$GEDI_shot %in% spatial_seeds_mat[spatial_seeds_mat$seed==seed,"GEDI_shot"],]
    if (nrow(tmp_seed)<10){next} 
    
    # lm quad
    my_lm<-lm(tmp_resp ~ tmp_pred, data = tmp_seed)
    results <- as.data.frame(summary(my_lm)$coefficients)
    names(results)[1:2]<- c("mCI","se")
    results$lCI<-results[,"mCI"]-results[,"se"]*1.96
    results$hCI<-results[,"mCI"]+results[,"se"]*1.96
    results2 <- results[2,c("lCI","mCI","hCI")]
    results3 <- cbind(results2[1,],summary(my_lm)$adj.r.squared)
    names(results3) <- c("lCI","mCI","hCI","adjR2")
    row.names(results3) <- NULL
    seed_list[[seed]] <- cbind(as.matrix(nrow(tmp_seed)),results3)
    
  }
  
  seeds_mat <-do.call(rbind,seed_list) 
  
  if (is.null(seeds_mat)) {seeds_mean3 <- NA
  } else {
    
    names(seeds_mat)[1] <- "nshot_seed_mean"
    seeds_mean <- apply(seeds_mat[,!names(seeds_mat) %ilike% "sig"],2,mean)
    seeds_mean2 <- as.data.frame(matrix(seeds_mean,1,length(seeds_mean)))
    names(seeds_mean2) <- names(seeds_mean)
    seeds_mean2$mCI_sig <-  "nonsig"
    if ( (seeds_mean2[1,"lCI"]<0 & seeds_mean2[1,"hCI"]<0) | (seeds_mean2[1,"lCI"]>0 & seeds_mean2[1,"hCI"]>0) ) {seeds_mean2$mCI_sig <- "sig"}
    seeds_mean2$nshot_seed_mean <- round(seeds_mean2$nshot_seed_mean)
  }
  return(seeds_mean2)
}

# blocked_bivar_lin_RSRII --------

blocked_bivar_lin_RSRII <- function(input_mat,spatial_seeds_mat,seeds=10){
  seed<-1;seed_list <- list()
  
  for (seed in 1:seeds){
    
    tmp_seed <- input_mat[input_mat$GEDI_shot %in% spatial_seeds_mat[spatial_seeds_mat$seed==seed,"GEDI_shot"],]
    tmp_seed <- tmp_seed[tmp_seed$fire%in%names(which(table(tmp_seed$fire)>500)),]
    if (nrow(tmp_seed)<10){next} 
    
    skip_to_next  <- FALSE
    tryCatch(re_lme <- lme(resp_var ~ strct_var*env_var, random=~strct_var*env_var|fire,data=tmp_seed,method="ML",control =list(msMaxIter = 1000, msMaxEval = 1000)), error = function(e) {skip_to_next <<- TRUE})
    if(skip_to_next) {
      re_lme <- lmer(resp_var ~ strct_var*env_var + (strct_var*env_var | fire), tmp_seed)
      results <- as.data.frame(summary(re_lme)$coefficients)
    } else {
      results <- as.data.frame(summary(re_lme)$tTable)
    }
    names(results)[1:2]<- c("mCI","se")
    results$lCI<-results[,"mCI"]-results[,"se"]*1.96
    results$hCI<-results[,"mCI"]+results[,"se"]*1.96
    results2 <- results[2:4,c("lCI","mCI","hCI")]
    results3 <- cbind(results2[1,],results2[2,],results2[3,],r.squaredGLMM(re_lme),cor(tmp_seed$resp_var,predict(re_lme))^2)
    names(results3) <- c("lCI_strct","mCI_strct","hCI_strct","lCI_env","mCI_env","hCI_env","lCI_int","mCI_int","hCI_int","R2m","R2c","R2alt")
    row.names(results3) <- NULL
    seed_list[[seed]] <- cbind(as.matrix(nrow(tmp_seed)),results3)
    
  }
  
  seeds_mat <-do.call(rbind,seed_list) 
  
  if (is.null(seeds_mat)) {seeds_mean3 <- NA
  } else {
    
    names(seeds_mat)[1] <- "nshot_seed_mean"
    seeds_mean <- apply(seeds_mat[,!names(seeds_mat) %ilike% "sig"],2,mean)
    seeds_mean2 <- as.data.frame(matrix(seeds_mean,1,length(seeds_mean)))
    names(seeds_mean2) <- names(seeds_mean)
    seeds_cv <- apply(seeds_mat[,names(seeds_mat) %ilike% "mCI"],2,function(x){sd(x)/mean(x)})
    seeds_cv2 <- as.data.frame(matrix(seeds_cv,1,length(seeds_cv)))
    names(seeds_cv2) <- paste0(names(seeds_cv),"_seed_cv")
    seeds_mean2$strct_sig <- seeds_mean2$env_sig <-  seeds_mean2$int_sig <- "nonsig"
    if ( (seeds_mean2[1,"lCI_strct"]<0 & seeds_mean2[1,"hCI_strct"]<0) | (seeds_mean2[1,"lCI_strct"]>0 & seeds_mean2[1,"hCI_strct"]>0) ) {seeds_mean2$strct_sig <- "sig"}
    if ( (seeds_mean2[1,"lCI_env"]<0 & seeds_mean2[1,"hCI_env"]<0) | (seeds_mean2[1,"lCI_env"]>0 & seeds_mean2[1,"hCI_env"]>0) ) {seeds_mean2$env_sig <- "sig"}
    if ( (seeds_mean2[1,"lCI_int"]<0 & seeds_mean2[1,"hCI_int"]<0) | (seeds_mean2[1,"lCI_int"]>0 & seeds_mean2[1,"hCI_int"]>0) ) {seeds_mean2$int_sig <- "sig"}
    
    #names(seeds_mean2) <- paste0(names(seeds_mean2),"_seed_mean")
    seeds_mean2$nshot_seed_mean <- round(seeds_mean2$nshot_seed_mean)
    seeds_mean3 <- cbind(seeds_mean2,seeds_cv2)
  }
  return(seeds_mean3)
}

# blocked_quad_ME_func_int ---------

blocked_quad_ME_func_int <- function(input_mat,spatial_seeds_mat,seeds=10){

   # seeds<-10
  # input_mat <- tmp_sc #c(tmp_resp,tmp_pred,fire,GEDI_shot)
  # spatial_seeds_mat <- spatial_seeds #c(fire,seed,GEDI_shot)
  
  seed<-1;seed_list <- list()
  for (seed in 1:seeds){
    
    tmp_seed <- input_mat[input_mat$GEDI_shot %in% spatial_seeds_mat[spatial_seeds_mat$seed==seed,"GEDI_shot"],]
        tmp_seed <- tmp_seed[tmp_seed$fire%in%names(which(table(tmp_seed$fire)>500)),]
        if (nrow(tmp_seed)<50){next} # for multi-fire ME
    
    # run MEs  

    tmp_lin_quad_ME <- ME_wrapper_lin_quad_int(tmp_resp,tmp_pred,fire,tmp_seed)
    seed_list[[seed]] <- cbind(as.matrix(nrow(tmp_seed)),tmp_lin_quad_ME)
    
  }
  
  seeds_mat <-do.call(rbind,seed_list) 
  
  if (is.null(seeds_mat)) {seeds_mean3 <- NA
  } else {
    
    names(seeds_mat)[1] <- "nshot_seed_mean"
    seeds_mean <- apply(seeds_mat[,!names(seeds_mat) %ilike% "sig"],2,mean)
    seeds_mean2 <- as.data.frame(matrix(seeds_mean,1,length(seeds_mean)))
    names(seeds_mean2) <- names(seeds_mean)
    seeds_cv <- apply(seeds_mat[,names(seeds_mat) %ilike% "mCI"],2,function(x){sd(x)/mean(x)})
    seeds_cv2 <- as.data.frame(matrix(seeds_cv,1,length(seeds_cv)))
    names(seeds_cv2) <- paste0(names(seeds_cv),"_seed_cv")
    
    seeds_mean2$int_mCI_sig <- seeds_mean2$quad_mCI_sig <- seeds_mean2$lin_mCI_sig <-  "nonsig"
    if ( (seeds_mean2[1,"int_lCI"]<0 & seeds_mean2[1,"int_hCI"]<0) | (seeds_mean2[1,"int_lCI"]>0 & seeds_mean2[1,"int_hCI"]>0) ) {seeds_mean2$int_mCI_sig <- "sig"}
    if ( (seeds_mean2[1,"lin_lCI"]<0 & seeds_mean2[1,"lin_hCI"]<0) | (seeds_mean2[1,"lin_lCI"]>0 & seeds_mean2[1,"lin_hCI"]>0) ) {seeds_mean2$lin_mCI_sig <- "sig"}
    if (all(!is.na(seeds_mat$quad_mCI))) {
      if ( (seeds_mean2[1,"quad_lCI"]<0 & seeds_mean2[1,"quad_hCI"]<0) | (seeds_mean2[1,"quad_lCI"]>0 & seeds_mean2[1,"quad_hCI"]>0) ) {seeds_mean2$quad_mCI_sig <- "sig"}
    } else {
      seeds_mean2$quad_mCI_sig <- NA
    }
    seeds_mean2$nshot_seed_mean <- round(seeds_mean2$nshot_seed_mean)
    seeds_mean3 <- cbind(seeds_mean2,seeds_cv2)
  }
  return(seeds_mean3)
}

# blocked_linquad_ME_func ------

blocked_linquad_ME_func <- function(input_mat,spatial_seeds_mat,seeds=10){
  seed<-1;seed_list <- list()
  for (seed in 1:seeds){
    
    tmp_seed <- input_mat[input_mat$GEDI_shot %in% spatial_seeds_mat[spatial_seeds_mat$seed==seed,"GEDI_shot"],]
    tmp_seed <- tmp_seed[tmp_seed$fire%in%names(which(table(tmp_seed$fire)>500)),]
    if (nrow(tmp_seed)<50){next} # for multi-fire ME
    
    # run MEs  
    tmp_lin_ME <- ME_wrapper_lin(tmp_resp,tmp_pred,fire,tmp_seed)
    tmp_lin_quad_ME <- ME_wrapper_lin_quad(tmp_resp,tmp_pred,fire,tmp_seed)
    
    seed_list[[seed]] <- cbind(as.matrix(nrow(tmp_seed)),tmp_lin_ME,tmp_lin_quad_ME)
    
  }
  
  seeds_mat <-do.call(rbind,seed_list) 
  
  if (is.null(seeds_mat)) {seeds_mean3 <- NA
  } else {
    
names(seeds_mat)[1] <- "nshot_seed_mean"
    seeds_mean <- apply(seeds_mat[,!names(seeds_mat) %ilike% "sig"],2,mean)
    seeds_mean2 <- as.data.frame(matrix(seeds_mean,1,length(seeds_mean)))
    names(seeds_mean2) <- names(seeds_mean)
    seeds_cv <- apply(seeds_mat[,names(seeds_mat) %ilike% "mCI"],2,function(x){sd(x)/mean(x)})
    seeds_cv2 <- as.data.frame(matrix(seeds_cv,1,length(seeds_cv)))
    names(seeds_cv2) <- paste0(names(seeds_cv),"_seed_cv")
    seeds_mean2$quad_mCI_sig <- seeds_mean2$lin_mCI_sig <- seeds_mean2$mCI_sig <-  "nonsig"
    if ( (seeds_mean2[1,"lCI"]<0 & seeds_mean2[1,"hCI"]<0) | (seeds_mean2[1,"lCI"]>0 & seeds_mean2[1,"hCI"]>0) ) {seeds_mean2$mCI_sig <- "sig"}
    if ( (seeds_mean2[1,"lin_lCI"]<0 & seeds_mean2[1,"lin_hCI"]<0) | (seeds_mean2[1,"lin_lCI"]>0 & seeds_mean2[1,"lin_hCI"]>0) ) {seeds_mean2$lin_mCI_sig <- "sig"}
    if (all(!is.na(seeds_mat$quad_mCI))) {
      if ( (seeds_mean2[1,"quad_lCI"]<0 & seeds_mean2[1,"quad_hCI"]<0) | (seeds_mean2[1,"quad_lCI"]>0 & seeds_mean2[1,"quad_hCI"]>0) ) {seeds_mean2$quad_mCI_sig <- "sig"}
    } else {
      seeds_mean2$quad_mCI_sig <- NA
    }
    seeds_mean2$nshot_seed_mean <- round(seeds_mean2$nshot_seed_mean)
    seeds_mean3 <- cbind(seeds_mean2,seeds_cv2)
  }
  return(seeds_mean3)
}

# norm_transform_func ------

norm_transform_func <- function(df,cols_used=c(1,2)){
  for (tmp_col in cols_used){
    tmp_sample<-sample(df[,tmp_col],min(length(df[,tmp_col]),5000))
    tmp_lambda <- transformTukey(tmp_sample,returnLambda=T,quiet=T,plotit=F)
    if (tmp_lambda >  0){df[,tmp_col] = df[,tmp_col] ^ tmp_lambda} 
    if (tmp_lambda == 0){df[,tmp_col] = log(df[,tmp_col])} 
    if (tmp_lambda <  0){df[,tmp_col] = -1 * df[,tmp_col] ^ tmp_lambda} 
  }
   df2 <- df[complete.cases(df),]
  return(df2)
}

# outlier_removal_cols_withheld ------

outlier_removal_cols_withheld <- function(df,z=4,cols_witheld=c(3,4)){
      tmp_2 <- df[,names(df)[-c(cols_witheld)] ]
      tmp_3 <- as.matrix(apply(as.matrix(tmp_2),2,scale))
      if (ncol(tmp_3)==1) {tmp_2[which(abs(tmp_3)>z)] <- NA}
      if (ncol(tmp_3)>1) {tmp_2[which(abs(tmp_3)>z,arr.ind=T)] <- NA}
      df[,names(df)[-c(cols_witheld)] ] <- tmp_2
      df2 <- df[complete.cases(df),]
      return(df2)
    }
    
# mycor2m_sig -------
    
mycor2m_sig<-function (p,r,method){
  p <- as.data.frame(p)
  cor_mat<-data.frame(matrix(nrow=dim(p)[2], ncol=dim(r)[2]))
  colnames(cor_mat)<-colnames(r)      
  row.names(cor_mat)<-colnames(p)  
  for (i in 1:dim(p)[2]){
    for (j in 1:dim(r)[2]){
      tmp_mat<-cbind(p[,i],r[,j])  
      tmp_mat<-tmp_mat[complete.cases(tmp_mat),]
      if (nrow(tmp_mat)==0){cor<-NA
      } else{
        cor<-cor.test(tmp_mat[,1],tmp_mat[,2],method=method,exact=FALSE)  
        if(!is.na(cor$p.value) & cor$p.value<0.05){cor_mat[i,j]<-cor$estimate
        } else {
          cor_mat[i,j]<-0
        }
      }
    }
  }
  return(cor_mat)
}
