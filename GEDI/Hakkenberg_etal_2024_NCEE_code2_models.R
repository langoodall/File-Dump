# code 2 models ----- Hakkenberg et. al. 2024 NCEE; chrishakkenberg@gmail.com

# Models for Fig. 3 analysis. Univariate GLMM - inter-fire (10min) -----------

  rm(list=ls())
  time1<-Sys.time()
  options(warn=1) 
  
  # input/output
    output_folder <- ""
    spatial_seeds <- read.csv("...6_prefire_sparse_seeds.csv")
    prefire<-fread("...5_prefire.csv", colClasses = list(factor=c("fire","GEDI_beam","region","typeForest")))
    prefire <- as.data.frame(prefire)
    prefire <- prefire[prefire$strct_RH_100>=5,]
    
  # parameters  
    strct_vars <- c("strct_tAGBD","strct_RH_98","strct_nmode","strct_mPAI_sparse_b10") 
    w_t_vars <- c("topo_slope","w_vpd_d_max","w_ET_b4_mean","w_vs_d_max") 
    preds <- c(strct_vars,w_t_vars)
    resps <- "mtbs_dnbrOW" 
    beamtypes <- "all"

    progress.bar <- create_progress_bar("text")
    progress.bar$init(length(preds)*length(resps)*length(beamtypes))    
    for (resp in resps){
      for (pred in preds){
        for (beamtype in beamtypes){
    
        # set up results mat  
          mat_cols <- c("resp","pred","beamtype","n_shots")
          results_mat <- as.data.frame(matrix(NA,1,length(mat_cols)))
          names(results_mat) <- mat_cols
          results_mat[,c("resp","pred","beamtype")] <- c(resp,pred,beamtype)
          
        # filter for current iteration
          tmp <- prefire[,c(resp,pred,"GEDI_beam","fire","GEDI_shot")]
          if (beamtype!="all"){ tmp <- tmp[which(tmp$GEDI_beam==beamtype),]}
          names(tmp) <- c("tmp_resp","tmp_pred","GEDI_beam","fire","GEDI_shot")
          tmp_no0s <- tmp[complete.cases(tmp),]

        # control sample sizes
          results_mat$n_shots <- nrow(tmp_no0s)
          if (nrow(tmp_no0s)<50){next}
          if ( length(unique(tmp_no0s$tmp_pred))<2  ) {
            print(paste("no variance in",pred))
            next
          }    
          
        # check data for normality, transform, and scale
          tmp_sc <- norm_transform_func(tmp_no0s,cols_used=c(1,2))
          tmp_sc[,c("tmp_resp","tmp_pred")] <- scale(tmp_sc[,c("tmp_resp","tmp_pred")])
          
        # remove outliers
          tmp_or <- outlier_removal_cols_withheld(tmp_sc,z=4,cols_witheld=c(3:5))
          
        # cor
          cor_spear <- cor.test(tmp_or$tmp_resp,tmp_or$tmp_pred,method="spearman",exact=FALSE)
          if (is.na(cor_spear$p.value)){next}
          if (cor_spear$p.value<=0.05){ results_mat$cor <- cor_spear$estimate } else { results_mat$cor <- 0 }
        
        # RF R2
          ranger_model <- ranger(tmp_resp~.,data=tmp_or[,c("tmp_resp","tmp_pred","fire")],importance="impurity",verbose=F)
          rf_lm<-lm(tmp_or$tmp_resp ~ ranger_model$predictions)
          results_mat$rfr2 <- summary(rf_lm)$adj.r.squared
    
        # run MEs  
            lin_ME <- ME_wrapper_lin(tmp_resp,tmp_pred,fire,tmp_sc)
            blocked_ME_results <- blocked_quad_ME_func_int(tmp_or,spatial_seeds,seeds=10)

          tmp_range <- as.data.frame(matrix(range(tmp_or$tmp_pred),1,2))
          names(tmp_range) <- c("pred_low","pred_high")   
          
        # tidy
          results_mat2 <- cbind(results_mat,lin_ME,blocked_ME_results,tmp_range)
          results_mat2[unlist(lapply(results_mat2, is.numeric))] <- round(results_mat2[unlist(lapply(results_mat2, is.numeric))],4)
          results_list[[index]] <- results_mat2
          index<-index+1
          progress.bar$step()
            
        } # end beamtype   
      } # end preds
    } # end resps

  results_all <-do.call(rbind.fill,results_list) 
  results_all2 <- results_all[,names(results_all)!="blocked_ME_results"]
  write.csv(results_all2,"...1_uni_inter_dnbrOW_int.csv",row.names=F)
  Sys.time()-time1

  
# Models for Fig. 4 analysis. Bivariate interaction GLMM: strct*env - interfire (4.2) -----------

  rm(list=ls())
  time1<-Sys.time()
  options(warn=1)
  
  # input
    spatial_seeds <- read.csv("...6_prefire_sparse_seeds.csv")
    prefire<-fread("...5_prefire.csv", colClasses = list(factor=c("fire","GEDI_beam","region","typeForest")))
    prefire <- as.data.frame(prefire)
    prefire <- prefire[prefire$strct_RH_100>=5,]

    strct_vars <- c("strct_tAGBD","strct_RH_98","strct_mPAI_sparse_b10","strct_nmode")
    env_vars <- c("topo_slope","w_ET_b4_mean","w_vs_d_max","w_vpd_d_max")

    resps <- "mtbs_dnbrOW"
    strct_var<-"strct_nmode";env_var<-"w_vs_d_max";resp<-"mtbs_dnbrOW";results_list<-list();index<-1
    
  # thread_loop inputs
    file.remove(list.files(tmp_folder)) 
    output_folder <- ""
    output_file_name <- "4_bi_inter.csv"
    threads<-1
    vars_to_split <- strct_vars # change all of these to "tmp_vars"

  # initialize threads
    dir.create(file.path(output_folder, gsub(".csv","_tmp",output_file_name)), showWarnings = FALSE)
    tmp_folder <- paste0(output_folder, gsub(".csv","_tmp",output_file_name),"/")
    setwd(tmp_folder)
    if ( file.exists("final_final.csv") ) { file.remove(list.files(tmp_folder)) }
    step <- round(length(vars_to_split)/threads) ; thread <- 1
    for (thread in 1:(threads+1)){
       thread_initial <- (length(list.files(path=tmp_folder,pattern="initial"))+1)
       
       if (thread_initial==threads) {
         tmp_vars <- vars_to_split[((step*thread_initial)-step+1):length(vars_to_split)]
       } else if (thread_initial<threads) { tmp_vars <- vars_to_split[((step*thread_initial)-step+1):(step*thread_initial)] }
       if (thread_initial<=threads) { write.csv(as.data.frame(matrix(tmp_vars,1,length(tmp_vars))),paste0(tmp_folder,"initial",thread_initial,".csv"),row.names=F)}
       if (thread_initial<=threads) {
         
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      
    progress.bar <- create_progress_bar("text")
    progress.bar$init(length(resps)*length(tmp_vars)*length(env_vars))    
    for (resp in resps){
      for (strct_var in tmp_vars){
        for (env_var in env_vars){
        
          # set up results mat  
            mat_cols <- c("resp_var","strct_var","env_var")
            results_mat <- as.data.frame(matrix(NA,1,length(mat_cols)))
            names(results_mat) <- mat_cols
            results_mat[,c("resp_var","strct_var","env_var")] <- c(resp,strct_var,env_var)

          # filter for current iteration
            tmp <- prefire[,c(resp,strct_var,env_var,"fire","GEDI_shot")]
            names(tmp)[1:3] <- c("resp_var","strct_var","env_var")
            tmp <- tmp[complete.cases(tmp),]
            tmp_no0s <- tmp#[which(tmp$strct_var!=0),] # remove 0s  @@@@@@@@@
            
          # control sample sizes
            results_mat$n_shots <- nrow(tmp_no0s)
            if (nrow(tmp_no0s)<50){next}
            if ( length(unique(tmp_no0s$env_var))<2  ) {
              print(paste("no variance in",env_var))
              next
            }    
            
          # check data for normality, transform, and scale
            tmp_no0s <- norm_transform_func(tmp_no0s,cols_used=c(1:3))
            tmp_no0s[,c("resp_var","strct_var","env_var")] <- scale(tmp_no0s[,c("resp_var","strct_var","env_var")])

          # remove outliers
            tmp_or <- outlier_removal_cols_withheld(tmp_no0s,z=4,cols_witheld=c(4:5))
            
          # run blocked MEs
            blocked_ME_results_RSI <- blocked_bivar_lin_RSRII(tmp_or,spatial_seeds,seeds=10)
            names(blocked_ME_results_RSI) <- paste0("RII_",names(blocked_ME_results_RSI))
            
          # tidy
            results_mat2 <- cbind(results_mat,blocked_ME_results_RSI)
            results_mat2[unlist(lapply(results_mat2, is.numeric))] <- round(results_mat2[unlist(lapply(results_mat2, is.numeric))],4)
            results_list[[index]] <- results_mat2
            index<-index+1
            progress.bar$step()
            
        } # end env_var
      } # end strct_var
    } # end resp_var
    results_all <-do.call(rbind.fill,results_list)
         
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         
   write.csv(results_all,paste0(tmp_folder,"final",thread_initial,".csv"),row.names=F)
   thread_final <- (length(list.files(path=tmp_folder,pattern="final"))+1) 
   
   }  else if (thread_final<=threads) {
     next
   } else if (thread_final>threads) {
     if (length(list.files(path=tmp_folder,pattern="final_final"))==1) {break}
     files_final <- list.files(path=tmp_folder,pattern="final")
     final_list<-list();index2<-1;file<-1
     for (file in 1:length(files_final)){
       final_list[[index2]] <- read.csv(files_final[file])
       index2 <- index2+1
     }
     
     results_final_final <- do.call(rbind.fill,final_list)
     write.csv(results_final_final,paste0(tmp_folder,"final_final.csv"),row.names=F)
     write.csv(results_final_final,paste0(output_folder,output_file_name),row.names=F)
     
    }
  }
  Sys.time()-time1  
   
  
# Models Fig. 5 analysis. Subsets: Univariate GLMM - topo and weather (21min) -----------

  rm(list=ls())
  time1<-Sys.time()

  options(warn=1)

  # input
    spatial_seeds <- read.csv("...6_prefire_sparse_seeds.csv")
    prefire<-fread("...5_prefire.csv", colClasses = list(factor=c("fire","GEDI_beam","region","typeForest")))
    prefire <- as.data.frame(prefire)
    prefire <- prefire[prefire$strct_RH_100>=5,]

  # parameters  
    fires <- unique(prefire$fire)
    strct_vars <- c("strct_tAGBD","strct_RH_98","strct_mPAI_sparse_b10","strct_nmode")
    subsets <- c("topo_slope","w_ET_b4_mean","w_vs_d_max","w_vpd_d_max")

    progress.bar <- create_progress_bar("text")
    progress.bar$init(length(subsets)*length(resps)*length(strct_vars))    
    for (resp in resps){
      for (strct_var in strct_vars){
        for (subset in subsets){
            
          # set up results mat  
            mat_cols <- c("tmp_resp","tmp_pred","subset","quant","n_shots")
            results_mat <- as.data.frame(matrix(NA,1,length(mat_cols)))
            names(results_mat) <- mat_cols
            results_mat[1,c("tmp_resp","tmp_pred","subset")] <- c(resp,strct_var,subset)
         
          # filter for current iteration
            tmp <- prefire[,c(resp,strct_var,subset,"fire","GEDI_shot")] 
            tmp <- tmp[complete.cases(tmp),]
            names(tmp)[1:3] <- c("tmp_resp","tmp_pred","tmp_subset")
            tmp_no0s <- tmp
            
          # remove outliers
            tmp_or <- outlier_removal_cols_withheld(tmp_no0s,z=3,cols_witheld=c(4:5))
            
          # control sample sizes
            results_mat$n_shots <- nrow(tmp_or)
            if (nrow(tmp_or)<50){next}
            if ( (length(unique(tmp_or$tmp_resp))<2) | (length(unique(tmp_or$tmp_pred))<2)  ) {next}
             
          # check data for normality, transform, and scale
            tmp_sc <- norm_transform_func(tmp_or,cols_used=c(1:2))
            tmp_sc[,c("tmp_resp","tmp_pred")] <- scale(tmp_sc[,c("tmp_resp","tmp_pred")])

            tmp_quants <- quantile(tmp_sc$tmp_subset,probs=seq(0,1,.33))
          
            for (quant in c("all","low","mid","high")){
              if (quant=="all"){ tmp_subset <- tmp_sc  }
              if (quant=="low"){ tmp_subset <- tmp_sc[which(tmp_sc$tmp_subset<tmp_quants[2]),] }
              if (quant=="mid"){ tmp_subset <- tmp_sc[which(tmp_sc$tmp_subset>=tmp_quants[2] & tmp_sc$tmp_subset<=tmp_quants[3] ),] }
              if (quant=="high"){ tmp_subset <- tmp_sc[which(tmp_sc$tmp_subset>tmp_quants[3]),] }
              results_mat$quant <- quant
              
              # run blocked MEs
                blocked_ME_results <- blocked_linquad_ME_func(tmp_subset,spatial_seeds,seeds=10)
                
              # tidy
                results_mat2 <- cbind(results_mat,blocked_ME_results)
                results_mat2[unlist(lapply(results_mat2, is.numeric))] <- round(results_mat2[unlist(lapply(results_mat2, is.numeric))],4)
                results_list[[index]] <- results_mat2
                index<-index+1
            }
  
        progress.bar$step()
      } # end subset
    } # end strct_var
  } # end resps
  
  results_all <-do.call(rbind.fill,results_list)
  results_all2 <- results_all[,names(results_all)!="blocked_ME_results"]
  write.csv(results_all2,"...3_uni_inter_subsets.csv",row.names=F)
  Sys.time()-time1
  
# Models for Fig 6. analysis Univariate GLMM - intra-fire (8.5min - 19 preds 2 resps) -----------

  rm(list=ls())
  time1<-Sys.time()

  # input
    spatial_seeds <- read.csv("...6_prefire_sparse_seeds.csv")
    prefire<-fread("...5_prefire.csv", colClasses = list(factor=c("fire","GEDI_beam","region","typeForest")))
    prefire <- as.data.frame(prefire)
    prefire <- prefire[prefire$strct_RH_100>=5,]

  # parameters
    fires <- unique(prefire$fire)
    strct_vars <- c("strct_tAGBD","strct_RH_98","strct_nmode","strct_mPAI_sparse_b10")
    preds <- strct_vars 
    resps <- "mtbs_dnbrOW" 
    results_list<-list();index<-1;pred<-"strct_VDR";resp<-"mtbs_dnbrOW";fire<-"CARMEL"

  for (fire in fires){
    print(fire)
    progress.bar <- create_progress_bar("text")
    progress.bar$init(length(preds)*length(resps))

    for (resp in resps){
      for (pred in preds){

        # set up results mat
          mat_cols <- c("fire","resp","pred","n_shots")
          results_mat <- as.data.frame(matrix(NA,1,length(mat_cols)))
          names(results_mat) <- mat_cols
          results_mat[,c("fire","resp","pred")] <- c(fire,resp,pred)

        # filter for current iteration
          tmp <- prefire[prefire$fire==fire,c(resp,pred,"GEDI_shot")] 
          tmp <- tmp[complete.cases(tmp),]
          names(tmp)[1:2] <- c("tmp_resp","tmp_pred")
          tmp_no0s <- tmp

        # remove outliers
          tmp_or <- outlier_removal_cols_withheld(tmp_no0s,z=3,cols_witheld=c(3))

        # control sample sizes
          results_mat$n_shots <- nrow(tmp_or)
          if (nrow(tmp_or)<20){next}
          if ( (length(unique(tmp_or$tmp_resp))<2) | (length(unique(tmp_or$tmp_pred))<2)  ) {next}

        # check data for normality, transform, and scale
          tmp_sc <- norm_transform_func(tmp_or,cols_used=c(1,2))
          tmp_sc[,c("tmp_resp","tmp_pred")] <- scale(tmp_sc[,c("tmp_resp","tmp_pred")])

        # cor
          cor_spear <- cor.test(tmp_sc$tmp_resp,tmp_sc$tmp_pred,method="spearman",exact=FALSE)
          if (is.na(cor_spear$p.value)){next}
          if (cor_spear$p.value<=0.05){ results_mat$cor <- cor_spear$estimate } else { results_mat$cor <- 0 }

        # RF R2
          ranger_model <- ranger(tmp_resp~tmp_pred,data=tmp_sc[,c("tmp_resp","tmp_pred")],importance="impurity")
          rf_lm<-lm(tmp_sc$tmp_resp ~ ranger_model$predictions)
          results_mat$rfr2 <- summary(rf_lm)$adj.r.squared

        # run blocked linquad (not MEs)
          blocked_ME_results <- blocked_lin_func(tmp_sc,spatial_seeds,seeds=10)

        # tidy
          results_mat2 <- cbind(results_mat,blocked_ME_results)
          results_mat2[unlist(lapply(results_mat2, is.numeric))] <- round(results_mat2[unlist(lapply(results_mat2, is.numeric))],4)
          results_list[[index]] <- results_mat2
          index<-index+1
          progress.bar$step()

      } # end pred
    }  # end resp
  } #end fire_csv

  results_all <-do.call(rbind.fill,results_list)
  results_all2 <- results_all[,names(results_all)!="blocked_ME_results"]
  write.csv(results_all2,"...2_uni_intra.csv",row.names=F)
  Sys.time()-time1
