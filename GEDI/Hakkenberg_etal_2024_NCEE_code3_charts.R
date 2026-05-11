# code 3 charts ----- Hakkenberg et. al. 2024 NCEE; chrishakkenberg@gmail.com

# Fig. 2 - strct profile distributions -----
    
    rm(list=ls())
    source("")
    options(warn=1)
    output_folder <- ""
    
    prefire<-fread("", colClasses = list(factor=c("fire","GEDI_beam","region","typeForest")))
    prefire <- as.data.frame(prefire)
    hist(prefire$strct_CBH)
    
    unique(prefire$fire)
    strct_vars <- c("strct_tAGBD","strct_RH_98","strct_mPAI_sparse_b10","strct_nmode")
    size_tmp<-9

    for (strct_var in strct_vars){
      tmp <- prefire[,c(strct_var,names(prefire) [names(prefire)%ilike% "_inc"])]
      tmp_quant <- quantile(tmp[,strct_var],probs = seq(0, 1, 0.33),na.rm=T)
      tmp[ which(tmp[,strct_var]<tmp_quant[2]),"Quantile"] <- "low"
      tmp[ which(tmp[,strct_var]>=tmp_quant[2] & tmp[,strct_var]<=tmp_quant[3]),"Quantile"] <- "mid"
      tmp[ which(tmp[,strct_var]>tmp_quant[3]),"Quantile"] <- "high"
      tmp$Quantile <- as.factor(tmp$Quantile)
      tmp$Quantile <- factor(tmp$Quantile, levels = rev(c("low","mid","high")))
      tmp2 <- tmp[which(!is.na(tmp$Quantile)),]
      tmp_mean <- as.data.frame(tmp2 %>% dplyr::group_by(Quantile) %>% dplyr::summarise_all(mean))
      tmp3 <- gather(tmp_mean, "height_class", "PAI_mean", strct_PAI_inc_00_05:strct_PAI_inc_100_105, factor_key=TRUE)
      tmp3$height <- rep(seq(5,105,5),each=3)
      tmp4 <- tmp3[tmp3$height<50,]
      tmp4$height <- as.factor( tmp4$height)
      tmp4$metric <- strct_var
      my_list[[index]] <- tmp4[,-2]
      index <- index +1
    }
    
    all_metrics <- do.call(rbind,my_list)
    all_metrics$metric <- as.factor(all_metrics$metric)
    all_metrics$metric <- factor(all_metrics$metric, levels = c("strct_tAGBD","strct_RH_98","strct_nmode","strct_mPAI_sparse_b10"))

    my_labeller = as_labeller(c(
      strct_tAGBD = "Biomass~(AGBD)", 
      strct_RH_98 = "Canopy~height~(RH98)", 
      strct_nmode = "Layering~(nmode)", 
      strct_mPAI_sparse_b10 = "Ladder~fuels~(mPAI['10m'])"), default = label_parsed)

    fig2_fill <- c("#440154FF" ,"#7FC97F" ,"#FDE725FF")
    fig2_colour <- c("#FDE725FF","black","#440154FF")
    
    dist_plot <- ggplot(data=all_metrics, aes(x=height, y=PAI_mean, colour=Quantile,fill=Quantile, alpha=Quantile)) +
                    geom_bar(stat="identity", position = position_dodge(width = -0.3), width = 3,size=.1)   + 
                    scale_colour_manual(values=fig2_colour) +
                    scale_fill_manual(values=fig2_fill) + 
                    facet_grid(cols = vars(metric),scales = "fixed",labeller = my_labeller) + 
                    scale_alpha_manual(values=c(1,1,1)) +
                    xlab("Height (m)") + ylab("PAI (mean)") +
                    coord_flip() + 
                    my_theme(mycolor_background="grey95",mycolor_text="black") +
                    theme(legend.position=c(.99,0.98),
                      legend.key.height = unit(1,"line"),
                      legend.key.size = unit(1,"line"),
                      legend.margin = ggplot2::margin(2,4,2,2),
                      legend.background = element_rect(colour = 'black',fill=alpha('grey98', 0.8),linewidth = .2),
                      legend.justification=c(1,1),
                      axis.text.x = element_text(angle = 45, vjust = 1, hjust=.9,size=size_tmp-2),
                      axis.text.y = element_text(vjust = -.5, hjust=1,size=size_tmp-2),
                      axis.title.x=element_text(family="Helvetica",size=size_tmp),
                      axis.title.y=element_text(family="Helvetica",size=size_tmp),
                      plot.title = element_text(size=size_tmp+1,vjust=1),
                      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
                      strip.text = element_text(size = size_tmp-1),
                      legend.text.align = 0,
                      legend.text=element_text(size=size_tmp-1),
                      legend.title=element_text(size=size_tmp),
                      legend.key = element_rect(linewidth = 1),
                      legend.spacing.x = unit(.1,"cm"))
    
    dist_plot
    
    pdf(paste0(output_folder,"fig2.pdf"),height=2.5,width=7)
      dist_plot
    dev.off()  

    
# Fig 3 - regression lines for strct--------
    
  rm(list=ls())
  time1<-Sys.time()
  source("")
  options(warn=1) 
  
  # input/output
    output_folder <- ""
    spatial_seeds <- read.csv("")
    prefire<-fread("", colClasses = list(factor=c("fire","GEDI_beam","region","typeForest")))
    prefire <- as.data.frame(prefire)
    prefire <- prefire[prefire$strct_RH_100>=5,]
    
  # strct_plot 
    strct_vars <- c("strct_tAGBD","strct_RH_98","strct_nmode","strct_mPAI_sparse_b10")
    prefire_subset <- prefire[,c("mtbs_dnbrOW",strct_vars)]
    tmp_sc <- norm_transform_func(prefire_subset,cols_used=c(2:ncol(prefire_subset)))
    tmp_sc[,strct_vars] <- scale(tmp_sc[,strct_vars])
    tmp_or <- outlier_removal_cols_withheld(tmp_sc,z=3,cols_witheld=c(1))
    tmp_sample <- tmp_or[sample(1:nrow(tmp_or),200000),]
    tmp_long <- tmp_sample %>% pivot_longer(!mtbs_dnbrOW,names_to='metric',values_to='value')
    tmp_long <- tmp_long[complete.cases(tmp_long),]
    tmp_long$metric <- as.factor(tmp_long$metric)
    tmp_long$metric <- factor(tmp_long$metric, levels = (strct_vars))
    
    fig3_colors<-c("#E31A1C","#39B600","gold1","#6A3D9A")
    
    size_tmp <- 8
    
    strct_GAM <- ggplot(tmp_long, aes(x = value,  y = mtbs_dnbrOW)) +
                    stat_smooth(aes(color = metric,group=metric,fill=metric),method="gam", alpha=0.1,formula = y ~ s(x, bs = "cs")) +
                    scale_color_manual(values=fig3_colors) +
                    scale_fill_manual(values=fig3_colors) +
                    annotate("text", x = -3, y = 470, label = "(a) GAM",hjust=0,size=size_tmp-5) +
                    xlab("") + #ylim(250,470) +
                    my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") + #ylab("Fire severity (dNBR)") +
                    theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_text(size=size_tmp-1,hjust=1),
                          axis.text.x = element_blank(),
                          plot.margin = unit(c(0.2,0.2,0,0.2), 'lines')) +
                      guides(colour="none",fill="none")

    strct_quad <- ggplot(tmp_long, aes(x = value,  y = mtbs_dnbrOW)) +
                    stat_smooth(aes(color = metric,group=metric,fill=metric),method = "lm", alpha=0.1,formula = y ~ x + I(x^2)) +
                    scale_color_manual(values=fig3_colors,labels = c("AGBD","RH98","nmode",expression(mPAI['10m']))) +
                    scale_fill_manual(values=fig3_colors,labels = c("AGBD","RH98","nmode",expression(mPAI['10m']))) +
                    annotate("text", x = -3, y = 470, label = "(b) spGLMM",hjust=0,size=size_tmp-5) +#ggtitle ("quadratic") + 
                    xlab("Fuel structure (std. values)") +
                    my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") + 
                    theme(legend.position=c(.99,0.01),
                      axis.text.x = element_text(size=size_tmp-1,vjust=1),
                      axis.text.y = element_text(size=size_tmp-1,hjust=1),
                      axis.title.x = element_text(hjust=0.5,size=size_tmp+1),
                      axis.title.y = element_blank(),
                      legend.title=element_blank(),
                      legend.key.height = unit(1,"line"),
                      legend.key.size = unit(2,"lines"),
                      legend.text=element_text(size=size_tmp),
                      legend.margin = ggplot2::margin(2,4,2,2),
                      legend.background = element_rect(colour = 'grey95',fill='grey95', linetype='solid'),
                      legend.key = element_rect(fill = NA, colour = "grey95"),
                      legend.spacing.x = unit(0.2, 'cm'),
                      legend.justification=c(1,0),
                      plot.margin = unit(c(0,0.2,0.2,0.2), 'lines'))

    pt_plot <- egg::ggarrange(strct_GAM,strct_quad,ncol = 1) 
    pt_plot
    
    pt_plot2 <- annotate_figure(pt_plot,left = text_grob("Fire severity (dNBR)",rot = 90, vjust=0.5,size=size_tmp+1))
    
    pdf(paste0(output_folder,"fig3.pdf"),height=5,width=3.25)
    pt_plot2
    dev.off() 
    
# Fig 4 - bivar catepillars -------      
 
  rm(list=ls())
  time1<-Sys.time()
  source("")
  options(warn=1)
  
  output_folder <- ""
  size_tmp<-9
  
  # input/output
    bivar_results <- read.csv("...4_bi_inter_multi.csv")
    strct_vars <- c("strct_tAGBD","strct_RH_98","strct_nmode","strct_mPAI_sparse_b10")
    env_vars <- c("topo_slope","w_vpd_d_max","w_ET_b4_mean","w_vs_d_max") 

    tmp <- bivar_results[bivar_results$resp_var=="mtbs_dnbrOW" & bivar_results$strct_var %in% strct_vars & bivar_results$env_var %in% env_vars,]
    tmp$strct_var <- as.factor(tmp$strct_var)
    tmp$strct_var <- factor(tmp$strct_var, levels = rev(strct_vars),labels=rev(c("Biomass","Canopy height","Layering","Ladder fuels"))) 
    tmp$env_var <- as.factor(tmp$env_var)
    levels(tmp$env_var) <- c("Slope","ET","VPD","Wind speed") 
 
    fg <- "RIs_"
    tmp3 <- tmp[,c("strct_var","env_var",names(tmp)[names(tmp)%ilike%fg])]
    names(tmp3) <- gsub(fg,"",names(tmp3))

    fig4_colors<-c("#4664DA","chartreuse3","darkorange1","#B71D02")
    
    strct_plot <- ggplot(tmp3,aes(x=mCI_strct,y=strct_var, shape=strct_sig,fill=strct_sig,group=env_var,color=env_var)) + 
                    geom_pointrange(aes(xmin = lCI_strct, xmax = hCI_strct),linewidth=.5,alpha=1,size=.2,position=position_dodge(width=-0.4)) +

                     scale_shape_manual(values = c(21,19)) +
                     scale_fill_manual(values = c('grey95',"red")) +
                     scale_colour_manual(values=rev(fig4_colors)) +
                     geom_vline(xintercept=0,linewidth=.3,linetype=2,color="grey20") + 
                     xlab("") + 
                    ggtitle("(a) Fuel structure") +
                    my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
                    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=.9,size = size_tmp-2),
                      axis.text.y = element_text(size = size_tmp,hjust=1,colour="black",face = "bold"), 
                      axis.title.y=element_blank(),
                      plot.title = element_text(size=size_tmp,vjust=1,hjust=0),
                      plot.margin = unit(c(0.2,-0.1,0.2,0.2), 'lines'),
                      legend.position='none',
                      legend.title=element_blank(),
                      strip.text = element_text(size = size_tmp),
                      legend.text.align = 0,
                      legend.key = element_rect(linewidth = 1),
                      legend.spacing.x = unit(.1,"cm"))
   
      env_plot <-   ggplot(tmp3,aes(x=mCI_env,y=strct_var, shape=env_sig,fill=env_sig,group=env_var,color=env_var)) + 
                      geom_pointrange(aes(xmin = lCI_env, xmax = hCI_env),linewidth=.5,alpha=1,size=.2,position=position_dodge(width=-0.4)) +
                      scale_shape_manual(values = c(21,19)) +
                      scale_fill_manual(values = c('grey95',"red")) +
                      scale_colour_manual(values=rev(fig4_colors)) +
                      geom_vline(xintercept=0,linewidth=.3,linetype=2,color="grey20") + 
                      xlab(expression("partial effect (std. "*beta[1]*")")) +
                      ggtitle("(b) Topo-weather") +
                      my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
                      theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=.9,size = size_tmp-2),
                        axis.text.y = element_blank(),
                        axis.title.y=element_blank(),
                        axis.title.x=element_text(family="Helvetica",size=size_tmp),
                        plot.title = element_text(size=size_tmp,vjust=1,hjust=0),
                        plot.margin = unit(c(0.2,-0.1,0.2,-0.1), 'lines'),
                        legend.position='none',
                        legend.title=element_blank(),
                        strip.text = element_text(size = size_tmp),
                        legend.text.align = 0,
                        legend.key = element_rect(linewidth = 1),
                        legend.spacing.x = unit(.1,"cm"))
      
      int_plot <-   ggplot(tmp3,aes(x=mCI_int,y=strct_var, shape=int_sig,fill=int_sig,group=env_var,color=env_var)) + 
                      geom_pointrange(aes(xmin = lCI_int, xmax = hCI_int),linewidth=.5,alpha=1,size=.2,position=position_dodge(width=-0.4)) +
                      scale_shape_manual(values = c(21,19)) +
                      scale_fill_manual(values = c('grey95',"red")) +
                      scale_colour_manual(values=rev(fig4_colors)) +
                      geom_vline(xintercept=0,linewidth=.3,linetype=2,color="grey20") + 
                      xlab("") + 
                      ggtitle("(c) Interaction") +
                      my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
                      guides(shape="none",fill = "none") + 
                      theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=.9,size = size_tmp-2),
                        axis.text.y = element_blank(),
                        axis.title.y=element_blank(),
                        axis.title.x=element_text(family="Helvetica",size=size_tmp),
                        plot.title = element_text(size=size_tmp,vjust=1,hjust=0),
                        plot.margin = unit(c(0.2,1,1,-0.1), 'lines'),
                        legend.margin=ggplot2::margin(1,1,1,0),
                        legend.title=element_blank(),
                        strip.text = element_text(size = size_tmp),
                        legend.text.align = 0,
                        legend.text=element_text(size=size_tmp),
                        legend.key.height = unit(1,"line"),
                        legend.key.size = unit(1.5,"lines"),
                        legend.key = element_rect(linewidth = 1),
                        legend.spacing.x = unit(.1,"cm"))



      pt_plot <- egg::ggarrange(strct_plot,env_plot,int_plot,ncol = 3)
      
    pdf(paste0(output_folder,"fig4.pdf"),height=2.5,width=8)
        print(pt_plot)
    dev.off() 


# Fig 5 - topo and weather subsets ------
  
  rm(list=ls())
  time1<-Sys.time()
  source("")
  options(warn=1)
  
  # input
  output_folder <- ""
  uni_inter <- read.csv("...3_uni_inter_subsets.csv")
  unique(uni_inter$subset)
  
  strct_vars <- c("strct_tAGBD","strct_RH_98","strct_nmode","strct_mPAI_sparse_b10")
  subsets <- c("topo_slope","w_vpd_d_max","w_ET_b4_mean","w_vs_d_max") #"w_th_d_mean"
  
  
  uni_inter <- uni_inter[uni_inter$tmp_pred%in%strct_vars & uni_inter$subset%in% subsets,]
  uni_inter$tmp_pred <- as.factor(uni_inter$tmp_pred)
  uni_inter$tmp_pred <- factor(uni_inter$tmp_pred, levels = rev(strct_vars),labels=rev(c("Biomass","Canopy height","Layering","Ladder fuels")))
  uni_inter$subset <- factor(uni_inter$subset, levels = rev(subsets),labels=rev(c("Slope","VPD","ET","Wind speed"))) # "Wind direction"
  
  size_tmp <- 9
  
  Fig5a_colors <- c(brewer.pal(n = 9, name = "Reds")[5:9],"black")
  Fig5b_colors <- c(brewer.pal(n = 9, name = "Blues")[5:9],"black")
  Fig5c_colors <- c(brewer.pal(n = 9, name = "Greens")[5:9],"black")
  Fig5d_colors <- c(brewer.pal(n = 9, name = "Purples")[5:9],"black")
  
  
  #5a - slopes
  
  tmp_slope <- uni_inter[uni_inter$subset=="Slope",]
  tmp_slope$quant <- factor(tmp_slope$quant, levels = c("low","mid-low","mid","mid-high","high","all"),labels=c("< 8","8 - 13","13 - 17","17 - 20","> 20","all")) 
  
  slopes_plot <- ggplot(tmp_slope,aes(x=mCI,y=tmp_pred, shape=mCI_sig,fill=mCI_sig,group=factor(quant),colour=factor(quant))) + 
                  geom_pointrange(aes(xmin = lCI, xmax = hCI),alpha=1,size=.3,position=position_dodge(width=0.5)) +
                  scale_shape_manual(values = c(21,19)) +
                  scale_colour_manual(values=Fig5a_colors) +
                  scale_fill_manual(values = c('grey98',"red")) +
                  geom_vline(xintercept=0,linetype=2,color="grey20") +
                  guides(shape = "none",fill="none",colour = guide_legend(reverse=T)) +
                  labs(colour="Slope°") +
                  xlab(expression("std. "*beta[1])) + ylab("") +
                  my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
                  theme(axis.text.x = element_text(angle = 45, vjust = .95, hjust=.9,size = size_tmp-1),
                        axis.text.y = element_text(size = size_tmp+3,face = "bold",hjust=1),
                        axis.title.y=element_blank(),
                        strip.background = element_rect(colour="white", fill="white"),
                        axis.title.x=element_text(family="Helvetica",size=size_tmp+3),
                        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
                        legend.margin=ggplot2::margin(0,1,0,0),
                        legend.title = element_text(hjust=0,size=size_tmp+3),
                        legend.position='right',
                        strip.text = element_text(size = size_tmp+1,face = "bold"),
                        legend.text = element_text(hjust=0,size=size_tmp+2),
                        legend.key = element_rect(linewidth = 1),
                        legend.spacing.x = unit(.1,"cm"))
  
  #5b - VPD
  
  tmp_VPD <- uni_inter[uni_inter$subset=="VPD",]
  tmp_VPD$quant <- factor(tmp_VPD$quant, levels = c("low","mid-low","mid","mid-high","high","all"),labels=c("< 2.6","2.6 - 2.9","2.9 - 3.2","3.2 - 3.8","> 3.8","all")) 
  
  VPD_plot <- ggplot(tmp_VPD,aes(x=mCI,y=tmp_pred, shape=mCI_sig,fill=mCI_sig,group=factor(quant),colour=factor(quant))) + 
    geom_pointrange(aes(xmin = lCI, xmax = hCI),alpha=1,size=.3,position=position_dodge(width=0.5)) +
    scale_shape_manual(values = c(21,19)) +
    scale_colour_manual(values=Fig5b_colors) +
    scale_fill_manual(values = c('grey98',"red")) +
    geom_vline(xintercept=0,linetype=2,color="grey20") +
    guides(shape = "none",fill="none",colour = guide_legend(reverse=T)) +
    labs(colour="VPD (kPa)") +
    xlab(expression("std. "*beta[1])) + ylab("") +
    my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
    theme(axis.text.x = element_text(angle = 45, vjust = .95, hjust=.9,size = size_tmp-1),
          axis.text.y = element_blank(),
          axis.title.y=element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          axis.title.x=element_text(family="Helvetica",size=size_tmp+3),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
          legend.margin=ggplot2::margin(0,1,0,0),
          legend.position='right',
          strip.text = element_text(size = size_tmp+1,face = "bold"),
          legend.title = element_text(hjust=0,size=size_tmp+3),
          legend.text = element_text(hjust=0,size=size_tmp+2),
          legend.key = element_rect(linewidth = 1),
          legend.spacing.x = unit(.1,"cm"))
  
  #5c - ET
  
  tmp_ET <- uni_inter[uni_inter$subset=="ET",]
  tmp_ET$quant <- factor(tmp_ET$quant, levels = c("low","mid-low","mid","mid-high","high","all"),labels=c("< .4",".4 - .6",".6 - .7",".7 - 1.0","> 1.0","all")) 
  
  ET_plot <- ggplot(tmp_ET,aes(x=mCI,y=tmp_pred, shape=mCI_sig,fill=mCI_sig,group=factor(quant),colour=factor(quant))) + 
    geom_pointrange(aes(xmin = lCI, xmax = hCI),alpha=1,size=.3,position=position_dodge(width=0.5)) +
    scale_shape_manual(values = c(21,19)) +
    scale_colour_manual(values=Fig5c_colors) +
    scale_fill_manual(values = c('grey98',"red")) +
    geom_vline(xintercept=0,linetype=2,color="grey20") +
    guides(shape = "none",fill="none",colour = guide_legend(reverse=T)) +
    labs(colour="ET (mm/8days)") +
    xlab(expression("std. "*beta[1])) + ylab("") +
    my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
    theme(axis.text.x = element_text(angle = 45, vjust = .95, hjust=.9,size = size_tmp-1),
          axis.text.y = element_blank(),
          axis.title.y=element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          axis.title.x=element_text(family="Helvetica",size=size_tmp+3),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
          legend.margin=ggplot2::margin(0,1,0,0),
          legend.position='right',
          strip.text = element_text(size = size_tmp+1,face = "bold"),
          legend.title = element_text(hjust=0,size=size_tmp+3),
          legend.text = element_text(hjust=0,size=size_tmp+2),
          legend.key = element_rect(linewidth = 1),
          legend.spacing.x = unit(.1,"cm"))
  
  #5d - wind
  
  tmp_ws <- uni_inter[uni_inter$subset=="Wind speed",]
  tmp_ws$quant <- factor(tmp_ws$quant, levels = c("low","mid-low","mid","mid-high","high","all"),labels=c("< 6.9","6.9 - 8.3","8.3 - 9.3","9.3 - 11.6","> 11.6","all")) 
  
  ws_plot <- ggplot(tmp_ws,aes(x=mCI,y=tmp_pred, shape=mCI_sig,fill=mCI_sig,group=factor(quant),colour=factor(quant))) + 
    geom_pointrange(aes(xmin = lCI, xmax = hCI),alpha=1,size=.3,position=position_dodge(width=0.5)) +
    scale_shape_manual(values = c(21,19)) +
    scale_colour_manual(values=Fig5d_colors) +
    scale_fill_manual(values = c('grey98',"red")) +
    geom_vline(xintercept=0,linetype=2,color="grey20") +
    guides(shape = "none",fill="none",colour = guide_legend(reverse=T)) +
    labs(colour="Wind Speed (m/s)") +
    xlab(expression("std. "*beta[1])) + ylab("") +
    my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
    theme(axis.text.x = element_text(angle = 45, vjust = .95, hjust=.9,size = size_tmp-1),
          axis.text.y = element_blank(),
          axis.title.y=element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          axis.title.x=element_text(family="Helvetica",size=size_tmp+3),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
          legend.margin=ggplot2::margin(0,1,0,0),
          legend.position='right',
          strip.text = element_text(size = size_tmp+1,face = "bold"),
          legend.title = element_text(hjust=0,size=size_tmp+3),
          legend.text = element_text(hjust=0,size=size_tmp+2),
          legend.key = element_rect(linewidth = 1),
          legend.spacing.x = unit(.1,"cm"))
  
  
  fig5 <- egg::ggarrange(slopes_plot,VPD_plot,ET_plot,ws_plot,ncol = 4)
  
  pdf(paste0(output_folder,"Fig5.pdf"),height=3,width=12)
  fig5
  dev.off()  


# Fig 6 cross-fire gradients -------
    
  # environment
    rm(list=ls())
    source("")
    options(warn=1)
    output_folder <-""
    
  # data
    univar_results <- read.csv("...2_uni_intra.csv")
    univar_results <- univar_results[univar_results$mCI_sig=="sig",]
    cwst_results <- read.csv("...fires_cwst.csv")
    cwst_results$TSeas_mean <- cwst_results$TSeas_mean/10
    cwst_results$AI_mean <- scale(cwst_results$AI_mean)
    strct_vars <- c("strct_tAGBD","strct_RH_98","strct_mPAI_sparse_b10","strct_nmode")

  # wrangling      
    univar_results2 <- univar_results[univar_results$pred%in% strct_vars & univar_results$mCI_sig=="sig",]

    size_tmp <- 8
    fig6_colors <- c("#6077c0","orange","chartreuse3","brown2")
    
  # panel a 
      strct_var<-"strct_tAGBD" ; gradient_var<-"topo_slope_mean"
      univar_results3 <- univar_results2[univar_results2$pred==strct_var,]
      test2<- cwst_results[,c("fire","Abbv",gradient_var)]
      names(test2)[3] <- "tmp_gradient"
      test3 <- merge(univar_results3,test2,by="fire")
      lmtest<-lm(mCI~tmp_gradient +I(tmp_gradient^2),data=test3)
      r2a <- paste0("~R^2 == ", round(summary(lmtest)$r.squared,2))
      
        plot_a <- ggplot(test3,aes(x=tmp_gradient,y=mCI)) + 
          geom_pointrange(aes(ymin = lCI, ymax = hCI),alpha=1,size=.3,position=position_dodge(width=0.5)) +
          scale_shape_manual(values = c(21,19)) +
          geom_label_repel(aes(label = Abbv),size = 2.5) + 
          stat_smooth(method = "lm",formula = y ~ x + I(x^2),fill=fig6_colors[1],colour=fig6_colors[1]) +
          geom_label(aes(x = -Inf, y = Inf, label = "(a)"),hjust=-0.2,label.size = NA,vjust=1,fill="grey95",colour="black",show.legend = FALSE) +
          my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
          theme(legend.position="none",
                axis.title.x = element_text(size=size_tmp+2),
                axis.title.y = element_text(margin = ggplot2::margin(0,0,0,0),size=size_tmp+2),
                axis.text.x = element_text(size=size_tmp-2,vjust=1),
                axis.text.y = element_text(hjust=1,size=size_tmp-2),
                plot.margin = unit(c(0.2, 0, 0.2, 0.2), "cm")) +
          xlab("Mean slope (°)") + ylab(expression('AGBD~dNBR'~~beta[1])) + 
          geom_label(aes(x = Inf, y = Inf, hjust=1.2,vjust=1.2,label = r2a), size=size_tmp-5,parse=T,colour="white",fill = fig6_colors[1])
        
        
        plot_a_marg <- ggMarginal(plot_a,  fill = fig6_colors[1],alpha=0.5,size=8)

    # panel b  
      strct_var<-"strct_tAGBD" ; gradient_var<-"TSeas_mean"
      univar_results3 <- univar_results2[univar_results2$pred==strct_var,]
      test2<- cwst_results[,c("fire","Abbv",gradient_var)]
      names(test2)[3] <- "tmp_gradient"
      test3 <- merge(univar_results3,test2,by="fire")
      lmtest<-lm(mCI~tmp_gradient +I(tmp_gradient^2),data=test3)
      r2b <- paste0("~R^2 == ", round(summary(lmtest)$r.squared,2)) 
      
      plot_b <- ggplot(test3,aes(x=tmp_gradient,y=mCI)) + 
          geom_pointrange(aes(ymin = lCI, ymax = hCI),alpha=1,size=.1,linewidth=.2,position=position_dodge(width=0.5)) +
          scale_shape_manual(values = c(21,19)) +
          stat_smooth(method = "lm",formula = y ~ x + I(x^2),fill=fig6_colors[2],colour=fig6_colors[2]) +
          geom_label(aes(x = -Inf, y = Inf, label = "(b)"),hjust=-0.2,label.size = NA,vjust=1,fill="grey95",colour="black",show.legend = FALSE) +
          my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
          theme(legend.position="none",
                axis.title.x = element_text(size=size_tmp+2),
                axis.title.y = element_text(margin = ggplot2::margin(0,0,0,0),size=size_tmp+2),
                axis.text.x = element_text(size=size_tmp-2,vjust=1),
                axis.text.y = element_text(hjust=1,size=size_tmp-2),
                plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
        xlab(expression(T["seas"])) + ylab(expression('AGBD~dNBR'~~beta[1])) + 
          geom_label(aes(x = Inf, y = Inf, hjust=1.1,vjust=1.1,label = r2b),size=size_tmp-5,parse=T,fill = fig6_colors[2])
        
      plot_b_marg <- ggMarginal(plot_b,  fill = fig6_colors[2],alpha=0.5,size=8)


    # panel c
      strct_var<-"strct_tAGBD" ; gradient_var<-"DGC_mean"
      univar_results3 <- univar_results2[univar_results2$pred==strct_var,]
      test2<- cwst_results[,c("fire","Abbv",gradient_var)]
      names(test2)[3] <- "tmp_gradient"
      test3 <- merge(univar_results3,test2,by="fire")
      lmtest<-lm(mCI~tmp_gradient +I(tmp_gradient^2),data=test3)
      r2c <- paste0("~R^2 == ", round(summary(lmtest)$r.squared,2)) 
        
      plot_c <- ggplot(test3,aes(x=tmp_gradient,y=mCI)) + 
          geom_pointrange(aes(ymin = lCI, ymax = hCI),alpha=1,size=.1,linewidth=.2,position=position_dodge(width=0.5)) +
          scale_shape_manual(values = c(21,19)) +
          stat_smooth(method = "lm",formula = y ~ x + I(x^2),fill=fig6_colors[3],colour=fig6_colors[3]) +
          geom_label(aes(x = -Inf, y = Inf, label = "(c)"),hjust=-0.2,label.size = NA,vjust=1,fill="grey95",colour="black",show.legend = FALSE) +
          my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
          theme(legend.position="none",
                axis.title.x = element_text(size=size_tmp+2),
                axis.title.y = element_text(margin = ggplot2::margin(0,0,0,0),size=size_tmp+2),
                axis.text.x = element_text(size=size_tmp-2,vjust=1),
                axis.text.y = element_text(hjust=1,size=size_tmp-2),
                plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
          xlab(expression(GPP['cum'])) + ylab(expression('AGBD~dNBR'~~beta[1])) + 
          geom_label(aes(x = Inf, y = Inf, hjust=1.1,vjust=1.1,label = r2c),size=size_tmp-5, parse=T,fill = fig6_colors[3])

      plot_c_marg <- ggMarginal(plot_c,  fill = fig6_colors[3],alpha=0.5,size=8)

            
    # panel d  
      strct_var<-"strct_tAGBD" ; gradient_var<-"AI_mean"
      univar_results3 <- univar_results2[univar_results2$pred==strct_var,]
      test2<- cwst_results[,c("fire","Abbv",gradient_var)]
      names(test2)[3] <- "tmp_gradient"
      test3 <- merge(univar_results3,test2,by="fire")
      lmtest<-lm(mCI~tmp_gradient +I(tmp_gradient^2),data=test3)
      r2d <- paste0("~R^2 == ", round(summary(lmtest)$r.squared,2)) 
      
      plot_d <- ggplot(test3,aes(x=tmp_gradient,y=mCI)) + 
          geom_pointrange(aes(ymin = lCI, ymax = hCI),alpha=1,size=.1,linewidth=.2,position=position_dodge(width=0.5)) +
          scale_shape_manual(values = c(21,19)) +
          stat_smooth(method = "lm",formula = y ~ x + I(x^2),fill=fig6_colors[4],colour=fig6_colors[4]) +
          geom_label(aes(x = -Inf, y = Inf, label = "(d)"),hjust=-0.2,label.size = NA,vjust=1,fill="grey95",colour="black",show.legend = FALSE) +
          my_theme(mycolor_background="grey95",mycolor_text="black",plot_background="white") +
          theme(legend.position="none",
                axis.title.x = element_text(size=size_tmp+2),
                axis.title.y = element_text(margin = ggplot2::margin(0,0,0,0),size=size_tmp+2),
                axis.text.x = element_text(size=size_tmp-2,vjust=1),
                axis.text.y = element_text(hjust=1,size=size_tmp-2),
                plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
                   xlab("AI") + ylab(expression('AGBD~dNBR'~~beta[1])) + 
          geom_label(aes(x = Inf, y = Inf, hjust=1.1,vjust=1.1,label = r2d),size=size_tmp-5, parse=T,fill = fig6_colors[4])
        
      plot_d_marg <- ggMarginal(plot_d, fill = fig6_colors[4],alpha=0.5,size=8)

            
    # gather subplots and print
        pt_plot <- ggpubr::ggarrange(plot_a_marg,vert_col,widths = c(3,1),ncol = 2)
    
        pdf(paste0(output_folder,"fig6.pdf"),height=5,width=8)
          pt_plot
        dev.off()  
      