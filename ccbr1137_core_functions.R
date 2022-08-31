create_top_list<-function(sub_df,mol_num){
  # find top values for each molecule type
  for (mol_id in mol_list){
    # subset for molecule
    sub_df2=subset(sub_df,Molecule.Type==mol_id)
    
    # sort by p value
    sorted_df=sub_df2[order(sub_df2$p.value.of.overlap),]
    
    # add top value to df
    if (exists("top_df")){
      top_df=full_join(top_df,sorted_df[mol_num,])
    } else{
      top_df=sorted_df[mol_num,]
    }
  }
  # add log10 of pvalue
  top_df$signA=top_df$Activation.z.score<0
  top_df$signA=gsub("TRUE",1,top_df$signA)
  top_df$signA=gsub("FALSE",-1,top_df$signA)
  top_df$logpvalue=log(top_df$p.value.of.overlap,100)*as.numeric(top_df$signA)
  return(top_df)
}

create_zscore_pvalue_plot<-function(top_df,mol_num){
  # create plot
  #https://stackoverflow.com/questions/40469757/change-the-shape-of-legend-key-for-geom-bar-in-ggplot2
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                 "#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
  subPalette=cbPalette[length(top_df$Molecule.Type)]
  scale_factor <-2
  
  #create plot without legend
  p1 = top_df %>%
    mutate(Upstream.Regulator = fct_reorder(Upstream.Regulator, Activation.z.score)) %>%
    ggplot(aes(x = Upstream.Regulator))+
    geom_bar(aes(y = Activation.z.score,fill=Molecule.Type),colour="red", stat="identity")+
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor,color=""),shape=20, cex=8)+ # create dummy points so pvalue legend symbol appears
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor),shape=20, cex=8,color="blue")+
    scale_color_manual(name = "log10(pvalue) x \nActivation Z Score",values = c(" " = "blue")) + # add legened for pvalue point
    scale_fill_brewer(direction = -1, palette = "Pastel2")+ # change palette for bars
    coord_flip()+ # flip x and y axis
    scale_y_continuous(sec.axis = sec_axis(~ .*1, #create a secondary "x" axis for the log10(pvalue data)
                                           labels = number_format(scale=scale_factor), 
                                           name="log10(pvalue) x sign of Activation Z Score"))+
    xlab("Upstream Regulator")+
    ylab("Activation Z Score")+
    theme_bw() + # remove grey background
    theme( axis.line.x.top = element_line(color = "blue"),axis.ticks.x.top= element_line(color = "blue"), # change top x axis color to match pvalues
           axis.text.x.top = element_text(color = "blue"),axis.title.x.top=element_text(colour="blue"))+
    theme( axis.line.x.bottom = element_line(color = "red"),axis.ticks.x.bottom= element_line(color = "red"), # change bottom x axis color to match z scores
           axis.text.x.bottom = element_text(color = "red"),axis.title.x.bottom=element_text(colour="red"))+
    ggtitle(paste0("Top ", max(mol_num), " Upstream Regulators by Molecule Type")) + # add title
    theme(plot.title = element_text(size=22, face="bold")) + #  increase title
    theme(legend.key=element_rect(fill=NA)) +# remove grey from pvalue background
    labs(color="log10(pvalue) x \nActivation Z Score",fill="Activation Z Score")+ # change the order of the legend so log10 is on top
    guides(fill  = guide_legend(order = 2),color = guide_legend(order = 1))+
    theme(legend.position = "none") # remove legend
  
  #legend with FILL only
  p2 = top_df %>%
    mutate(Upstream.Regulator = fct_reorder(Upstream.Regulator, Activation.z.score)) %>%
    ggplot(aes(x = Upstream.Regulator))+
    geom_bar(aes(y = Activation.z.score,fill=Molecule.Type),colour="red", stat="identity")+
    scale_fill_brewer(direction = -1, palette = "Pastel2")+ # change palette for bars
    coord_flip()+ # flip x and y axis
    scale_y_continuous(sec.axis = sec_axis(~ .*1, #create a secondary "x" axis for the log10(pvalue data)
                                           labels = number_format(scale=scale_factor), 
                                           name="log10(pvalue) x sign of Activation Z Score"))+
    xlab("Upstream Regulator")+
    ylab("Activation Z Score")+
    theme_bw() + # remove grey background
    theme( axis.line.x.top = element_line(color = "blue"),axis.ticks.x.top= element_line(color = "blue"), # change top x axis color to match pvalues
           axis.text.x.top = element_text(color = "blue"),axis.title.x.top=element_text(colour="blue"))+
    theme( axis.line.x.bottom = element_line(color = "red"),axis.ticks.x.bottom= element_line(color = "red"), # change bottom x axis color to match z scores
           axis.text.x.bottom = element_text(color = "red"),axis.title.x.bottom=element_text(colour="red"))+
    theme(legend.key=element_rect(fill=NA)) +# remove grey from pvalue background
    theme(legend.direction = "vertical", legend.box = "horizontal")+
    theme(legend.title = element_text(size=12,color="red",face="bold"))
  
  # legend for COLOR only
  p3 = top_df %>%
    mutate(Upstream.Regulator = fct_reorder(Upstream.Regulator, Activation.z.score)) %>%
    ggplot(aes(x = Upstream.Regulator))+
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor,color=""),shape=20, cex=8)+ # create dummy points so pvalue legend symbol appears
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor),shape=20, cex=8,color="blue")+
    scale_color_manual(name = "log10(pvalue) x \nsign of Activation Z Score",values = c(" " = "blue")) + # add legened for pvalue point
    coord_flip()+ # flip x and y axis
    scale_y_continuous(sec.axis = sec_axis(~ .*1, #create a secondary "x" axis for the log10(pvalue data)
                                           labels = number_format(scale=scale_factor), 
                                           name="log10(pvalue) x sign of Activation Z Score"))+
    theme_bw() + # remove grey background
    theme(legend.key=element_rect(fill=NA)) +# remove grey from pvalue background
    labs(color="log10(pvalue) x \nActivation Z Score",fill="Activation Z Score")+ # change the order of the legend so log10 is on top
    theme(legend.direction = "vertical", legend.box = "horizontal")+
    theme(legend.title = element_text(size=12,color="blue",face="bold"))
  
  # get two legends
  leg2=get_legend(p2)
  leg3=get_legend(p3)
  
  # combine legends
  leg32 = plot_grid(leg3, leg2,
                    nrow = 2)
  
  final_p = plot_grid(p1,
                      leg32,
                      nrow = 1,
                      align = "h",
                      axis = "t",
                      rel_widths = c(1, 0.3))
  
  
  # save and print
  fpath=paste0(output_dir,sample_id,"_top_",max(mol_num),"_molecules.png")
  ggsave(fpath,final_p)
  print(final_p)
}

main_function<-function(sample_id,merged_df,nmol_include){
  #subset for sample
  sub_df=subset(merged_df,sample==sample_id)
  
  # create top mol list
  mol_num=c(1:nmol_include)
  top_df=create_top_list(sub_df,mol_num)
  
  # create plot
  create_zscore_pvalue_plot(top_df,mol_num)
  
  #create DT
  select_cols=c("Molecule.Type","Upstream.Regulator","Activation.z.score","p.value.of.overlap","logpvalue")
  out_df=top_df[,select_cols]
  out_df$logpvalue=signif(out_df$logpvalue, 4)
  colnames(out_df)=c("Molecule.Type","Upstream.Regulator","Activation.z.score","pvalue","log10.pvalue")
  DT::datatable(out_df,caption = paste0("Top ", max(mol_num), " Upstream Regulators by Molecule Type"))
}