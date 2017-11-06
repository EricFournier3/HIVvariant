#executer avec la commande 
#Rscript MakePlot.R "Complexity Shannon Diversity" "across gp120 region" "stat_results.txt"


args=commandArgs(TRUE)

param_list=list(ev_param=args[1],titre=args[2],data=args[3])

#print(param_list$ev_param)
#print(param_list$titre)
#print(param_list$data)


library(ggplot2)
library(data.table)

plot_f=function(evol_param,Title,data_source){
	print(evol_param)	

	my_data=read.table(data_source,header=T)
	my_data=data.table(my_data)

	g=ggplot(my_data,aes(color=Type))
	g=g+geom_density(aes_string(evol_param),data=subset(my_data,Type=="RECENT"),alpha=0.3)
	g=g+geom_density(aes_string(evol_param),data=subset(my_data,Type=="CHRONIC"),alpha=0.3)
	g=g+facet_wrap(~ Region, ncol=3,scales="free")
	g=g+labs(title=Title,y="Density")
	g=g+scale_color_discrete(name="Infection level") + theme(axis.text.x=element_text(angle=45,hjust=1))
	print(g)
	ggsave(paste(evol_param,".pdf",sep=""))
	
}

title_suffix=param_list$titre

for(i in strsplit(param_list$ev_param," ")[[1]]){
	print(paste("Make plot for",i,sep=":"))	

	if(i=="Shannon"){
		title=paste("Sahnnon index",title_suffix)
	}
	else if(i=="Diversity"){
		title=paste("Nucleotide diversity",title_suffix)
	}

        else{
		title=paste(i,title_suffix)

	}

	plot_f(i,title,param_list$data)
	
}


warnings()
