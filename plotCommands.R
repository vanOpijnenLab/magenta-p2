#Scatter plot with conditional labeling

ggplot(dat, aes(x=mean.gluc, y=mean.dapto)) +
+     geom_point(shape=1,color="purple")+ geom_smooth(method=lm,se=TRUE) +geom_text(aes(label=ifelse(sub[,3]>.2,as.character(row.names(sub)),''),check_overlap=TRUE),size=3,hjust=-.2,vjust=-.4)+ggtitle("Gene fitness for TIGR4: Glucose v. Daptomycin")+xlab("Mean fitness for genes --glucose")+ylab("Mean fitness for genes --daptomycin")