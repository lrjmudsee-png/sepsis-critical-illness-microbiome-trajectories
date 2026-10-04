args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
figure_rows <- read.csv(args[1],check.names=FALSE)
stopifnot(nrow(figure_rows)==5L, all(c('label','estimate','low','high','type') %in% names(figure_rows)))
# Pure layout repair from frozen figure-source values; no statistical rerun.
pdf(args[2],width=10,height=5.5)
par(mar=c(5,12,3,2))
table <- figure_rows[rev(seq_len(nrow(figure_rows))),,drop=FALSE]
limits <- range(c(table$low,table$high,0),finite=TRUE)
pad <- diff(limits)*0.12
plot(table$estimate,seq_len(nrow(table)),xlim=limits+c(-pad,pad),
     ylim=c(0.5,nrow(table)+0.5),yaxt='n',ylab='',
     xlab='Mean early-to-late delta_B',
     pch=ifelse(table$type=='Pooled',18,19),
     cex=ifelse(table$type=='Pooled',1.4,1),
     main='Direction balance change across cohorts')
abline(v=0,lty=2,col='grey50')
segments(table$low,seq_len(nrow(table)),table$high,seq_len(nrow(table)),lwd=2)
axis(2,at=seq_len(nrow(table)),labels=table$label,las=1,cex.axis=0.85)
legend('bottomright',legend=c('Natural cohort','External primary','Pooled'),
       pch=c(19,19,18),col='black',bty='n',cex=0.85)
dev.off()
