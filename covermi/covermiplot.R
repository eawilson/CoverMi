library(lattice)

coverageplot<- function(dataframe, flip=TRUE, shrink=TRUE, ...) {

    df<-dataframe
    df$component<-factor(df$component, levels=c("minimum", "coverage", "amplicon", "exon", "exon_number", "variants", "other"), ordered=TRUE)
    df$name<-factor(df$name, levels=sort(levels(df$name)), ordered=TRUE)
    df$y[df$y<1]<-1
    if (shrink) {
        df$x<-df$x-df$adjustment
        yisactuallyx<-df$component=="amplicon" | df$component=="exon" | df$component=="other"
        df$y[yisactuallyx]<-df$y[yisactuallyx]-df$adjustment[yisactuallyx]
    }
    if (flip) {
        negative<-df$strand=="-"
        df$x[negative]<-df$x[negative]*(-1)
        df<-df[order(df$x, df$component),]
    }

    xyplot(y~x|name, df, group=component, scale=list(x=list(relation="free", draw=!shrink), y=list(alternating=FALSE, at=c(0, 1, 2, 3, 4), labels=c("0", "10", "100", "1000", "10,000"))), xlab="Exon", ylab="Read Depth (log scale)", as.table=TRUE, ..., 
        panel=panel.superpose,
        panel.groups=function(x, y, ..., group.number) {
            if (group.number==1) {
                #minimum
                panel.abline(h=log10(y), lty="dotted")
            }
            if (group.number==2) {
                #coverage
                panel.xyplot(abs(x), log10(y), type="s", lwd=2)
            }
            if (group.number==3) {
                #amplicon
                panel.rect(abs(x[seq(1, length(x), 2)]), 0, y[seq(1, length(y), 2)], -0.04, fill="black") 
                if (length(x) > 1) {
                    panel.rect(abs(x[seq(2, length(x), 2)]), 0.04, y[seq(2, length(y), 2)], 0, fill="black")
                }
            }
            if (group.number==4) {
                #exon
                panel.abline(h=-0.2, lwd=2, col="black")
                panel.rect(abs(x), -0.3, abs(y), -0.1, fill="dimgray", border="black", lwd=1.5)
            }
            if (group.number==5) {
                #exon number
                step<-seq(1, length(x), (length(x)%/%20)+1)
                panel.text( abs(x)[step], -0.3, labels=y[step], pos=1 )
            }
            if (group.number==6) {
                #variant
                panel.xyplot(abs(x), log10(y), type="p", col="black", pch=4)
            }
            if (group.number==7) {
                #other
                panel.rect(abs(x), 0.8, abs(y), 1)
            }

        },
        prepanel=function(x, y, ...) list(xlim=c(abs(min(x)), abs(max(x))), ylim=c(-0.5, 4)),
    )
}     

df<-read.delim("INPUT")
pdf(file="OUTPUT", paper="a4r", width=11, height=7)
print(coverageplot(df, main="NAME", layout=c(1,2)))
dev.off()
