

library(ggplot2)
library(viridis) #color pallete
library(reshape2)
library(hrbrthemes)
library(ggsci)


#---------
# for heritability curves (figures a & b)

vG <- c(1,2,5,10,20)
vG <- c(1,2,5,10,20)
vEe <- 5
c_seq <- rev(seq(0,20,0.01))

#compression <- seq(0,1,0.01)
#c_seq <- 20 * rev(compression)
h2A <-  vG[1] / (vG[1] + vEe + c_seq^2/12 )
h2B <-  vG[2] / (vG[2] + vEe + c_seq^2/12 )
h2C <-  vG[3] / (vG[3] + vEe + c_seq^2/12 )
h2D <-  vG[4] / (vG[4] + vEe + c_seq^2/12 )
h2E <-  vG[5] / (vG[5] + vEe + c_seq^2/12 )
dat <- data.frame(c_seq, h2A, h2B, h2C, h2D, h2E)
dat.long <- melt(dat, id="c_seq")
levels(dat.long$variable)  <-  as.character(vG)
colnames(dat.long) <- c("culture_range","V_g","heritability")


plotA <- ggplot(dat.long) +
geom_path(aes(x= culture_range, y= heritability, color=V_g)) +
theme_ipsum()+
scale_x_reverse()+
ylim(0,1)+
scale_color_lancet()  
plotA

ggsave("h2curves_vEe0.eps" ,plotA, device="eps", width=115, height=90, units="mm")
ggsave("h2curves_vEe5.eps" ,plotA, device="eps", width=115, height=90, units="mm")



#---------
# for heritability heatmap (figure c)

vG <- c(1,2,5,10,20)
vG <- c(1,2,5,10,20)
vEe <- 5
c_seq <- rev(seq(0,20,0.01))

ab_seq <- seq(0,20,1)
ab <- matrix(ab_seq, length(ab_seq), length(ab_seq))
ab_mat <- t(ab) - ab
ab_mat[lower.tri(ab_mat)] <- NA
ab_var <- ab_mat^2/12
ab_H2 <- vG[1] / (vG[1] + vEe + ab_var)
ab_H2_long <-  melt(ab_H2)
colnames(ab_H2_long) <- c("trait_min","trait_max","heritability")
ab_H2_long$trait_min <- (ab_H2_long$trait_min - 1)/ ((length(ab_seq)-1)/max(ab_seq)) 
ab_H2_long$trait_max <- (ab_H2_long$trait_max - 1)/ ((length(ab_seq)-1)/max(ab_seq))
ab_H2_long <- ab_H2_long[!is.na(ab_H2_long$heritability),]

plotB <- ggplot(data = ab_H2_long, aes(x =trait_min, y=trait_max, fill=heritability)) +
geom_raster() +
scale_fill_viridis(discrete=F, direction = -1, option="D")+
theme_ipsum()
plotB

ggsave("h2heatmap_vG1_vEe0.bmp" ,plotB, device="bmp", width=140, height=80, units="mm")







