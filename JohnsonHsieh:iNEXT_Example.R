# library(devtools)

# install_github('JohnsonHsieh/iNEXT')

library(iNEXT)
library(ggplot2)
library(gridExtra)

forests<-list(
Dip = c(8, 1, 2, 1, 1, 1),
HillDip = c(50, 3, 6, 5, 2, 1,2),
UppHill = c(31, 2, 3, 4, 1, 1),
OakLau = c(1,1, 1, 1),
Montane = c(8, 2, 2, 1, 3, 1))

out<-iNEXT(forests, q=0, datatype = 'abundance', nboot = 5000)

p1<-ggiNEXT(out, type=1, color.var = 'site')
p2<-ggiNEXT(out, type=2, color.var = 'site')
#p3<-ggiNEXT(out, type=3, facet.var="site")

grid.arrange(p1,p2, ncol = 2)

par(mfrow = c(1,2))
plot(out, type = 1)
plot(out, type = 2)
