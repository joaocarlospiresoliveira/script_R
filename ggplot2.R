library(ggplot2)
library(extrafont)
# font_import()
loadfonts(device="win") 


# Color<-factor(c(1:52), labels=c(family))
{(p <- ggplot(wc, aes(x = gcm, y = means)) +
  geom_point(stat = "identity",color="red", fill="black") +# Parametros basicos do tipo de grafico
    geom_errorbar(aes(ymax = means + sds, ymin = means - sds),
                  position = "dodge") +
    theme_test(                         # Tema do grafico
    base_size = 10,
    base_family = '',
    base_line_size = 10,
    base_rect_size = 1) +
  # geom_label(aes(label = NULL), vjust = 0.5, size = 4)) + #
  # scale_fill_manual(values = rep(1, 28), c('black')) +
  ggtitle('"Correlação entre GCMs - Cluster 1"') +
    theme(plot.title = element_text(hjust=0.5))+
  # ggtitle('"Families of Indicator Species of the ACT"')+
  labs(x = 'GCM', y = "rho Médio") +
  # labs(x = 'Family', y = "Species Number") +
  theme(
    plot.title = element_text(family="",
    size = rel(2),
    vjust = 2,
    lineheight = -3),
    legend.position = "bottom") +
  theme(axis.title  =  element_text(size = 25)) +
  theme(text = element_text(family = "" )) +
  theme(axis.text.y = element_text(size = 14,lineheight = 8 ))+
    theme(axis.text.x = element_text(size = 19,lineheight = 8 )))# windows(w=100, h=100)
(grafico<-p + geom_vline(xintercept = 4, linetype = "dashed", colour = 'red'))}

jpeg(filename = 'test.jpg', height = 25, width =45, res = 400, units = 'cm')
grafico
dev.off()
