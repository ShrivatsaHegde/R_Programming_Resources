setwd(R"(C:\Dessertation\metabric_analysis)")
list.files()

Sam = read.table('HER2Sample.txt', sep='\t',header = 1)

library(ggpubr)
library(RColorBrewer)

j <- ggboxplot(Sam, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2",
               color = "INFERRED_MENOPAUSAL_STATE",
               fill = "INFERRED_MENOPAUSAL_STATE",
               add = "jitter", shape = "INFERRED_MENOPAUSAL_STATE") +
  stat_compare_means(method = "t.test") +
  scale_fill_brewer(palette = "Set2", type = "qual", na.value = "transparent", name = "")
# Save the plot
ggsave("PTK2HER2Ttestboxplot.png",
       plot = j,
       width = 10,
       height = 10,
       units = "in",
       dpi = 500)

f=ggboxplot(Sam, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2",
            color = "INFERRED_MENOPAUSAL_STATE", palette = "jco")+
  stat_compare_means(method = "anova", label.y = 13)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")
ggsave("PTK2boxplotANOVAHER2.png",
       plot = f,
       width = 10,
       height = 10,
       units = "in",
       dpi = 400)
HR=read.table('HR+Sample.txt', sep='\t',header = 1)

b=ggboxplot(HR, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2",
            color = "INFERRED_MENOPAUSAL_STATE",
            fill = "INFERRED_MENOPAUSAL_STATE",
            add = "jitter", shape = "INFERRED_MENOPAUSAL_STATE") +
  stat_compare_means(method = "t.test") +
  scale_fill_brewer(palette = "Set2", type = "qual", na.value = "transparent", name = "")
# Save the plot
ggsave("PTK2HR+Ttestboxplot.png",
       plot = b,
       width = 10,
       height = 10,
       units = "in",
       dpi = 500)
m=ggboxplot(HR, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2",
            color = "INFERRED_MENOPAUSAL_STATE", palette = "jco")+
  stat_compare_means(method = "anova", label.y = 13)+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")
ggsave("PPTK2HR+boxplotANOVA.png",
       plot = m,
       width = 10,
       height = 10,
       units = "in",
       dpi = 500)



TNBC=read.table('TNBCSample.txt', sep='\t',header = 1)


o=ggboxplot(TNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2",
            color = "INFERRED_MENOPAUSAL_STATE",
            fill = "INFERRED_MENOPAUSAL_STATE",
            add = "jitter", shape = "INFERRED_MENOPAUSAL_STATE") +
  stat_compare_means(method = "t.test") +
  scale_fill_brewer(palette = "Set2", type = "qual", na.value = "transparent", name = "")
# Save the plot
ggsave("PTK2TNBCTtestboxplot.png",
       plot = o,
       width = 10,
       height = 10,
       units = "in",
       dpi = 500)
x=ggboxplot(TNBC, x = "INFERRED_MENOPAUSAL_STATE", y = "PTK2",
            color = "INFERRED_MENOPAUSAL_STATE", palette = "jco")+
  stat_compare_means(method = "anova", label.y = 13)+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")
ggsave("PPTK2TNBCboxplotANOVA.png",
       plot = x,
       width = 10,
       height = 10,
       units = "in",
       dpi = 500)