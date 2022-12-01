# @lpmor22 | https://lpmor22.github.io/

if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

library("tidyverse")

gantt <- read.csv("Gantt_Input.csv", h=T)

acts <- c("Ricardo Khouri", "Valdomiro Moitinho", "Maria da Purificacao", "Laise de Moraes", "Lia Barbara Arruda", "Cibele Orge", "Luciane Amorim", "Lucas Campos", "Camila Brito", "Felipe Torres", "Lucas Gentil", "Lais Cambui", "Marcos Braz", "Isabele Carvalho", "Thaline Mabel Santos", "Victoria Ribeiro", "Kleverton Ribeiro", "Debora Magnavita", "Jessica Duarte", "Raphael Borges", "Marina Cucco", "Ellen Pimentel", "Luisa Pedrosa", "Lilian Gomes", "Raphaela Andrade", "Jose Adriano Silva", "Agatha Morais", "Ana Carolina Lima", "Beatriz Vasconcelos", "Giulia Botelho", "Igor Silveira", "Joyce Silva", "Icaro Strobel", "Clara Porto", "Dhara Isabella Silva", "Laura Barreto", "Jessica Lais Santos", "Caio Galiano")

els <- c("PI", "Research Associate", "Ph.D. student", "pre-Ph.D. student", "M.Sc. student", "pre-M.Sc. student", "Undergraduated student", "Technician")

g.gantt <- gather(gantt, "state", "date", 4:5) %>% mutate(date = as.Date(date, "%Y.%m.%d"), a2=factor(a2, acts[length(acts):1]), a3=factor(a3, els))

actcols <- c( "#fa4749", "#ef8cc8", "#79e16f", "#39a638", "#2c872d", "#19adff", "#193bff", "#907aff")

pdf(file="Gantt_Output.pdf")
ggplot(g.gantt, aes(date, a2, colour = a3)) +
  geom_line(size = 3) +
  scale_color_manual(values=actcols) +
  theme_classic() +
  scale_x_date(limits = c(as.Date("2015-01-01"), NA), expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype="dashed", size=0.5),
        panel.grid.minor = element_line(linetype="dashed", size=0.5),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.21, 0.3),
        legend.title = element_blank(),
        legend.background = element_rect(size=0.5, linetype="solid", colour="grey")) +
  labs(x=NULL, y=NULL)
dev.off()