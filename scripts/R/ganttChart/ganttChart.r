# https://lpmor22.github.io/ | 2025-06-03

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", dependencies = TRUE)
library("pacman")

p_load(ggplot2, this.path, tidyr)

setwd(dirname(this.path()))

input <- "KhouriousTeam2025.csv"

df <- read.csv(input, header = TRUE)

acts <- c("Ricardo Khouri", "Valdomiro Moitinho", "Maria da Purificacao",
          "Laise de Moraes", "Lia Barbara Arruda", "Cibele Orge",
          "Luciane Amorim", "Lucas Campos", "Camila Brito", "Felipe Torres",
          "Lucas Gentil","Lais Cambui", "Marcos Braz", "Isabele Carvalho",
          "Thaline Mabel Santos", "Victoria Ribeiro", "Kleverton Ribeiro",
          "Debora Magnavita", "Jessica Duarte", "Raphael Borges", "Marina Cucco",
          "Ellen Pimentel", "Luisa Pedrosa", "Lilian Gomes", "Raphaela Andrade",
          "Jose Adriano Silva", "Agatha Morais", "Ana Carolina Lima",
          "Beatriz Vasconcelos", "Giulia Botelho", "Igor Silveira",
          "Joyce Silva", "Icaro Strobel", "Clara Porto", "Dhara Isabella Silva",
          "Laura Barreto", "Jessica Lais Santos", "Caio Galiano"
)

els <- c("PI", "Research Associate", "Ph.D. student", "pre-Ph.D. student",
         "M.Sc. student", "pre-M.Sc. student", "Undergraduated student",
         "Technician"
)

actcols <- c("#FA4749", "#EF8CC8", "#79E16F", "#39A638", "#2C872D",
             "#19ADFF", "#193BFF", "#907AFF")

gantt <- gather(df, "state", "date", 4:5) %>%
    mutate(date = as.Date(date, "%Y.%m.%d"),
           a2 = factor(a2, acts[length(acts):1]),
           a3 = factor(a3, els))

gantt_chart <- ggplot(gantt, aes(date, a2, colour = a3)) +
  geom_line(linewidth = 3) +
  scale_color_manual(values = actcols) +
  theme_classic() +
  scale_x_date(limits = c(as.Date("2015-01-01"), NA),
               expand = c(0, 0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", size = .5),
        panel.grid.minor = element_line(linetype = "dashed", size = .5),
        axis.line = element_line(colour = "#000000"),
        legend.position = c(.21, .3),
        legend.title = element_blank(),
        legend.background = element_rect(size = .5,
                                         linetype = "solid",
                                         colour = "#BEBEBE")) +
  labs(x = NULL, y = NULL)

save_plot(sub("\\.csv$", ".pdf", input), gantt_chart, base_height = 5, base_width = 10)
save_plot(sub("\\.csv$", ".png", input), gantt_chart, base_height = 5, base_width = 10)