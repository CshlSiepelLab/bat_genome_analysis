library(dplyr)
library(ggplot2)

# Main figures

total <- read.delim("data/final_all_landscape_total.tab",sep = "\t")

sum <- total %>% group_by(class,species) %>%
  summarise(Frequency = sum(frequency_bp))

sum$species <- factor(sum$species, levels=c("mymy", "meso", "shon","drot","phydis","ajam","jam"))
sum <- sum[sum$Frequency>1000000,]

ggplot(sum, aes(y = species, x = Frequency/1000000, fill = class)) +
  geom_bar(stat="identity",position="stack") +
 # geom_vline(xintercept=6.6, linetype="dashed",
 #            color = "red", size=0.5)+
  xlab("Bases (Mb)") +
  ylab("Bat Species") +
  labs(fill = "Repeat Class")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 22)) +
  ggtitle("Repeat landscape in noctilionid bats") 

ggsave("repeat_landscape.pdf")

div <- read.csv("data/final_all_landscape_div.tab", sep = "\t", header = TRUE)
div <- div %>% mutate(Total = select(., X0percent:X6percent) %>% rowSums(na.rm = TRUE))
divtor <- select(div, c("species","class","Total"))
# sum up
divtor <- divtor %>% group_by(species,class) %>%
  summarise(classtot = sum(Total))

divtor$species <- factor(divtor$species, levels=c("mymy", "meso", "shon","drot","phydis","ajam","jam"))
divtor <- divtor[divtor$classtot>0.01,]
ggplot(divtor, aes(y = species, x = classtot, fill = class)) +
  geom_bar(stat="identity",position="stack") +
  xlab("% of genome") +
  ylab("Bat Species") +
  labs(fill = "Repeat Class")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 22)) +
  ggtitle("Recent repeat expansions in noctilionid bats")

ggsave("repeat_expansions.pdf")
