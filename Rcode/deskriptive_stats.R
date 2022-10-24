### Trait-based microbial community

## Descriptive Analysis of Traitdata Table
colSums(!is.na(x[,-c(1:2)]))
sum(!is.na(x[,-c(1:2)]))
# statistical summary
summary(x[,-c(1:3)]) %>% 
  kable() %>% 
  kable_styling()
# calculate mode
calcmode <- function(a) {  
  vector <- unique(a)  
  vector[which.max(tabulate(match(a, vector)))]  
   } 
# replace NA with absent and values with present 
stackplot <- x
stackplot <- stackplot %>% 
  mutate(Aggreg_score = ifelse(is.na(Aggregation_score),"absent","present"),
         B_vit = ifelse(is.na(B_vitamins),"absent","present"),
         Copies_16S = ifelse(is.na(Copies_16S),"absent","present"),
         GC_ = ifelse(is.na(GC_content),"absent","present"),
         Gene_num = ifelse(is.na(Gene_number),"absent","present"),
         Genome_ = ifelse(is.na(Genome_Mb),"absent","present"),
         Gram = ifelse(is.na(Gram_positive),"absent","present"),
         IgA = ifelse(is.na(IgA),"absent","present"),
         Length = ifelse(is.na(Length),"absent","present"),
         Motility = ifelse(is.na(Motility),"absent","present"),
         O2_tol = ifelse(is.na(Oxygen_tolerance),"absent","present"),
         pH_opt = ifelse(is.na(pH_optimum),"absent","present"),
         Salt_opt = ifelse(is.na(Salt_optimum),"absent","present"),
         Spore = ifelse(is.na(Sporulation),"absent","present"),
         Temp_opt = ifelse(is.na(Temp_optimum),"absent","present"),
         Width = ifelse(is.na(Width),"absent","present"))
# filter out unnecessary/ duplicated col
stackplot <- stackplot[,-c(3:4,6:9,13:17)] 
stackplot <- stackplot %>% rename("16S"="Copies_16S","Aggreg"="Aggreg_score","Gene_"="Gene_num")
stackplot <- stackplot %>% mutate(Species = paste0(Genus, " ",Species))
## stacked bar chart
st.plot <- stackplot %>% select(2,3:18) %>% 
  gather(trait, val, -Species) %>% 
  mutate(val = factor(val, levels = c("absent", "present")))

# grouping common traits - long format table
st.plot <- st.plot %>% group_by(trait) %>%
  mutate(count=length(which(val =="present"))/length(which(val =="present"| val == "absent")),
         count = round(count, 4))  
 
# create stacked bar chart
bar <- ggplot(st.plot, aes(x = trait, fill = val)) +
  geom_bar( position = "fill") +
  scale_y_continuous(labels = scales::percent) + 
  labs(fill = 'Condition') +
  ggtitle("Bacterial strain traits") +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "Traits",y="Relative abundance")
bar

# add percentage text label to plot
hist <- bar+geom_text(aes(x = trait,y=count,label=percent(count)),vjust = 0.5,hjust=0.5, size = 2.8) 

# Histogram
library(tidyverse)
hist <- x %>% 
  gather(Attributes, value, 3:18) %>% 
  ggplot(aes(x=value,)) + 
  geom_histogram(fill = "#E69F00", color = "black", bins = 15) +
  geom_density(alpha=.2, fill="#FF6666")+
  ggtitle("Distribution of each trait across 23,058 bacterial species")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Attributes,scales = "free") + #scales = "free_x"
  labs(x = "Trait score", y = "number of bacterial species")
hist


# Pie chart 

install.packages("plotrix")
library(plotrix)

lab <- paste0(round(DF$count/sum(DF$count) * 100, 2), "%")

lp <- pie3D(DF$count,
      col = hcl.colors(length(DF$trait), "Spectral"),
      main="Trait distribution",
      labels = paste0(DF$trait,"  ",lab),
      labelcex = 0.6)
      #explode = 0.25)
lp

