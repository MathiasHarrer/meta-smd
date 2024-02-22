library(readxl)
library(xlsx)
library(metapsyTools)
library(tidyverse)

# Import original depression database
data = read_excel("data.xlsx")

# Import attrition data database
data.attr = read_excel("data.attr.xlsx")

# Make sure that IDs are created the same way
vars.for.id = c("study", "instrument", "time", "outcome_type")
d2 = calculateEffectSizes(data, vars.for.id = vars.for.id)
d2.attr = calculateEffectSizes(data.attr, vars.for.id = vars.for.id)

# Select attr variable from both datasets
d2 %>% select(.id, attr_arm1, attr_arm2) -> d2
d2.attr %>% select(.id, attr_arm1, attr_arm2) %>% 
  rename(attr_arm1_new = attr_arm1, 
         attr_arm2_new = attr_arm2) -> d2.attr

# Join attr data with big dataset
left_join(d2, d2.attr, by = ".id") -> d2.join
within(d2.join, {
  attr_arm1 = attr_arm1 %>% parse_number(c("", "NA", "NR"))
  attr_arm2 = attr_arm2 %>% parse_number(c("", "NA", "NR"))
  attr_arm1_new = attr_arm1_new %>% parse_number(c("", "NA", "NR"))
  attr_arm2_new = attr_arm2_new %>% parse_number(c("", "NA", "NR"))
}) -> d2.join

# Take available data from both datasets; if both are available and 
# diverge, take the maximum amount.
data.frame(
  attr_arm1 = apply(cbind(d2.join$attr_arm1, d2.join$attr_arm1_new),
                    1,max,na.rm=T) %>% ifelse(.==-Inf, NA, .),
  attr_arm2 = apply(cbind(d2.join$attr_arm2, d2.join$attr_arm2_new),
                    1,max,na.rm=T) %>% ifelse(.==-Inf, NA, .)
  ) -> df



# Make sure that IDs are created the same way
vars.for.id = c("study", "instrument", "time", "outcome_type")
d2 = calculateEffectSizes(data, vars.for.id = vars.for.id)
d2.attr = calculateEffectSizes(data.attr, vars.for.id = vars.for.id)

# Select attr variable from both datasets
d2 %>% select(.id, rand_arm1, rand_arm2) -> d2
d2.attr %>% select(.id, rand_arm1, rand_arm2) %>% 
  rename(rand_arm1_new = rand_arm1, 
         rand_arm2_new = rand_arm2) -> d2.attr

# Join attr data with big dataset
left_join(d2, d2.attr, by = ".id") -> d2.join
within(d2.join, {
  rand_arm1_new = rand_arm1_new %>% parse_number(c("", "NA", "NR"))
  rand_arm2_new = rand_arm2_new %>% parse_number(c("", "NA", "NR"))
}) -> d2.join

# Take available data from both datasets; if both are available and 
# diverge, take the maximum amount.
data.frame(
  rand_arm1 = apply(cbind(d2.join$rand_arm1, d2.join$rand_arm1_new),
                    1,max,na.rm=T) %>% ifelse(.==-Inf, NA, .),
  rand_arm2 = apply(cbind(d2.join$rand_arm2, d2.join$rand_arm2_new),
                    1,max,na.rm=T) %>% ifelse(.==-Inf, NA, .)
) -> df2


# Save data
cbind(data %>% select(-attr_arm1, -attr_arm2, -rand_arm1, -rand_arm2), df, df2) %>% 
  write.xlsx("data.attr.added.xlsx", showNA = FALSE)






