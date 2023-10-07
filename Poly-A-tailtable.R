rm(list=ls())

library(gt)

name_ <- c(
  "010306.1",

  "015668.1",

 "015874.1",

  "038296.1",

  "040361.1",

  "040534.1",

  "040711.1",

  "055538.1",

  "076908.1"

)

# Create a dataframe with 3 columns and 3 rows
name_df <- as.data.frame(matrix(name_, ncol = 3, byrow = TRUE))
colnames(name_df) <- c("Spalte 1", "Spalte 2", "Spalte 3")

# Create the table using gt
table <- name_df %>%
  gt() %>%
  tab_header(title = "Viren mit Poly-A-Schwanz",
             subtitle = md("**NC-Nummern** der Viren mit Poly-A-Schwanz")) %>%
  opt_align_table_header(align = "left")%>%
  tab_options(column_labels.hidden = TRUE)%>%
  opt_table_font(
    font= google_font("Monderrat"),
    weight=600,style = "italic")

print(table)
