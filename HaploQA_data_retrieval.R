library(rvest)
library(httr)
library(dplyr)
library(tidyr)
library(rstudioapi)
library(logr)

### set up environment
# get the directory this file is stored in
root <- dirname(getSourceEditorContext()$path)
# create data subdirectory, if not exists
data_dir <- file.path(root, 'haploqa_miniMUGA')
dir.create(data_dir, showWarnings = FALSE)

### log-in to HaploQA
# domain
url_domain <- 'https://haploqa.jax.org' # domain

# username 
username <- 'Christine.Lin@jax.org'
password <- 'christine.lin'

# start a session
url_session <- session(url_domain)
# fetch login form
login_form <- html_form(url_session)[[2]]
# fill credentials
filled_form <- html_form_set(login_form, 'email'=username, 'password'=password)
filled_form$action <- url_domain
# submit to session
login_session <- session_submit(url_session, filled_form)
# get html with submitted/logged in session
session_jump_to(url_session, url_domain) %>% html_nodes(".table") %>% html_table()



### params
# URL to HaploQA sample page
haploqa_cc_html <- read_html('https://haploqa.jax.org/tag/MiniMUGA.html')

# get the main html table first
html_table <- haploqa_cc_html %>% html_nodes(".table") %>% html_table()
sum_table <- html_table[[1]][-1]
sum_table$`ID (Secondary IDs)` <- gsub("\\s+", "", sum_table$`ID (Secondary IDs)`)
sum_table <- as.data.frame(sum_table, row.names = F)



url <-"https://www.silversanz.com/es/account/login"

#create a web session with the desired login address
pgsession<-html_session(url)
pgform<-html_form(pgsession)
