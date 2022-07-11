Error in file(file, "rt") : cannot open the connection
 
Traceback: 
[1] "2: file(file, \"rt\")"                                      
[2] "1: read.table(paste0(data_dir, \"/\", file), header = TRUE)"

[1] "Error in -data_files <- dir(data_dir, pattern = \"\\\\.txt$\", full.names = TRUE) : \n  could not find function \"-<-\"\n"
Traceback: 
No traceback available  

[1] "Error in UseMethod(\"html_table\") : \n  no applicable method for 'html_table' applied to an object of class \"character\"\n"
Traceback: 
[1] "2: html_table(.)"                                                               
[2] "1: haploqa_cc_html %>% html_nodes(\".col-md-12\") %>% html_attr(\"class\") %>% "
[3] "       html_table()"                                                            

[1] "Error in UseMethod(\"html_table\") : \n  no applicable method for 'html_table' applied to an object of class \"character\"\n"
Traceback: 
[1] "2: html_table(.)"                                                           
[2] "1: haploqa_cc_html %>% html_nodes(\".table\") %>% html_attr(\"class\") %>% "
[3] "       html_table()"                                                        

[1] "Error in as.data.frame.data.frame(sum_table, row.names = F) : \n  invalid 'row.names', length 1 for a data frame with 66 rows\n"
Traceback: 
[1] "3: stop(sprintf(ngettext(nr, \"invalid 'row.names', length %d for a data frame with %d row\", "
[2] "       \"invalid 'row.names', length %d for a data frame with %d rows\"), "                    
[3] "       length(row.names), nr), domain = NA)"                                                   
[4] "2: as.data.frame.data.frame(sum_table, row.names = F)"                                         
[5] "1: as.data.frame(sum_table, row.names = F)"                                                    

Error in paste0(data.dir, "/collaborative_cross_summary.csv") : 
  object 'data.dir' not found
 
Traceback: 
[1] "6: paste0(data.dir, \"/collaborative_cross_summary.csv\")"                                
[2] "5: utils::write.table(sum_table, paste0(data.dir, \"/collaborative_cross_summary.csv\"), "
[3] "       row.names = FALSE, col.names = TRUE, sep = \",\", dec = \".\", "                   
[4] "       qmethod = \"double\")"                                                             
[5] "4: eval(expr, p)"                                                                         
[6] "3: eval(expr, p)"                                                                         
[7] "2: eval.parent(Call)"                                                                     
[8] "1: write.csv(sum_table, paste0(data.dir, \"/collaborative_cross_summary.csv\"), "         
[9] "       row.names = FALSE)"                                                                

Error in paste0(data.dir, "/collaborative_cross_summary.csv") : 
  object 'data.dir' not found
 
Traceback: 
1: paste0(data.dir, "/collaborative_cross_summary.csv") 

[1] "Error in fread(paste0(data_dir, \"collaborative_cross_summary.csv\")) : \n  File '/Users/linc/Documents/GitHub/HaploQA_qtl2_comp/haploqa_collaborative_crosscollaborative_cross_summary.csv' does not exist or is non-readable. getwd()=='/Users/linc'\n"
Traceback: 
[1] "2: stop(\"File '\", file, \"' does not exist or is non-readable. getwd()=='\", "
[2] "       getwd(), \"'\")"                                                         
[3] "1: fread(paste0(data_dir, \"collaborative_cross_summary.csv\"))"                

Error in `[.data.frame`(df, df$snp_id == "UNC318619") : 
  undefined columns selected
 
Traceback: 
[1] "3: stop(\"undefined columns selected\")"          
[2] "2: `[.data.frame`(df, df$snp_id == \"UNC318619\")"
[3] "1: df[df$snp_id == \"UNC318619\"]"                

Error: unexpected '=' in "test <- df[df$snp_id ="
 
Traceback: 
No traceback available  

