#-----------------------------------------------------------------------------------------------------#
#		Block 00		GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#

# Copyright statement comment:
#   All rights reserved.
#
# Author comment:
#   Rick Reijnders
#   Script version: 23-09-2019

# File description:
#	Name
#	  Query.R
#
#   Purpose
#     
#   Inputs
#	  


#-----------------------------------------------------------------------------------------------------#
#		Block 00		TODO
#-----------------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------------------#
#		Block 01		(Install &) Load packages
#-----------------------------------------------------------------------------------------------------#
# Install packages if needed and load, or just load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = F)

# Insert all packages in requiredpackages
requiredpackages <-
  c("WikidataQueryServiceR","ggplot2",
	"crayon")
	
for (i in requiredpackages) {
	if (!requireNamespace(i, quietly = TRUE))
		BiocManager::install(i, ask = F, dependencies = TRUE)
	require(as.character(i), character.only = TRUE)
	print(i)
}

# Set message styles
Confirm.style <- green 
Warning.style <- yellow 
Error.style <- red


cat(Warning.style("Done running block 01\n"))
cat(Confirm.style("Starting...\n"))

#-----------------------------------------------------------------------------------------------------#
#		Block 02		GENERAL Settings
#-----------------------------------------------------------------------------------------------------#
# Insert settings

options(stringsAsFactors 	= F)
Verbose 					= 3		#0-3:	0= no feedback  1= prints results  2=prints results + feedback  3= prints all

if(Verbose>=3) cat(Warning.style("Done running block 02a\n"))


#-----------------------------------------------------------------------------------------------------#
#		Block 03		
#-----------------------------------------------------------------------------------------------------#
endpoint = "https://query.wikidata.org/bigdata/namespace/wdq/sparql"
query = 'SELECT ?item WHERE {
?item wdt:P31 wd:Q2934.
}'

results = query_wikidata(query)


SELECT ?anes ?anesLabel ?CC ?bp
WHERE {
?anes wdt:P31/wdt:P279* wd:Q41581 .
?anes wdt:P233  ?CC .
?anes wdt:P2102 ?bp .
SERVICE wikibase:label { bd:serviceParam wikibase:language 'en' }
}




