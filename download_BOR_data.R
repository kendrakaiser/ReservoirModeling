library(XML)
setwd("~/Documents/GitRepos/ReservoirModeling/")

# pu: cumulative water year precip; qv: power discharge; px: observed daily totatl precip; fb: reservoir water surface elvation; af: reservoir water storage, acre-feet; ID: computed reservoir inflow

# Assign URL - Arrowrock
url_ark <- "Data/ark_daily.htm" #(10/01/1997 through 09/12/2018)
# Read the HTML table
data_ark <- readHTMLTable(url_ark, header=TRUE, as.data.frame = TRUE, stringsAsFactors=FALSE)
#convert to data frame
ark=data_ark[[1]]


# Assign URL - Anderson Ranch
url_and <- "Data/and_daily.htm"
data_and <- readHTMLTable(url_and, header=TRUE, as.data.frame = TRUE, stringsAsFactors=FALSE)
and=data_and[[1]]

#Lucky Peak
url_luc <- "Data/luc_daily.htm"
data_luc <- readHTMLTable(url_luc, header=TRUE, as.data.frame = TRUE, stringsAsFactors=FALSE)
luc=data_luc[[1]]


res<-data.frame("Date"= as.Date(and[,1]), "andAF"=as.numeric(and$and_af), "arkAF" = as.numeric(ark$ark_af), "lucAF" = as.numeric(luc$luc_af), "in_computed"=as.numeric(luc$luc_id), "in_unreg"= as.numeric(luc$luc_qu), "qo" = as.numeric(luc$luc_qd))

write.csv(res, file="BRB_reservoir_data_1997-2018.csv")

lowell<-"Data/lowell_daily.htm"
data_lowell <- readHTMLTable(lowell, header=TRUE, as.data.frame = TRUE, stringsAsFactors=FALSE)
low=data_lowell[[1]]
write.csv(low, file="Data/Lowell_data.csv")

nyc<-"Data/NY_Canal.htm"
data_nyc <- readHTMLTable(nyc, header=TRUE, as.data.frame = TRUE, stringsAsFactors=FALSE)
ny_canal=data_nyc[[1]]
write.csv(ny_canal, file="Data/NY_canal_data.csv")
