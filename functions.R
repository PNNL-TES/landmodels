
#*****************************************************************************************************************
# prepare data 
#*****************************************************************************************************************
# read file
read_file <- function(x) read.csv(file.path(DATA_DIR, x), comment.char = "#", stringsAsFactors = FALSE)

# function to get Rh data from wieder according to latitude, longitude, year, and day